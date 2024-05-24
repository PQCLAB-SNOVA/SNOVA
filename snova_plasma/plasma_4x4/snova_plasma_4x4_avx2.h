/*
Use uint64_t to store a GF16 4x4 matrix (with each half-byte storing one element of GF16)."
Optimizing GF16 4x4 matrix multiplication using AVX2.
*/

#ifndef PLASMA_4x4_AVX2_H
#define PLASMA_4x4_AVX2_H

#include <immintrin.h>
#include <stdint.h>

#include "../../deriv_params.h"
#include "../../gf16.h"
#include "../../snova_kernel.h"
#include "../snova_plasma_avx2.h"
// uint64_t to store a GF16 4x4 matrix
#include "snova_plasma_4x4_data.h"

static __m256i a_mask[4] __attribute__((aligned(32))) = {0};

int init_4x4() {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                set_gf16m_u64(&(S_u64[i]), j, k, get_gf16m(S[i], j, k));
            }
        }
    }

    init_avx_table();

    // 4x4 array mask, and gf16m_u64_mul_4x4_4way_b4
    a_mask[0] = _mm256_set1_epi64x(0x1000010000100001ull);
    a_mask[1] = _mm256_set1_epi64x(0x2000020000200002ull);
    a_mask[2] = _mm256_set1_epi64x(0x4000040000400004ull);
    a_mask[3] = _mm256_set1_epi64x(0x8000080000800008ull);

    return 1;
}

/**
 * a11 x x x  mul  b11 b12 b13 b14  =>  (a11 x b11) (a11 x b12) (a11 x b13) (a11 x b14)
 * x a22 x x  mul  b21 b22 b23 b24  =>  (a22 x b21) (a22 x b22) (a22 x b23) (a22 x b24)   =>  add(xor) to sum
 * x x a33 x  mul  b31 b32 b33 b34  =>  (a33 x b31) (a33 x b32) (a33 x b33) (a33 x b34)
 * x x x a44  mul  b41 b42 b43 b44  =>  (a44 x b41) (a44 x b42) (a44 x b43) (a44 x b44)
 */
static inline void diag_1x1_mul_1x4_4way_add(__m256i at_256, __m256i* b4_256, __m256i* sum) {
    __m256i mask_256;
    for (int j = 0; j < 4; j++) {
        mask_256 = _mm256_cmpeq_epi16(_mm256_and_si256(at_256, a_mask[j]), a_mask[j]);  // a_mask is init from "init_4x4()""
        *sum = _mm256_xor_si256(*sum, _mm256_and_si256(b4_256[j], mask_256));
    }
}

/**
 * a11 a12 a13 a14      a21 a22 a23 a24
 * a21 a22 a23 a24  =>  a31 a32 a33 a34
 * a31 a32 a33 a34      a41 a42 a43 a44
 * a41 a42 a43 a44      a11 a12 a13 a14
 */
static inline __m256i row_cyclic_shift_4x4_4way(__m256i a) {
    // 0123 -> 3012 -> 2301 -> 1230 -> 0123
    __m256i shf_temp_mask = _mm256_setr_epi8(6, 7, 0, 1, 2, 3, 4, 5,         //
                                             14, 15, 8, 9, 10, 11, 12, 13,   //
                                             6, 7, 0, 1, 2, 3, 4, 5,         //
                                             14, 15, 8, 9, 10, 11, 12, 13);  //
    return _mm256_shuffle_epi8(a, shf_temp_mask);
}

/**
 * The difference from "gf16m_u64_mul_4x4_4way" is that "bt4_256" precomputes.
 * ps. bt4_256[4] = { b * t^0, b * t^1, b * t^2, b * t^3 }
 */
static inline void gf16m_u64_mul_4x4_4way_b4(__m256i at_256, __m256i* b4_256, __m256i* c_256) {
    __m256i temp_256 = _mm256_setzero_si256();

    diag_1x1_mul_1x4_4way_add(at_256, b4_256, &temp_256);
    temp_256 = row_cyclic_shift_4x4_4way(temp_256);
    for (int i = 0; i < 3; i++) {
        at_256 = row_cyclic_shift_4x4_4way(at_256);
        diag_1x1_mul_1x4_4way_add(at_256, b4_256, &temp_256);
        temp_256 = row_cyclic_shift_4x4_4way(temp_256);
    }

    *c_256 = temp_256;
}

/**
 * A[0:3] x B[0:3] = C[0:3]
 */
static inline void gf16m_u64_mul_4x4_4way(__m256i at_256, __m256i b_256, __m256i* c_256) {
    __m256i bt4_256[4];
    bt4_256[0] = b_256;
    bt4_256[1] = mul2_avx_r(b_256);
    bt4_256[2] = mul4_avx_r(b_256);
    bt4_256[3] = mul8_avx_r(b_256);
    gf16m_u64_mul_4x4_4way_b4(at_256, bt4_256, c_256);
}

/**
 * A x B[0:3] = C[0:3]
 */
void gf16m_u64_mul_4x4_4way_14(uint64_t a, uint64_t* b, uint64_t* c) {
    __m256i a_256 = _mm256_set1_epi64x(a);
    __m256i b_256 = _mm256_setr_epi64x(b[0], b[1], b[2], b[3]);
    gf16m_u64_mul_4x4_4way(a_256, b_256, (__m256i*)c);
}
void gf16m_256_mul_4x4_4way_14(uint64_t a, __m256i b_256, __m256i* c_256) {
    __m256i a_256 = _mm256_set1_epi64x(a);
    gf16m_u64_mul_4x4_4way(a_256, b_256, c_256);
}

/**
 * C[0:3] ^= A x B[0:3]
 */
void gf16m_u64_mul_4x4_4way_14_add(uint64_t a, uint64_t* b, __m256i* c_256) {
    __m256i temp_256;
    __m256i a_256 = _mm256_set1_epi64x(a);
    __m256i b_256 = _mm256_setr_epi64x(b[0], b[1], b[2], b[3]);
    gf16m_u64_mul_4x4_4way(a_256, b_256, &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}

/**
 * A[0:3] x B = C[0:3]
 */
void gf16m_u64_mul_4x4_4way_41(uint64_t* a, uint64_t b, uint64_t* c) {
    __m256i a_256 = _mm256_setr_epi64x(a[0], a[1], a[2], a[3]);
    __m256i b_256 = _mm256_set1_epi64x(b);
    gf16m_u64_mul_4x4_4way(a_256, b_256, (__m256i*)c);
}

/**
 * C[0:3] ^= A[0:3] x B
 */
void gf16m_u64_mul_4x4_4way_41_add(uint64_t* a, uint64_t b, __m256i* c_256) {
    __m256i temp_256;
    __m256i a_256 = _mm256_setr_epi64x(a[0], a[1], a[2], a[3]);
    __m256i b_256 = _mm256_set1_epi64x(b);
    gf16m_u64_mul_4x4_4way(a_256, b_256, &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}

/**
 * A[0:3] x B[0:3] = C[0:3]
 */
void gf16m_u64_mul_4x4_4way_44(uint64_t* a, uint64_t* b, uint64_t* c) {
    __m256i a_256 = _mm256_setr_epi64x(a[0], a[1], a[2], a[3]);
    __m256i b_256 = _mm256_setr_epi64x(b[0], b[1], b[2], b[3]);
    gf16m_u64_mul_4x4_4way(a_256, b_256, (__m256i*)c);
}
void gf16m_256_mul_4x4_4way_44(__m256i a_256, __m256i b_256, __m256i* c_256) { gf16m_u64_mul_4x4_4way(a_256, b_256, c_256); }

void gf16m_u64_mul_4x4_4way_14_b4(uint64_t a, __m256i* bt4_256, __m256i* c_256) {
    __m256i at_256 = _mm256_set1_epi64x(a);
    gf16m_u64_mul_4x4_4way_b4(at_256, bt4_256, c_256);
}

/**
 * C[0:3] ^= A x B[0:3]
 */
void gf16m_u64_mul_4x4_4way_14_b4_add(uint64_t a, __m256i* bt4_256, __m256i* c_256) {
    __m256i temp_256;
    __m256i at_256 = _mm256_set1_epi64x(a);
    gf16m_u64_mul_4x4_4way_b4(at_256, bt4_256, &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}

/**
 * C[0:3] ^= A[0:3] x B[0:3]
 */
void gf16m_u64_mul_4x4_4way_44_add(uint64_t* a, uint64_t* b, __m256i* c_256) {
    __m256i a_256 = _mm256_set_epi64x(a[0], a[1], a[2], a[3]);
    __m256i b_256 = _mm256_set_epi64x(b[0], b[1], b[2], b[3]);
    __m256i temp_256;
    gf16m_u64_mul_4x4_4way(a_256, b_256, &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}
void gf16m_256_mul_4x4_4way_44_add(__m256i a_256, __m256i b_256, __m256i* c_256) {
    __m256i temp_256;
    gf16m_u64_mul_4x4_4way(a_256, b_256, &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}

/**
 * Generate private key (F part), use avx2 vtl
 * @param map2 - output: F11 F12 F21
 * @param map1 - input: P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param T12 - input
 */
void gen_F_4x4_vtl(map_group2* map2, map_group1* map1, T12_t T12) {
    uint8_t p11_8[mvl_SNOVA * vl_SNOVA] __attribute__((aligned(32)));
    uint8_t t12_8[mvl_SNOVA * l_SNOVA] __attribute__((aligned(32)));
    uint8_t res8[mvl_SNOVA * o_SNOVA * l_SNOVA] __attribute__((aligned(32))) = {0};

    __m256i* p11_256 = (__m256i*)p11_8;
    __m256i* res256 = (__m256i*)res8;

    memcpy(map2->F11, map1->P11, m_SNOVA * v_SNOVA * v_SNOVA * lsq_SNOVA);
    memcpy(map2->F12, map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * lsq_SNOVA);
    memcpy(map2->F21, map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * lsq_SNOVA);

    __m128i transpose4x4 = _mm_set_epi64x(0x0f0b07030e0a0602ll, 0x0d0905010c080400ll);
    uint32_t l_mask32 = 0x0f0f0f0f;

    // F12
    for (int dk = 0; dk < v_SNOVA; ++dk)
        for (int dj = 0; dj < o_SNOVA; ++dj)
            for (int k1 = 0; k1 < l_SNOVA; k1 += 2)
            {
                uint64_t val = *(uint64_t*)&T12[dk][dj][k1 * l_SNOVA];
                *(uint16_t*)&t12_8[(dk * l_SNOVA + k1) * 2 * o_SNOVA + 2 * dj] = (val & 0x0f0f) ^ ((val >> 12) & 0xf0f0);
                *(uint16_t*)&t12_8[(dk * l_SNOVA + k1 + 1) * 2 * o_SNOVA + 2 * dj] = ((val >> 32) & 0x0f0f) ^ ((val >> 44) & 0xf0f0);
            }

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int di = 0; di < v_SNOVA; di++)
            for (int dk = 0; dk < v_SNOVA; dk++)
            {
                uint32_t* pp11 = (uint32_t*)map1->P11[mi][di][dk];
                __m128i val128 = _mm_setr_epi32(pp11[0], pp11[1], pp11[2], pp11[3]);
                uint32_t tres[4];
                _mm_storeu_si128((__m128i_u*)&tres, _mm_shuffle_epi8(val128, transpose4x4));
                uint32_t* pres = (uint32_t*)&p11_8[dk * l_SNOVA * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + di * l_SNOVA];

                for (int idx = 0; idx < l_SNOVA; idx++)
                    pres[idx * mvl_SNOVA / 4] = tres[idx];
            }

    for (int dk_k1 = 0; dk_k1 < v_SNOVA * l_SNOVA; ++dk_k1)
        for (int dj_j1 = 0; dj_j1 < 2 * o_SNOVA; dj_j1++)
        {
            __m256i k_lh = vtl_ct_multtab(t12_8[dk_k1 * 2 * o_SNOVA + dj_j1] & 0xf) ^
                _mm256_slli_epi16(vtl_ct_multtab(t12_8[dk_k1 * 2 * o_SNOVA + dj_j1] >> 4), 4);
            for (int mi_di_i1 = 0; mi_di_i1 < mvl_SNOVA32; mi_di_i1++)
                res256[dj_j1 * mvl_SNOVA32 + mi_di_i1] ^= _mm256_shuffle_epi8(k_lh, p11_256[dk_k1 * mvl_SNOVA32 + mi_di_i1]);
        }

    for (int di = 0; di < v_SNOVA; ++di)
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int dj = 0; dj < o_SNOVA; ++dj)
            {
                uint32_t* res32 = (uint32_t*)&res8[dj * 2 * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + di * l_SNOVA];
                __m128i val128 = _mm_setr_epi32(res32[0] & l_mask32, res32[mvl_SNOVA / 4] & l_mask32,
                                                (res32[0] >> 4) & l_mask32, (res32[mvl_SNOVA / 4] >> 4) & l_mask32);

                *(__m128i*)map2->F12[mi][di][dj] ^= _mm_shuffle_epi8(val128, transpose4x4);
            }

    // Same for F21
    memset(res256, 0, sizeof(res8));

    for (int dj = 0; dj < v_SNOVA; ++dj)
        for (int dk = 0; dk < o_SNOVA; ++dk)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
            {
                t12_8[(dj * l_SNOVA + j1) * 2 * o_SNOVA + 2 * dk + 0] = T12[dj][dk][0 * l_SNOVA + j1] ^ (T12[dj][dk][2 * l_SNOVA + j1] << 4);
                t12_8[(dj * l_SNOVA + j1) * 2 * o_SNOVA + 2 * dk + 1] = T12[dj][dk][1 * l_SNOVA + j1] ^ (T12[dj][dk][3 * l_SNOVA + j1] << 4);
            }

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int di = 0; di < v_SNOVA; di++)
            for (int dk = 0; dk < v_SNOVA; dk++)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    *(uint32_t*)&p11_8[(di * l_SNOVA + i1) * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + dk * l_SNOVA] =
                        *(uint32_t*)&map1->P11[mi][di][dk][i1 * l_SNOVA];

    for (int dk_k1 = 0; dk_k1 < v_SNOVA * l_SNOVA; ++dk_k1)
        for (int dj_j1 = 0; dj_j1 < 2 * o_SNOVA; dj_j1++)
        {
            __m256i k_lh = vtl_ct_multtab(t12_8[dk_k1 * 2 * o_SNOVA + dj_j1] & 0xf) ^
                _mm256_slli_epi16(vtl_ct_multtab(t12_8[dk_k1 * 2 * o_SNOVA + dj_j1] >> 4), 4);
            for (int mi_di_i1 = 0; mi_di_i1 < mvl_SNOVA32; mi_di_i1++)
                res256[dj_j1 * mvl_SNOVA32 + mi_di_i1] ^= _mm256_shuffle_epi8(k_lh, p11_256[dk_k1 * mvl_SNOVA32 + mi_di_i1]);
        }

    // Shuffle back
    for (int dj = 0; dj < o_SNOVA; ++dj)
        for (int di = 0; di < v_SNOVA; ++di)
            for (int mi = 0; mi < m_SNOVA; ++mi) {
                uint32_t* res32 = (uint32_t*)&res8[dj * 2 * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + di * l_SNOVA];
                __m128i val128 = _mm_setr_epi32(res32[0] & l_mask32, res32[mvl_SNOVA / 4] & l_mask32,
                                                (res32[0] >> 4) & l_mask32, (res32[mvl_SNOVA / 4] >> 4 ) & l_mask32);
                *(__m128i*)map2->F21[mi][dj][di] ^= val128;
            }
}

#endif