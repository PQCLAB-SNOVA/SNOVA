/*
Use uint16_t to store a GF16 2x2 matrix (with each half-byte storing one element of GF16)."
Optimizing GF16 2x2 matrix multiplication using AVX2.
*/

#ifndef PLASMA_2x2_AVX2_H
#define PLASMA_2x2_AVX2_H

#include <immintrin.h>
#include <stdint.h>

#include "../../deriv_params.h"
#include "../../gf16.h"
#include "../../snova_kernel.h"
#include "../snova_plasma_avx2.h"
// uint16_t to store a GF16 2x2 matrix
#include "snova_plasma_2x2_data.h"

static __m256i a_mask[4] __attribute__((aligned(32))) = {0};

void init_2x2() {
    for (int i = 0; i < rank; i++) {
        for (int j = 0; j < rank; j++) {
            for (int k = 0; k < rank; k++) {
                set_gf16m_u16(&(S_u16[i]), j, k, get_gf16m(S[i], j, k));
            }
        }
    }

    init_avx_table();

    // 4x4 array mask, used gf16m_u64_mul_4x4_4way and gf16m_u64_mul_4x4_4way_b4
    // bit1
    a_mask[0] = _mm256_set1_epi64x(0x1001100110011001ull);
    a_mask[1] = _mm256_set1_epi64x(0x2002200220022002ull);
    a_mask[2] = _mm256_set1_epi64x(0x4004400440044004ull);
    a_mask[3] = _mm256_set1_epi64x(0x8008800880088008ull);
}

/**
 * a11 x  mul  b11 b12 =>  (a11 x b11) (a11 x b12)
 * x a22  mul  b21 b22 =>  (a22 x b21) (a22 x b22)   =>  add(xor) to sum
 */
static inline void diag_1x1_mul_1x2_4way_add(__m256i at_256, __m256i* b4_256, __m256i* sum) {
    __m256i mask_256;
    for (int j = 0; j < 4; j++) {
        mask_256 = _mm256_cmpeq_epi8(_mm256_and_si256(at_256, a_mask[j]), a_mask[j]);
        *sum = _mm256_xor_si256(*sum, _mm256_and_si256(b4_256[j], mask_256));
    }
}

/**
 * a11 a12      a21 a22
 * a21 a22  =>  a11 a12
 */
static inline __m256i row_cyclic_shift_2x2_16way(__m256i a) {
    // 01 -> 10 -> 01
    __m256i shf_temp_mask = _mm256_setr_epi8(1, 0, 3, 2, 5, 4, 7, 6,         //
                                             9, 8, 11, 10, 13, 12, 15, 14,   //
                                             1, 0, 3, 2, 5, 4, 7, 6,         //
                                             9, 8, 11, 10, 13, 12, 15, 14);  //
    return _mm256_shuffle_epi8(a, shf_temp_mask);
}

static inline void gf16m_u64_mul_2x2_16way_b4(__m256i at_256, __m256i* b4_256, __m256i* c_256) {
    __m256i temp_256 = _mm256_setzero_si256();

    diag_1x1_mul_1x2_4way_add(at_256, b4_256, &temp_256);

    temp_256 = row_cyclic_shift_2x2_16way(temp_256);
    at_256 = row_cyclic_shift_2x2_16way(at_256);
    diag_1x1_mul_1x2_4way_add(at_256, b4_256, &temp_256);

    temp_256 = row_cyclic_shift_2x2_16way(temp_256);

    *c_256 = temp_256;
}

/**
 * return _mm256_setr_epi16(u16[0], u16[0], u16[0], u16[0], u16[1], u16[1], u16[1], u16[1], u16[2], u16[2], u16[2], u16[2],
 * u16[3], u16[3], u16[3], u16[3]);
 */
static inline __m256i u64_2_u256_4_16(uint64_t u64) {
    return _mm256_shuffle_epi8(_mm256_set1_epi64x(u64), _mm256_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3, 2, 3, 2, 3, 4, 5,
                                                                         4, 5, 4, 5, 4, 5, 6, 7, 6, 7, 6, 7, 6, 7));
}

/**
 * A[0:15] x B[0:15] = C[0:15]
 */
static inline void gf16m_u64_mul_2x2_16way(__m256i a_256, __m256i b_256, __m256i* c_256) {
    __m256i bt4_256[4];
    bt4_256[0] = b_256;
    bt4_256[1] = mul2_avx_r(b_256);
    bt4_256[2] = mul4_avx_r(b_256);
    bt4_256[3] = mul8_avx_r(b_256);
    gf16m_u64_mul_2x2_16way_b4(a_256, bt4_256, c_256);
}

void gf16m_u64_mul_2x2_16way_4_16(uint64_t* a_256, __m256i* b_256, __m256i* c_256) {
    gf16m_u64_mul_2x2_16way(u64_2_u256_4_16(*a_256), *b_256, c_256);
}

void gf16m_u64_mul_2x2_16way_16_4(__m256i* a_256, uint64_t* b_256, __m256i* c_256) {
    gf16m_u64_mul_2x2_16way(*a_256, u64_2_u256_4_16(*b_256), c_256);
}

void gf16m_u64_mul_2x2_16way_16_16(__m256i* a_256, __m256i* b_256, __m256i* c_256) {
    gf16m_u64_mul_2x2_16way(*a_256, *b_256, c_256);
}

void gf16m_u64_mul_2x2_16way_4_16_b4_add(uint64_t* a_256, __m256i* bt4_256, __m256i* c_256) {
    __m256i temp_256;
    gf16m_u64_mul_2x2_16way_b4(u64_2_u256_4_16(*a_256), bt4_256, &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}

void gf16m_u64_mul_2x2_16way_16_4_add(__m256i* a_256, uint64_t* b_256, __m256i* c_256) {
    __m256i temp_256;
    gf16m_u64_mul_2x2_16way(*a_256, u64_2_u256_4_16(*b_256), &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}

void gf16m_u64_mul_2x2_16way_4_16_add(uint64_t* a_256, __m256i* b_256, __m256i* c_256) {
    __m256i temp_256;
    gf16m_u64_mul_2x2_16way(u64_2_u256_4_16(*a_256), *b_256, &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}

void gf16m_u64_mul_2x2_16way_16_16_add(__m256i* a_256, __m256i* b_256, __m256i* c_256) {
    __m256i temp_256;
    gf16m_u64_mul_2x2_16way(*a_256, *b_256, &temp_256);
    *c_256 = _mm256_xor_si256(*c_256, temp_256);
}

static inline uint16_t xor_u16_elements(__m256i data) {
    __m256i temp = _mm256_xor_si256(data, _mm256_srli_epi64(data, 32));
    temp = _mm256_xor_si256(temp, _mm256_srli_epi32(temp, 16));
    return _mm256_extract_epi16(temp, 0) ^ _mm256_extract_epi16(temp, 4) ^ _mm256_extract_epi16(temp, 8) ^
           _mm256_extract_epi16(temp, 12);
}

void gen_F_2x2(map_group2* map2, map_group1* map1, T12_t T12) {
    // memcpy(map2->F11, map1->P11, m_SNOVA * v_SNOVA * v_SNOVA * sq_rank);
    // memcpy(map2->F12, map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * sq_rank);
    // memcpy(map2->F21, map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * sq_rank);

    uint16_t T12_u16[v_SNOVA][o_SNOVA] __attribute__((aligned(32))) = {0};
    uint16_t T12_u16_ov[o_SNOVA][v_SNOVA_mult16] __attribute__((aligned(32))) = {0};
    uint16_t P11_u16[m_SNOVA][v_SNOVA][v_SNOVA_mult16] __attribute__((aligned(32))) = {0};
    uint16_t P11_u16_vr[m_SNOVA][v_SNOVA][v_SNOVA_mult16] __attribute__((aligned(32))) = {0};
    uint16_t P12_u16[m_SNOVA][v_SNOVA][o_SNOVA] __attribute__((aligned(32)));
    uint16_t P21_u16[m_SNOVA][o_SNOVA][v_SNOVA] __attribute__((aligned(32)));

    uint16_t F11_u16[m_SNOVA][v_SNOVA][v_SNOVA] __attribute__((aligned(32)));
    uint16_t F12_u16[m_SNOVA][v_SNOVA][o_SNOVA] __attribute__((aligned(32)));
    uint16_t F21_u16[m_SNOVA][o_SNOVA][v_SNOVA] __attribute__((aligned(32)));

    convert_GF16s_to_bytes((uint8_t*)T12_u16, (uint8_t*)T12, v_SNOVA * o_SNOVA * sq_rank);
    for (int i = 0; i < m_SNOVA; i++) {
        for (int j = 0; j < v_SNOVA; ++j) {
            convert_GF16s_to_bytes((uint8_t*)(P11_u16[i][j]), (uint8_t*)(map1->P11[i][j]), v_SNOVA * sq_rank);
            memcpy((uint8_t*)(F11_u16[i][j]), (uint8_t*)(P11_u16[i][j]), v_SNOVA * sq_rank / 2);
        }
    }

    convert_GF16s_to_bytes((uint8_t*)P12_u16, (uint8_t*)map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * sq_rank);
    convert_GF16s_to_bytes((uint8_t*)P21_u16, (uint8_t*)map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * sq_rank);

    memcpy((uint8_t*)F12_u16, (uint8_t*)P12_u16, m_SNOVA * v_SNOVA * o_SNOVA * sq_rank / 2);
    memcpy((uint8_t*)F21_u16, (uint8_t*)P21_u16, m_SNOVA * o_SNOVA * v_SNOVA * sq_rank / 2);

    for (int i = 0; i < o_SNOVA; ++i) {
        for (int j = 0; j < v_SNOVA; ++j) {
            T12_u16_ov[i][j] = T12_u16[j][i];
        }
    }

    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < v_SNOVA; ++j) {
            for (int k = 0; k < v_SNOVA; ++k) {
                P11_u16_vr[i][j][k] = P11_u16[i][k][j];
            }
        }
    }

    // 6960 4x4 mul
    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < v_SNOVA; ++j) {
            for (int k = 0; k < o_SNOVA; ++k) {
                __m256i sum = _mm256_setzero_si256();
                for (int index = 0; index < v_SNOVA; index += 16) {
                    gf16m_u64_mul_2x2_16way_16_16_add((__m256i*)(P11_u16[i][j] + index), (__m256i*)(T12_u16_ov[k] + index),
                                                      &sum);
                }
                F12_u16[i][j][k] ^= xor_u16_elements(sum);
            }
        }

        for (int j = 0; j < o_SNOVA; ++j) {
            for (int k = 0; k < v_SNOVA; ++k) {
                __m256i sum = _mm256_setzero_si256();
                for (int index = 0; index < v_SNOVA; index += 16) {
                    gf16m_u64_mul_2x2_16way_16_16_add((__m256i*)(T12_u16_ov[j] + index), (__m256i*)(P11_u16_vr[i][k] + index),
                                                      &sum);
                }
                F21_u16[i][j][k] ^= xor_u16_elements(sum);
            }
        }
    }

    convert_bytes_to_GF16s((uint8_t*)F11_u16, (uint8_t*)map2->F11, m_SNOVA * v_SNOVA * v_SNOVA * sq_rank);
    convert_bytes_to_GF16s((uint8_t*)F12_u16, (uint8_t*)map2->F12, m_SNOVA * v_SNOVA * o_SNOVA * sq_rank);
    convert_bytes_to_GF16s((uint8_t*)F21_u16, (uint8_t*)map2->F21, m_SNOVA * o_SNOVA * v_SNOVA * sq_rank);
}

/**
 * Generate private key (F part)
 */
void gen_F_2x2_vtl(map_group2 *map2, map_group1 *map1, T12_t T12)
{
    __m256i p11_256[mvl_SNOVA32 * vl_SNOVA] = {0};
    __m256i t12_256[mvl_SNOVA32 * l_SNOVA] = {0};
    __m256i res256[mvl_SNOVA32 * o_SNOVA] = {0};

    uint8_t *p11_8 = (uint8_t *)p11_256;
    uint8_t *t12_8 = (uint8_t *)t12_256;
    uint8_t *res8 = (uint8_t *)res256;

    memcpy(map2->F11, map1->P11, m_SNOVA * v_SNOVA * v_SNOVA * lsq_SNOVA);
    memcpy(map2->F12, map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * lsq_SNOVA);
    memcpy(map2->F21, map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * lsq_SNOVA);

    for (int k1 = 0; k1 < l_SNOVA; ++k1)
        for (int dj = 0; dj < o_SNOVA; ++dj)
            for (int dk = 0; dk < v_SNOVA; ++dk)
                t12_8[(dk * l_SNOVA + k1) * o_SNOVA + dj] = T12[dk][dj][k1 * l_SNOVA + 0] ^ (T12[dk][dj][k1 * l_SNOVA + 1] << 4);

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int dk = 0; dk < v_SNOVA; dk++)
            for (int di = 0; di < v_SNOVA; di++)
            {
                uint32_t val = *(uint32_t *)map1->P11[mi][di][dk];
                *(uint16_t *)&p11_8[(dk * l_SNOVA + 0) * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + di * l_SNOVA] = ((val >> 0) & 0xf) ^ ((val >> 8) & 0xf00);
                *(uint16_t *)&p11_8[(dk * l_SNOVA + 1) * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + di * l_SNOVA] = ((val >> 8) & 0xf) ^ ((val >> 16) & 0xf00);
            }

    for (int dk_k1 = 0; dk_k1 < v_SNOVA * l_SNOVA; ++dk_k1)
        for (int dj_j1 = 0; dj_j1 < o_SNOVA; dj_j1++)
        {
            __m256i k_lh = vtl_ct_multtab(t12_8[dk_k1 * o_SNOVA + dj_j1] & 0xf) ^
                _mm256_slli_epi16(vtl_ct_multtab(t12_8[dk_k1 * o_SNOVA + dj_j1] >> 4), 4);
            for (int mi_di_i1 = 0; mi_di_i1 < mvl_SNOVA32; mi_di_i1++)
                res256[dj_j1 * mvl_SNOVA32 + mi_di_i1] ^= _mm256_shuffle_epi8(k_lh, p11_256[dk_k1 * mvl_SNOVA32 + mi_di_i1]);
        }

    for (int di = 0; di < v_SNOVA; ++di)
        for (int dj = 0; dj < o_SNOVA; ++dj)
            for (int mi = 0; mi < m_SNOVA; ++mi)
                {
                    uint32_t val = *(uint32_t *)&res8[dj * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + di * l_SNOVA];
                    *(uint32_t *)map2->F12[mi][di][dj] ^= 
                        (val & 0xf) ^ ((val << 4) & 0xf00) ^ ((val << 8) & 0xf0000) ^ ((val << 12) & 0xf000000);
                }

    // Same for F21
    memset(res256, 0, sizeof(res256));

    for (int dj = 0; dj < v_SNOVA; ++dj)
        for (int dk = 0; dk < o_SNOVA; ++dk)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                t12_8[(dj * l_SNOVA + j1) * o_SNOVA + dk] = T12[dj][dk][j1] ^ (T12[dj][dk][l_SNOVA + j1] << 4);

    for (int di = 0; di < v_SNOVA; di++)
        for (int mi = 0; mi < m_SNOVA; mi++)
            for (int dk = 0; dk < v_SNOVA; dk++)
            {
                uint32_t val = *(uint32_t *)map1->P11[mi][di][dk];
                *(uint16_t *)&p11_8[(di * l_SNOVA + 0) * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + dk * l_SNOVA] = val & 0x0f0f;
                *(uint16_t *)&p11_8[(di * l_SNOVA + 1) * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + dk * l_SNOVA] = (val >> 16) & 0x0f0f;
            }

    for (int dk_k1 = 0; dk_k1 < v_SNOVA * l_SNOVA; ++dk_k1)
        for (int dj_j1 = 0; dj_j1 < o_SNOVA; dj_j1++)
        {
            __m256i k_lh = vtl_ct_multtab(t12_8[dk_k1 * o_SNOVA + dj_j1] & 0xf) ^
                _mm256_slli_epi16(vtl_ct_multtab(t12_8[dk_k1 * o_SNOVA + dj_j1] >> 4), 4);
            for (int mi_di_i1 = 0; mi_di_i1 < mvl_SNOVA32; mi_di_i1++)
                res256[dj_j1 * mvl_SNOVA32 + mi_di_i1] ^= _mm256_shuffle_epi8(k_lh, p11_256[dk_k1 * mvl_SNOVA32 + mi_di_i1]);
        }

    // Shuffle back
    for (int mi = 0; mi < m_SNOVA; ++mi)
        for (int dj = 0; dj < o_SNOVA; ++dj)
            for (int di = 0; di < v_SNOVA; ++di)
            {
                uint32_t val = *(uint16_t *)&res8[dj * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + di * l_SNOVA];
                *(uint32_t *)map2->F21[mi][dj][di] ^= (val ^ (val << 12)) & 0x0f0f0f0f;
            }
}

#endif