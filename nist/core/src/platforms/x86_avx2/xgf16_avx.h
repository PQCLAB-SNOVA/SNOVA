/**
 * @file xgf16_avx.h
 */
#ifndef M2_XGF16_AVX_H
#define M2_XGF16_AVX_H

#include <stdint.h>
#include <immintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

extern __m256i vtl_mt[16];

extern __m256i mtk2_16[256];

void vtl_init(void);

extern __m256i vtl_basis_t1, vtl_basis_t2, vtl_basis_t4, vtl_basis_t8;
extern __m256i vtl_basis_m1, vtl_basis_m2, vtl_basis_m4, vtl_basis_m8;
extern __m256i vtl_basis_zero;

static inline __m256i vtl_ct_multtab(uint8_t k) {
    __m256i v = _mm256_set1_epi32(k & 0x0F);
    __m256i r = _mm256_setzero_si256();
    r = _mm256_xor_si256(r, _mm256_and_si256(vtl_basis_t1,
            _mm256_cmpgt_epi32(_mm256_and_si256(v, vtl_basis_m1), vtl_basis_zero)));
    r = _mm256_xor_si256(r, _mm256_and_si256(vtl_basis_t2,
            _mm256_cmpgt_epi32(_mm256_and_si256(v, vtl_basis_m2), vtl_basis_zero)));
    r = _mm256_xor_si256(r, _mm256_and_si256(vtl_basis_t4,
            _mm256_cmpgt_epi32(_mm256_and_si256(v, vtl_basis_m4), vtl_basis_zero)));
    r = _mm256_xor_si256(r, _mm256_and_si256(vtl_basis_t8,
            _mm256_cmpgt_epi32(_mm256_and_si256(v, vtl_basis_m8), vtl_basis_zero)));
    return r;
}

static inline __m256i vtl_ct_multtab_pair(uint8_t k_lo, uint8_t k_hi) {
    __m256i t_lo = vtl_ct_multtab(k_lo);
    __m256i t_hi = vtl_ct_multtab(k_hi);
    return _mm256_or_si256(t_lo, _mm256_slli_epi16(t_hi, 4));
}

/**
 * @return c.
 */
static inline __m256i gf16_32_mul_k(__m256i a, uint8_t k) {
    return _mm256_shuffle_epi8(vtl_mt[k & 0x0F], a);
}

/**
 * @param a 32 GF(16) elements. @param k Scalar. @param acc Accumulator.
 */
static inline __m256i gf16_32_mul_k_add(__m256i a, uint8_t k, __m256i acc) {
    return _mm256_xor_si256(acc, _mm256_shuffle_epi8(vtl_mt[k & 0x0F], a));
}

/**
 * @param a 32 GF(16) elements. @param b 32 GF(16) elements.
 * @return c[i] = a[i] * b[i].
 */
__m256i gf16_32_mul_32(__m256i a, __m256i b);

/**
 * @param a 32 GF(16) elements. @param b 32 GF(16) elements. @param acc Accumulator.
 */
__m256i gf16_32_mul_32_add(__m256i a, __m256i b, __m256i acc);

#ifdef __cplusplus
}
#endif

#endif
