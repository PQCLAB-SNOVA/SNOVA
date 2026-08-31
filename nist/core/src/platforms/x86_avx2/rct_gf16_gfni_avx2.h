#ifndef RCT_GF16_GFNI_AVX2_H
#define RCT_GF16_GFNI_AVX2_H

#if RCT_USE_GFNI
#define RCT_BC(s)    _mm256_set1_epi8((char)(s))
#define RCT_SV(bc, v) RCT_GFMUL256((bc), (v))
#define RCT_BC128(s)     _mm_set1_epi8((char)(s))
#define RCT_SV128(bc, v) RCT_GFMUL128((bc), (v))
#define RCT_BC_SEC(s)    RCT_BC(s)
#define RCT_BC128_SEC(s) RCT_BC128(s)

static inline void rct_gauss_row_axpy(gf_t *dst, const gf_t *src, gf_t s, int k0, int kend) {
    __m256i sv = _mm256_set1_epi8((char)s);
    for (int k = k0; k < kend; k += 32) {
        __m256i p = RCT_GFMUL256(sv, _mm256_loadu_si256((const __m256i *)(src + k)));
        __m256i d = _mm256_loadu_si256((const __m256i *)(dst + k));
        _mm256_storeu_si256((__m256i *)(dst + k), _mm256_xor_si256(d, rct_gfni_cleanup256(p)));
    }
}
static inline void rct_gauss_row_scale(gf_t *row, gf_t s, int k0, int kend) {
    __m256i sv = _mm256_set1_epi8((char)s);
    for (int k = k0; k < kend; k += 32) {
        __m256i p = RCT_GFMUL256(sv, _mm256_loadu_si256((const __m256i *)(row + k)));
        _mm256_storeu_si256((__m256i *)(row + k), rct_gfni_cleanup256(p));
    }
}

#if SNOVA_l == 4
static inline void rct_gf4_matmul_add(gf_t *acc, const gf_t *A, const gf_t *B) {
    static const _Alignas(16) uint8_t AMASK[4][16] = {
        {0,0,0,0, 4,4,4,4, 8,8,8,8, 12,12,12,12},
        {1,1,1,1, 5,5,5,5, 9,9,9,9, 13,13,13,13},
        {2,2,2,2, 6,6,6,6, 10,10,10,10, 14,14,14,14},
        {3,3,3,3, 7,7,7,7, 11,11,11,11, 15,15,15,15}};
    static const _Alignas(16) uint8_t BMASK[4][16] = {
        {0,1,2,3, 0,1,2,3, 0,1,2,3, 0,1,2,3},
        {4,5,6,7, 4,5,6,7, 4,5,6,7, 4,5,6,7},
        {8,9,10,11, 8,9,10,11, 8,9,10,11, 8,9,10,11},
        {12,13,14,15, 12,13,14,15, 12,13,14,15, 12,13,14,15}};
    __m128i av = _mm_loadu_si128((const __m128i *)A);
    __m128i bv = _mm_loadu_si128((const __m128i *)B);
    __m128i prod = _mm_setzero_si128();
    for (int k = 0; k < 4; ++k) {
        __m128i ak = _mm_shuffle_epi8(av, _mm_load_si128((const __m128i *)AMASK[k]));
        __m128i bk = _mm_shuffle_epi8(bv, _mm_load_si128((const __m128i *)BMASK[k]));
        prod = _mm_xor_si128(prod, RCT_GFMUL128(ak, bk));
    }
    __m128i r = rct_gfni_cleanup128(prod);
    __m128i cur = _mm_loadu_si128((const __m128i *)acc);
    _mm_storeu_si128((__m128i *)acc, _mm_xor_si128(cur, r));
}

static inline void rct_matmul_l4rows(gf_t *C, const gf_t *A, const gf_t *B, int ad, int bd) {
    for (int i = 0; i < ad; i++) {
        __m128i acc = _mm_setzero_si128();
        for (int k = 0; k < bd; k++)
            acc = _mm_xor_si128(acc, RCT_GFMUL128(
                _mm_set1_epi8((char)A[i * bd + k]),
                _mm_loadu_si128((const __m128i *)&B[k * 4])));
        int32_t v = _mm_cvtsi128_si32(rct_gfni_cleanup128(acc));
        memcpy(&C[i * 4], &v, 4);
    }
}

static inline void rct_matmul_l4rows_add(gf_t *C, const gf_t *A, const gf_t *B, int ad, int bd) {
    for (int i = 0; i < ad; i++) {
        __m128i acc = _mm_setzero_si128();
        for (int k = 0; k < bd; k++)
            acc = _mm_xor_si128(acc, RCT_GFMUL128(
                _mm_set1_epi8((char)A[i * bd + k]),
                _mm_loadu_si128((const __m128i *)&B[k * 4])));
        int32_t v = _mm_cvtsi128_si32(rct_gfni_cleanup128(acc));
        int32_t c;
        memcpy(&c, &C[i * 4], 4);
        c ^= v;
        memcpy(&C[i * 4], &c, 4);
    }
}
#endif

#define rct_gauss_row_scale_sec rct_gauss_row_scale
#define rct_gauss_row_axpy_sec  rct_gauss_row_axpy
#if SNOVA_l == 4
#define rct_gf4_matmul_add_sec  rct_gf4_matmul_add
#define rct_matmul_l4rows_sec   rct_matmul_l4rows
#endif
#endif

#endif
