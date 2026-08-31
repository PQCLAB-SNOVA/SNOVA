#ifndef RCT_GF16_PSHUFB_AVX2_H
#define RCT_GF16_PSHUFB_AVX2_H

#if RCT_USE_PSHUFB
static __m256i rct_vtl[16];
static void rct_build_vtl(void) {
    _Alignas(32) uint8_t buf[32];
    for (int k = 0; k < 16; ++k) {
        for (int j = 0; j < 16; ++j) { buf[j] = rct_multtab[k * SNOVA_q + j]; buf[16 + j] = buf[j]; }
        rct_vtl[k] = _mm256_loadu_si256((const __m256i *)buf);
    }
}
static inline __m128i rct_vtl128(uint8_t k) { return _mm256_castsi256_si128(rct_vtl[k & 0x0F]); }

#define RCT_BC(s)    rct_vtl[(uint8_t)(s) & 0x0F]
#define RCT_SV(bc, v) _mm256_shuffle_epi8((bc), (v))
#define RCT_BC128(s)     rct_vtl128((uint8_t)(s))
#define RCT_SV128(bc, v) _mm_shuffle_epi8((bc), (v))

static inline __m256i rct_bc_sec(uint8_t s) {
    const __m256i sv = _mm256_set1_epi8((char)(s & 0x0F));
    const __m256i b1 = _mm256_set1_epi8(1), b2 = _mm256_set1_epi8(2),
                  b4 = _mm256_set1_epi8(4), b8 = _mm256_set1_epi8(8);
    __m256i r;
    r = _mm256_and_si256(rct_vtl[1], _mm256_cmpeq_epi8(_mm256_and_si256(sv, b1), b1));
    r = _mm256_xor_si256(r, _mm256_and_si256(rct_vtl[2], _mm256_cmpeq_epi8(_mm256_and_si256(sv, b2), b2)));
    r = _mm256_xor_si256(r, _mm256_and_si256(rct_vtl[4], _mm256_cmpeq_epi8(_mm256_and_si256(sv, b4), b4)));
    r = _mm256_xor_si256(r, _mm256_and_si256(rct_vtl[8], _mm256_cmpeq_epi8(_mm256_and_si256(sv, b8), b8)));
    return r;
}
static inline __m128i rct_bc128_sec(uint8_t s) { return _mm256_castsi256_si128(rct_bc_sec(s)); }
#define RCT_BC_SEC(s)    rct_bc_sec((uint8_t)(s))
#define RCT_BC128_SEC(s) rct_bc128_sec((uint8_t)(s))

static inline __m128i rct_gf16_mul128_sec(__m128i av, __m128i bv) {
    const __m128i b1 = _mm_set1_epi8(1), b2 = _mm_set1_epi8(2),
                  b4 = _mm_set1_epi8(4), b8 = _mm_set1_epi8(8);
    __m128i r;
    r = _mm_and_si128(_mm_shuffle_epi8(rct_vtl128(1), bv), _mm_cmpeq_epi8(_mm_and_si128(av, b1), b1));
    r = _mm_xor_si128(r, _mm_and_si128(_mm_shuffle_epi8(rct_vtl128(2), bv), _mm_cmpeq_epi8(_mm_and_si128(av, b2), b2)));
    r = _mm_xor_si128(r, _mm_and_si128(_mm_shuffle_epi8(rct_vtl128(4), bv), _mm_cmpeq_epi8(_mm_and_si128(av, b4), b4)));
    r = _mm_xor_si128(r, _mm_and_si128(_mm_shuffle_epi8(rct_vtl128(8), bv), _mm_cmpeq_epi8(_mm_and_si128(av, b8), b8)));
    return r;
}

static inline void rct_gauss_row_scale(gf_t *row, gf_t s, int k0, int kend) {
    __m256i t = rct_vtl[s & 0x0F];
    for (int k = k0; k < kend; k += 32) {
        __m256i v = _mm256_loadu_si256((const __m256i *)(row + k));
        _mm256_storeu_si256((__m256i *)(row + k), _mm256_shuffle_epi8(t, v));
    }
}
static inline void rct_gauss_row_axpy(gf_t *dst, const gf_t *src, gf_t s, int k0, int kend) {
    __m256i t = rct_vtl[s & 0x0F];
    for (int k = k0; k < kend; k += 32) {
        __m256i p = _mm256_shuffle_epi8(t, _mm256_loadu_si256((const __m256i *)(src + k)));
        __m256i d = _mm256_loadu_si256((const __m256i *)(dst + k));
        _mm256_storeu_si256((__m256i *)(dst + k), _mm256_xor_si256(d, p));
    }
}

#if SNOVA_l == 4
static inline void rct_gf4_matmul_add(gf_t *acc, const gf_t *A, const gf_t *B) {
    static const _Alignas(16) int8_t BM[4][16] = {
        {0,1,2,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {4,5,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {8,9,10,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}};
    __m128i bv = _mm_loadu_si128((const __m128i *)B);
    __m128i brow[4];
    for (int k = 0; k < 4; ++k) brow[k] = _mm_shuffle_epi8(bv, _mm_load_si128((const __m128i *)BM[k]));
    for (int i = 0; i < 4; ++i) {
        __m128i a = _mm_setzero_si128();
        for (int k = 0; k < 4; ++k)
            a = _mm_xor_si128(a, _mm_shuffle_epi8(rct_vtl128(A[i * 4 + k]), brow[k]));
        uint32_t cur; memcpy(&cur, &acc[i * 4], 4);
        cur ^= (uint32_t)_mm_cvtsi128_si32(a);
        memcpy(&acc[i * 4], &cur, 4);
    }
}
static inline void rct_matmul_l4rows(gf_t *C, const gf_t *A, const gf_t *B, int ad, int bd) {
    for (int i = 0; i < ad; ++i) {
        __m128i acc = _mm_setzero_si128();
        for (int k = 0; k < bd; ++k)
            acc = _mm_xor_si128(acc, _mm_shuffle_epi8(rct_vtl128(A[i * bd + k]),
                _mm_loadu_si128((const __m128i *)&B[k * 4])));
        uint32_t r = (uint32_t)_mm_cvtsi128_si32(acc);
        memcpy(&C[i * 4], &r, 4);
    }
}
static inline void rct_matmul_l4rows_add(gf_t *C, const gf_t *A, const gf_t *B, int ad, int bd) {
    for (int i = 0; i < ad; ++i) {
        __m128i acc = _mm_setzero_si128();
        for (int k = 0; k < bd; ++k)
            acc = _mm_xor_si128(acc, _mm_shuffle_epi8(rct_vtl128(A[i * bd + k]),
                _mm_loadu_si128((const __m128i *)&B[k * 4])));
        uint32_t cur; memcpy(&cur, &C[i * 4], 4);
        cur ^= (uint32_t)_mm_cvtsi128_si32(acc);
        memcpy(&C[i * 4], &cur, 4);
    }
}

static inline void rct_gf4_matmul_add_sec(gf_t *acc, const gf_t *A, const gf_t *B) {
    static const _Alignas(16) int8_t BM[4][16] = {
        {0,1,2,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {4,5,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {8,9,10,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}};
    __m128i bv = _mm_loadu_si128((const __m128i *)B);
    __m128i brow[4];
    for (int k = 0; k < 4; ++k) brow[k] = _mm_shuffle_epi8(bv, _mm_load_si128((const __m128i *)BM[k]));
    for (int i = 0; i < 4; ++i) {
        __m128i a = _mm_setzero_si128();
        for (int k = 0; k < 4; ++k)
            a = _mm_xor_si128(a, rct_gf16_mul128_sec(_mm_set1_epi8((char)A[i * 4 + k]), brow[k]));
        uint32_t cur; memcpy(&cur, &acc[i * 4], 4);
        cur ^= (uint32_t)_mm_cvtsi128_si32(a);
        memcpy(&acc[i * 4], &cur, 4);
    }
}
static inline void rct_matmul_l4rows_sec(gf_t *C, const gf_t *A, const gf_t *B, int ad, int bd) {
    for (int i = 0; i < ad; ++i) {
        __m128i acc = _mm_setzero_si128();
        for (int k = 0; k < bd; ++k)
            acc = _mm_xor_si128(acc, rct_gf16_mul128_sec(_mm_set1_epi8((char)A[i * bd + k]),
                _mm_loadu_si128((const __m128i *)&B[k * 4])));
        uint32_t r = (uint32_t)_mm_cvtsi128_si32(acc);
        memcpy(&C[i * 4], &r, 4);
    }
}
#endif

static inline void rct_gauss_row_scale_sec(gf_t *row, gf_t s, int k0, int kend) {
    __m256i t = rct_bc_sec(s);
    for (int k = k0; k < kend; k += 32) {
        __m256i v = _mm256_loadu_si256((const __m256i *)(row + k));
        _mm256_storeu_si256((__m256i *)(row + k), _mm256_shuffle_epi8(t, v));
    }
}
static inline void rct_gauss_row_axpy_sec(gf_t *dst, const gf_t *src, gf_t s, int k0, int kend) {
    __m256i t = rct_bc_sec(s);
    for (int k = k0; k < kend; k += 32) {
        __m256i p = _mm256_shuffle_epi8(t, _mm256_loadu_si256((const __m256i *)(src + k)));
        __m256i d = _mm256_loadu_si256((const __m256i *)(dst + k));
        _mm256_storeu_si256((__m256i *)(dst + k), _mm256_xor_si256(d, p));
    }
}
#endif

#endif
