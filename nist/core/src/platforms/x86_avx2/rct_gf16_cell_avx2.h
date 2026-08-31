#ifndef RCT_GF16_CELL_AVX2_H
#define RCT_GF16_CELL_AVX2_H

#if RCT_USE_SIMD
static inline uint8_t rct_gfni_cleanup(uint8_t v) {
    return (uint8_t)((v ^ ((v & 0xf0) >> 3) ^ (v >> 4)) & 0x0f);
}
static inline __m256i rct_gfni_cleanup256(__m256i v) {
    const __m256i m0f = _mm256_set1_epi8(0x0f);
    __m256i vhi = _mm256_and_si256(v, _mm256_set1_epi8((char)0xf0));
    __m256i a = _mm256_srli_epi16(vhi, 3);
    __m256i b = _mm256_srli_epi16(v, 4);
    return _mm256_and_si256(_mm256_xor_si256(_mm256_xor_si256(v, a), b), m0f);
}
#if SNOVA_l == 4
static inline __m128i rct_gfni_cleanup128(__m128i v) {
    const __m128i m0f = _mm_set1_epi8(0x0f);
    __m128i vhi = _mm_and_si128(v, _mm_set1_epi8((char)0xf0));
    __m128i a = _mm_srli_epi16(vhi, 3);
    __m128i b = _mm_srli_epi16(v, 4);
    return _mm_and_si128(_mm_xor_si128(_mm_xor_si128(v, a), b), m0f);
}
#if RCT_USE_GFNI && (SNOVA_r == 4)
static inline __m256i rct_cm_cellmm256(__m256i A, __m256i B,
                                       const __m256i *SAY, const __m256i *SBY) {
    __m256i acc = _mm256_setzero_si256();
    for (int k = 0; k < SNOVA_l; ++k)
        acc = _mm256_xor_si256(acc, RCT_GFMUL256(
            _mm256_shuffle_epi8(A, SAY[k]), _mm256_shuffle_epi8(B, SBY[k])));
    return rct_gfni_cleanup256(acc);
}
#endif

static inline __m256i gf16_expand_u16x16(__m256i x) {
    __m256i v = _mm256_or_si256(_mm256_or_si256(x, _mm256_slli_epi16(x, 3)),
                                _mm256_or_si256(_mm256_slli_epi16(x, 6), _mm256_slli_epi16(x, 9)));
    return _mm256_and_si256(v, _mm256_set1_epi16(0x1111));
}
static inline __m256i gf16_compress_u16x16(__m256i a) {
    const __m256i m0f = _mm256_set1_epi16(0x000f);
    __m256i val = _mm256_xor_si256(
        _mm256_xor_si256(_mm256_and_si256(a, m0f),
                         _mm256_srli_epi16(_mm256_and_si256(a, _mm256_set1_epi16(0x00f0)), 3)),
        _mm256_xor_si256(_mm256_srli_epi16(_mm256_and_si256(a, _mm256_set1_epi16(0x0f00)), 6),
                         _mm256_srli_epi16(_mm256_and_si256(a, _mm256_set1_epi16((short)0xf000)), 9)));
    val = _mm256_xor_si256(_mm256_xor_si256(val,
              _mm256_srli_epi16(_mm256_and_si256(val, _mm256_set1_epi16(0x00f0)), 3)),
              _mm256_srli_epi16(val, 4));
    return _mm256_and_si256(val, m0f);
}
static inline __m128i gf16_pack_u16_to_bytes(__m256i c) {
    return _mm_packus_epi16(_mm256_castsi256_si128(c), _mm256_extracti128_si256(c, 1));
}
#define RCT_GF16_MM0 _mm256_setr_epi8(0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9, 0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9)
#define RCT_GF16_MM1 _mm256_setr_epi8(2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11, 2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11)
#define RCT_GF16_MM2 _mm256_setr_epi8(4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13, 4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13)
#define RCT_GF16_MM3 _mm256_setr_epi8(6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15, 6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15)
static inline void rct_gf16_cperm_exp(__m256i cw_raw, __m256i *cperm) {
    __m256i e = gf16_expand_u16x16(cw_raw);
    cperm[0] = _mm256_permute4x64_epi64(e, 0x00);
    cperm[1] = _mm256_permute4x64_epi64(e, 0x55);
    cperm[2] = _mm256_permute4x64_epi64(e, 0xAA);
    cperm[3] = _mm256_permute4x64_epi64(e, 0xFF);
}
static inline __m256i rct_gf16_mm4_bc(__m256i bw_raw, const __m256i *cperm_exp) {
    const __m256i m0 = RCT_GF16_MM0, m1 = RCT_GF16_MM1, m2 = RCT_GF16_MM2, m3 = RCT_GF16_MM3;
    __m256i a = _mm256_mullo_epi16(_mm256_shuffle_epi8(bw_raw, m0), cperm_exp[0]);
    a = _mm256_xor_si256(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw_raw, m1), cperm_exp[1]));
    a = _mm256_xor_si256(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw_raw, m2), cperm_exp[2]));
    a = _mm256_xor_si256(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw_raw, m3), cperm_exp[3]));
    return a;
}
static inline void rct_gf16_bshuf_exp(__m256i bw_raw, __m256i *bsh) {
    const __m256i m0 = RCT_GF16_MM0, m1 = RCT_GF16_MM1, m2 = RCT_GF16_MM2, m3 = RCT_GF16_MM3;
    __m256i e = gf16_expand_u16x16(bw_raw);
    bsh[0] = _mm256_shuffle_epi8(e, m0);
    bsh[1] = _mm256_shuffle_epi8(e, m1);
    bsh[2] = _mm256_shuffle_epi8(e, m2);
    bsh[3] = _mm256_shuffle_epi8(e, m3);
}
static inline __m256i rct_gf16_mm4_bs(const __m256i *bsh_exp, __m256i cw_raw) {
    __m256i a = _mm256_mullo_epi16(bsh_exp[0], _mm256_permute4x64_epi64(cw_raw, 0x00));
    a = _mm256_xor_si256(a, _mm256_mullo_epi16(bsh_exp[1], _mm256_permute4x64_epi64(cw_raw, 0x55)));
    a = _mm256_xor_si256(a, _mm256_mullo_epi16(bsh_exp[2], _mm256_permute4x64_epi64(cw_raw, 0xAA)));
    a = _mm256_xor_si256(a, _mm256_mullo_epi16(bsh_exp[3], _mm256_permute4x64_epi64(cw_raw, 0xFF)));
    return a;
}

#if defined(RCT_COEF_MULLO) && (RCT_COEF_MULLO + 0) \
    && (SNOVA_q == 16) && (SNOVA_r != SNOVA_l)
#define RCT_CM_ACTIVE 1
#endif
#if !defined(RCT_CM_ACTIVE) && defined(RCT_SQ_CM_S3) && (RCT_SQ_CM_S3 + 0) \
    && RCT_USE_SIMD && (SNOVA_q == 16) && (SNOVA_r == SNOVA_l) \
    && (SNOVA_L == 4)
#define RCT_CMS3_ONLY 1
#endif
#if defined(RCT_CM_ACTIVE) || defined(RCT_CMS3_ONLY)
static inline uint16_t rct_cm_exp(uint8_t a) {
    return (uint16_t)((a | ((uint16_t)a << 3) | ((uint16_t)a << 6) | ((uint16_t)a << 9)) & 0x1111);
}
static inline uint16_t rct_cm_cmp(uint16_t a) {
    uint16_t v = (uint16_t)((a & 0xf) ^ ((a & 0xf0) >> 3) ^ ((a & 0xf00) >> 6) ^ ((a & 0xf000) >> 9));
    return (uint16_t)((v ^ ((v & 0xf0) >> 3) ^ (v >> 4)) & 0xf);
}
static inline void rct_cm_expand_arr(uint16_t *dst, const uint8_t *src, int n) {
    int i = 0;
    for (; i + 16 <= n; i += 16)
        _mm256_storeu_si256((__m256i *)(dst + i),
            gf16_expand_u16x16(_mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(src + i)))));
    for (; i < n; ++i) dst[i] = rct_cm_exp(src[i]);
}
#endif
#endif
#endif

#endif
