#ifndef RCT_VERIFY_KERNELS_AVX2_H
#define RCT_VERIFY_KERNELS_AVX2_H

#if RCT_VF_MTK2 || RCT_VF_JOG
static _Alignas(64) uint8_t rct_mtk2[256][16];
static void rct_build_mtk2(void) {
    for (int idx = 0; idx < 256; ++idx)
        for (int x = 0; x < 16; ++x)
            rct_mtk2[idx][x] = (uint8_t)(rct_multtab[(idx & 0x0F) * SNOVA_q + x] |
                                         (rct_multtab[(idx >> 4) * SNOVA_q + x] << 4));
}
static inline __m256i rct_mtk2t(uint8_t idx) {
    return _mm256_broadcastsi128_si256(_mm_load_si128((const __m128i *)rct_mtk2[idx]));
}
static inline __m256i rct_mtk2t16(uint16_t idx16) {
    return _mm256_broadcastsi128_si256(
        _mm_load_si128((const __m128i *)((const uint8_t *)rct_mtk2 + idx16)));
}
static inline __m256i rct_nib_lo(__m256i v) { return _mm256_and_si256(v, _mm256_set1_epi8(0x0f)); }
static inline __m256i rct_nib_hi(__m256i v) {
    return _mm256_and_si256(_mm256_srli_epi16(v, 4), _mm256_set1_epi8(0x0f));
}
static inline void rct_vf_expand_sig(gf_t *out, const uint8_t *in, size_t num) {
    const __m128i m0f = _mm_set1_epi8(0x0f);
    size_t nb = num / 2, i = 0;
    for (; i + 16 <= nb; i += 16) {
        __m128i b = _mm_loadu_si128((const __m128i *)(in + i));
        __m128i lo = _mm_and_si128(b, m0f);
        __m128i hi = _mm_and_si128(_mm_srli_epi16(b, 4), m0f);
        _mm_storeu_si128((__m128i *)(out + 2 * i), _mm_unpacklo_epi8(lo, hi));
        _mm_storeu_si128((__m128i *)(out + 2 * i + 16), _mm_unpackhi_epi8(lo, hi));
    }
    for (; i < nb; ++i) {
        out[2 * i] = (gf_t)(in[i] & 0x0F);
        out[2 * i + 1] = (gf_t)(in[i] >> 4);
    }
}
#endif

#if RCT_VF_MTK2
static inline __m128i rct_vf_pack_pair32(__m256i v, __m256i pl, __m256i ph) {
    __m256i lo = _mm256_shuffle_epi8(v, pl);
    __m256i hi = _mm256_shuffle_epi8(v, ph);
    __m256i pk = _mm256_or_si256(lo, _mm256_slli_epi16(hi, 4));
    return _mm256_castsi256_si128(_mm256_permute4x64_epi64(pk, 0x08));
}
static const _Alignas(32) uint8_t RCT_VF_PPL[32] = {
    0, 8, 1, 9, 2, 10, 3, 11, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
    0, 8, 1, 9, 2, 10, 3, 11, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80};
static const _Alignas(32) uint8_t RCT_VF_PPH[32] = {
    4, 12, 5, 13, 6, 14, 7, 15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
    4, 12, 5, 13, 6, 14, 7, 15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80};
static const _Alignas(32) uint8_t RCT_VF_WPL[32] = {
    0, 2, 4, 6, 8, 10, 12, 14, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
    0, 2, 4, 6, 8, 10, 12, 14, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80};
static const _Alignas(32) uint8_t RCT_VF_WPH[32] = {
    1, 3, 5, 7, 9, 11, 13, 15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
    1, 3, 5, 7, 9, 11, 13, 15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80};

#if RCT_USE_GFNI
static _Alignas(32) uint8_t rct_whipM[SNOVA_l2][32];
static _Alignas(32) uint8_t rct_whipR[SNOVA_l][32];
static _Alignas(32) uint8_t rct_vf_dwexp[32];
static void rct_build_vf_gfni(void) {
    memset(rct_whipM, 0, sizeof(rct_whipM));
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int k1 = 0; k1 < SNOVA_l; k1++)
            for (int ab = 0; ab < SNOVA_l; ab++)
                for (int j = 0; j < SNOVA_r; j++)
                    rct_whipM[i1 * SNOVA_l + k1][ab * SNOVA_r + j] =
                        rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1];
    for (int k1 = 0; k1 < SNOVA_l; k1++)
        for (int b = 0; b < 32; b++)
            rct_whipR[k1][b] = (b < SNOVA_lr) ? (uint8_t)((k1 & 1) * SNOVA_r + (b % SNOVA_r)) : 0x80;
    for (int b = 0; b < 32; b++)
        rct_vf_dwexp[b] = (b / 4 < SNOVA_r) ? (uint8_t)(b / 4) : 0x80;
}
static inline __m256i rct_vf_mm_dw(const uint8_t *Acm, const gf_t *B) {
    const __m256i pat = _mm256_load_si256((const __m256i *)rct_vf_dwexp);
    __m256i acc = _mm256_setzero_si256();
    for (int k = 0; k < SNOVA_r; k++) {
        int64_t aw;
        memcpy(&aw, Acm + k * 8, 8);
        __m256i acol = _mm256_shuffle_epi8(_mm256_set1_epi64x(aw), pat);
        int32_t bw;
        memcpy(&bw, B + k * SNOVA_l, 4);
        acc = _mm256_xor_si256(acc, _mm256_gf2p8mul_epi8(acol, _mm256_set1_epi32(bw)));
    }
    return acc;
}
static inline void rct_vf_tr8(uint8_t *dst, const uint8_t *src, int rstride) {
    __m128i l0 = _mm_loadl_epi64((const __m128i *)(src + 0 * rstride));
    __m128i l1 = _mm_loadl_epi64((const __m128i *)(src + 1 * rstride));
    __m128i l2 = _mm_loadl_epi64((const __m128i *)(src + 2 * rstride));
    __m128i l3 = _mm_loadl_epi64((const __m128i *)(src + 3 * rstride));
    __m128i l4 = _mm_loadl_epi64((const __m128i *)(src + 4 * rstride));
    __m128i l5 = _mm_loadl_epi64((const __m128i *)(src + 5 * rstride));
    __m128i l6 = _mm_loadl_epi64((const __m128i *)(src + 6 * rstride));
    __m128i x0 = _mm_unpacklo_epi64(l0, l1);
    __m128i x1 = _mm_unpacklo_epi64(l2, l3);
    __m128i x2 = _mm_unpacklo_epi64(l4, l5);
    __m128i x3 = _mm_unpacklo_epi64(l6, _mm_setzero_si128());
    __m128i u0 = _mm_unpacklo_epi8(x0, x1);
    __m128i u1 = _mm_unpackhi_epi8(x0, x1);
    __m128i u2 = _mm_unpacklo_epi8(x2, x3);
    __m128i u3 = _mm_unpackhi_epi8(x2, x3);
    __m128i v0 = _mm_unpacklo_epi8(u0, u1);
    __m128i v1 = _mm_unpackhi_epi8(u0, u1);
    __m128i v2 = _mm_unpacklo_epi8(u2, u3);
    __m128i v3 = _mm_unpackhi_epi8(u2, u3);
    _mm_store_si128((__m128i *)(dst + 0), _mm_unpacklo_epi32(v0, v2));
    _mm_store_si128((__m128i *)(dst + 16), _mm_unpackhi_epi32(v0, v2));
    _mm_store_si128((__m128i *)(dst + 32), _mm_unpacklo_epi32(v1, v3));
    _mm_store_si128((__m128i *)(dst + 48), _mm_unpackhi_epi32(v1, v3));
}

static __m256i rct_vf_thv[4];
static __m256i rct_vf_aqko;
static __m256i rct_vf_aqlo;
static __m256i rct_vf_madp;
__attribute__((unused)) static void rct_build_vf_aq(void) {
    for (int t = 0; t < 4; ++t) {
        uint8_t colb[8];
        for (int j = 0; j < 8; ++j) {
            uint8_t v = (uint8_t)(1u << (j & 3));
            uint8_t rh = 0;
            for (int k = 0; k < 4; ++k)
                rh |= (uint8_t)(((rct_multtab[v * SNOVA_q + (1u << k)] >> t) & 1) << k);
            colb[j] = (j < 4) ? rh : (uint8_t)(rh << 4);
        }
        uint64_t qw = 0;
        for (int i = 0; i < 8; ++i) {
            uint8_t row = 0;
            for (int j = 0; j < 8; ++j) row |= (uint8_t)(((colb[j] >> i) & 1) << j);
            qw |= (uint64_t)row << (8 * (7 - i));
        }
        rct_vf_thv[t] = _mm256_set1_epi64x((long long)qw);
    }
    _Alignas(32) uint8_t ko[32], lo[32];
    static const uint8_t ord8[8] = {2, 0, 6, 4, 3, 1, 7, 5};
    for (int b = 0; b < 32; ++b) {
        ko[b] = (uint8_t)((b & 8) | ord8[b & 7]);
        lo[b] = (uint8_t)((b & 16) | ((b & 15) ^ 1));
    }
    rct_vf_aqko = _mm256_load_si256((const __m256i *)ko);
    rct_vf_aqlo = _mm256_load_si256((const __m256i *)lo);
    rct_vf_madp = _mm256_set1_epi16(0x1001);
}
static inline void rct_vf_aq_tree(__m256i keys, __m256i out[4]) {
    __m256i R0 = _mm256_gf2p8affine_epi64_epi8(keys, rct_vf_thv[0], 0);
    __m256i R1 = _mm256_gf2p8affine_epi64_epi8(keys, rct_vf_thv[1], 0);
    __m256i R2 = _mm256_gf2p8affine_epi64_epi8(keys, rct_vf_thv[2], 0);
    __m256i R3 = _mm256_gf2p8affine_epi64_epi8(keys, rct_vf_thv[3], 0);
    __m256i u32lo = _mm256_unpacklo_epi8(R3, R2);
    __m256i u10lo = _mm256_unpacklo_epi8(R1, R0);
    __m256i u32hi = _mm256_unpackhi_epi8(R3, R2);
    __m256i u10hi = _mm256_unpackhi_epi8(R1, R0);
    out[0] = _mm256_unpacklo_epi16(u32lo, u10lo);
    out[1] = _mm256_unpackhi_epi16(u32lo, u10lo);
    out[2] = _mm256_unpacklo_epi16(u32hi, u10hi);
    out[3] = _mm256_unpackhi_epi16(u32hi, u10hi);
}
static inline void rct_vf_aq_quad(const uint8_t *pc, uint8_t *dst) {
    __m256i c01 = _mm256_loadu_si256((const __m256i *)pc);
    __m256i c23 = _mm256_loadu_si256((const __m256i *)(pc + 32));
    __m256i k01 = _mm256_maddubs_epi16(c01, rct_vf_madp);
    __m256i k23 = _mm256_maddubs_epi16(c23, rct_vf_madp);
    __m256i keys = _mm256_shuffle_epi8(_mm256_packus_epi16(k01, k23), rct_vf_aqko);
    __m256i q[4];
    rct_vf_aq_tree(keys, q);
    _mm256_store_si256((__m256i *)(dst + 0), q[0]);
    _mm256_store_si256((__m256i *)(dst + 32), q[1]);
    _mm256_store_si256((__m256i *)(dst + 64), q[2]);
    _mm256_store_si256((__m256i *)(dst + 96), q[3]);
}
#define RCT_VF_AQ_LOFF(p) (32 * (((p) & 7) >> 1) + 16 * ((p) >> 3) + 8 * ((p) & 1))
static inline __m256i rct_vf_aq_bq(const uint8_t *base, int off) {
    int64_t w;
    memcpy(&w, base + off, 8);
    return _mm256_set1_epi64x(w);
}
#endif

#endif

#endif
