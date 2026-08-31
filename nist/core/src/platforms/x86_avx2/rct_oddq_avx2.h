#ifndef RCT_ODDQ_AVX2_H
#define RCT_ODDQ_AVX2_H

#if defined(__AVX2__) && (SNOVA_Q != 16) && (SNOVA_L == 4) && !defined(RCT_Q_SIMD_OFF)
#include <immintrin.h>
#define RCT_Q_SIMD 1
#define RCT_Q_LR16 (((SNOVA_lr) + 15) / 16)
#define RCT_Q_LRP  (RCT_Q_LR16 * 16)
static inline void rct_q_matmul4_add(uint16_t *acc16, const gf_t *b, const gf_t *c) {
    const __m256i bw = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)b));
    const __m256i cw = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)c));
    const __m256i m0 = _mm256_setr_epi8(0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9, 0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9);
    const __m256i m1 = _mm256_setr_epi8(2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11, 2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11);
    const __m256i m2 = _mm256_setr_epi8(4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13, 4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13);
    const __m256i m3 = _mm256_setr_epi8(6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15, 6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15);
    __m256i a = _mm256_loadu_si256((const __m256i *)acc16);
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m0), _mm256_permute4x64_epi64(cw, 0x00)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m1), _mm256_permute4x64_epi64(cw, 0x55)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m2), _mm256_permute4x64_epi64(cw, 0xAA)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m3), _mm256_permute4x64_epi64(cw, 0xFF)));
    _mm256_storeu_si256((__m256i *)acc16, a);
}
static inline void rct_q_matmul4T_add(uint16_t *acc16, const gf_t *b, const gf_t *c) {
    const __m256i bw = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)b));
    const __m256i cw = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)c));
    const __m256i tm = _mm256_setr_epi8(0,1,0,1,0,1,0,1,2,3,2,3,2,3,2,3, 4,5,4,5,4,5,4,5,6,7,6,7,6,7,6,7);
    __m256i a = _mm256_loadu_si256((const __m256i *)acc16);
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(_mm256_permute4x64_epi64(bw, 0x00), tm), _mm256_permute4x64_epi64(cw, 0x00)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(_mm256_permute4x64_epi64(bw, 0x55), tm), _mm256_permute4x64_epi64(cw, 0x55)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(_mm256_permute4x64_epi64(bw, 0xAA), tm), _mm256_permute4x64_epi64(cw, 0xAA)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(_mm256_permute4x64_epi64(bw, 0xFF), tm), _mm256_permute4x64_epi64(cw, 0xFF)));
    _mm256_storeu_si256((__m256i *)acc16, a);
}

#if   SNOVA_q == 11
#define RCT_Q_MAGIC_M 47663
#define RCT_Q_MAGIC_S 3
#define RCT_Q_HAVE_MAGIC 1
#elif SNOVA_q == 13
#define RCT_Q_MAGIC_M 20165
#define RCT_Q_MAGIC_S 2
#define RCT_Q_HAVE_MAGIC 1
#elif SNOVA_q == 19
#define RCT_Q_MAGIC_M 55189
#define RCT_Q_MAGIC_S 4
#define RCT_Q_HAVE_MAGIC 1
#else
#define RCT_Q_HAVE_MAGIC 0
#endif

#if RCT_Q_HAVE_MAGIC
static inline __m256i rct_q_barrett16(__m256i x) {
    __m256i quo = _mm256_srli_epi16(_mm256_mulhi_epu16(x, _mm256_set1_epi16((short)RCT_Q_MAGIC_M)),
                                    RCT_Q_MAGIC_S);
    return _mm256_sub_epi16(x, _mm256_mullo_epi16(quo, _mm256_set1_epi16((short)SNOVA_q)));
}

#endif

static inline __m256i rct_q_mm4_u16(__m256i bw, __m256i cw) {
    const __m256i m0 = _mm256_setr_epi8(0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9, 0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9);
    const __m256i m1 = _mm256_setr_epi8(2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11, 2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11);
    const __m256i m2 = _mm256_setr_epi8(4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13, 4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13);
    const __m256i m3 = _mm256_setr_epi8(6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15, 6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15);
    __m256i a = _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m0), _mm256_permute4x64_epi64(cw, 0x00));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m1), _mm256_permute4x64_epi64(cw, 0x55)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m2), _mm256_permute4x64_epi64(cw, 0xAA)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m3), _mm256_permute4x64_epi64(cw, 0xFF)));
    return a;
}

static inline __m256i rct_q_mm4_bc(__m256i bw, const __m256i *cperm) {
    const __m256i m0 = _mm256_setr_epi8(0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9, 0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9);
    const __m256i m1 = _mm256_setr_epi8(2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11, 2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11);
    const __m256i m2 = _mm256_setr_epi8(4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13, 4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13);
    const __m256i m3 = _mm256_setr_epi8(6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15, 6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15);
    __m256i a = _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m0), cperm[0]);
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m1), cperm[1]));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m2), cperm[2]));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(_mm256_shuffle_epi8(bw, m3), cperm[3]));
    return a;
}
static inline void rct_q_cperm(__m256i cw, __m256i *cperm) {
    cperm[0] = _mm256_permute4x64_epi64(cw, 0x00);
    cperm[1] = _mm256_permute4x64_epi64(cw, 0x55);
    cperm[2] = _mm256_permute4x64_epi64(cw, 0xAA);
    cperm[3] = _mm256_permute4x64_epi64(cw, 0xFF);
}
static inline __m256i rct_q_mm4_bs(const __m256i *bsh, __m256i cw) {
    __m256i a = _mm256_mullo_epi16(bsh[0], _mm256_permute4x64_epi64(cw, 0x00));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(bsh[1], _mm256_permute4x64_epi64(cw, 0x55)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(bsh[2], _mm256_permute4x64_epi64(cw, 0xAA)));
    a = _mm256_add_epi16(a, _mm256_mullo_epi16(bsh[3], _mm256_permute4x64_epi64(cw, 0xFF)));
    return a;
}
static inline void rct_q_bshuf(__m256i bw, __m256i *bsh) {
    const __m256i m0 = _mm256_setr_epi8(0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9, 0,1,0,1,0,1,0,1,8,9,8,9,8,9,8,9);
    const __m256i m1 = _mm256_setr_epi8(2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11, 2,3,2,3,2,3,2,3,10,11,10,11,10,11,10,11);
    const __m256i m2 = _mm256_setr_epi8(4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13, 4,5,4,5,4,5,4,5,12,13,12,13,12,13,12,13);
    const __m256i m3 = _mm256_setr_epi8(6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15, 6,7,6,7,6,7,6,7,14,15,14,15,14,15,14,15);
    bsh[0] = _mm256_shuffle_epi8(bw, m0);
    bsh[1] = _mm256_shuffle_epi8(bw, m1);
    bsh[2] = _mm256_shuffle_epi8(bw, m2);
    bsh[3] = _mm256_shuffle_epi8(bw, m3);
}

#if RCT_Q_HAVE_MAGIC && (SNOVA_r <= 8) && !defined(RCT_Q_MADD_OFF)
#define RCT_Q_MADD 1
#else
#define RCT_Q_MADD 0
#endif

#if RCT_Q_MADD
_Static_assert(2 * (SNOVA_q - 1) * (SNOVA_q - 1) < 32768, "i16 accumulation overflow guard");
static _Alignas(32) uint8_t rct_qv_wpair[SNOVA_n * 2 * 2 * RCT_Q_LRP];
static _Alignas(32) uint8_t rct_qv_s0all[SNOVA_n * 2 * 2 * RCT_Q_LRP];
static _Alignas(32) uint16_t rct_qv_sseg[SNOVA_l2][RCT_Q_LRP];
static _Alignas(32) uint8_t rct_qv_sigpat[RCT_Q_LR16][32];
static int rct_qv_built = 0;
static void rct_qv_build(void) {
    if (rct_qv_built) return;
    memset(rct_qv_wpair, 0, sizeof(rct_qv_wpair));
    memset(rct_qv_sseg, 0, sizeof(rct_qv_sseg));
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int k1 = 0; k1 < SNOVA_l; k1++)
            for (int ab = 0; ab < SNOVA_l; ab++)
                for (int j1 = 0; j1 < SNOVA_r; j1++)
                    rct_qv_sseg[i1 * SNOVA_l + k1][ab * SNOVA_r + j1] =
                        rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1];
    for (int g = 0; g < RCT_Q_LR16; g++)
        for (int t = 0; t < 16; t++) {
            int lane = g * 16 + t;
            int pos = (t & 7) * 2 + (t >> 3) * 16;
            rct_qv_sigpat[g][pos] = (lane < SNOVA_lr) ? (uint8_t)(lane % SNOVA_r) : 0x80;
            rct_qv_sigpat[g][pos + 1] = 0x80;
        }
    rct_qv_built = 1;
}
static inline __m128i rct_qv_pack16(__m256i v) {
    return _mm256_castsi256_si128(_mm256_permute4x64_epi64(
        _mm256_packus_epi16(v, _mm256_setzero_si256()), 0xD8));
}
#if RCT_Q_LR16 == 2
static inline __m256i rct_qv_pack32(__m256i v0, __m256i v1) {
    return _mm256_permute4x64_epi64(_mm256_packus_epi16(v0, v1), 0xD8);
}
static inline void rct_qv_ilv32(uint8_t *dst, __m256i A, __m256i B) {
    __m256i t0 = _mm256_unpacklo_epi8(A, B), t1 = _mm256_unpackhi_epi8(A, B);
    _mm256_store_si256((__m256i *)dst, _mm256_permute2x128_si256(t0, t1, 0x20));
    _mm256_store_si256((__m256i *)(dst + 32), _mm256_permute2x128_si256(t0, t1, 0x31));
}
#endif

#if defined(__AVXVNNI__) && !defined(RCT_Q_VNNI_OFF)
#define RCT_Q_VNNI 1
#else
#define RCT_Q_VNNI 0
#endif
#if RCT_Q_VNNI
static _Alignas(32) uint8_t rct_qv_wquad[SNOVA_n * 4 * RCT_Q_LRP];
static inline __m256i rct_qv_low16(__m256i a, __m256i b) {
    const __m256i p = _mm256_setr_epi8(0, 1, 4, 5, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
                                       0, 1, 4, 5, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1);
    __m256i x = _mm256_unpacklo_epi64(_mm256_shuffle_epi8(a, p), _mm256_shuffle_epi8(b, p));
    return _mm256_permute4x64_epi64(x, 0xD8);
}
#endif

#if RCT_Q_MADD && (SNOVA_r != SNOVA_l) && (SNOVA_r <= 7) && (SNOVA_r2 <= 64)
#define RCT_Q_EMM 1
#else
#define RCT_Q_EMM 0
#endif

#if RCT_Q_EMM
_Static_assert((SNOVA_alpha % 2) == 0, "SNOVA_alpha must be even");
static _Alignas(64) uint8_t rct_qv_s1p8[SNOVA_m1 * 8 * 128];
static _Alignas(64) uint8_t rct_qv_amt[SNOVA_o * SNOVA_alpha * 64];
static inline void rct_qv_tr8(uint8_t *dst, const uint8_t *src, int rstride) {
    __m128i l0 = _mm_loadl_epi64((const __m128i *)(src + 0 * rstride));
    __m128i l1 = _mm_loadl_epi64((const __m128i *)(src + 1 * rstride));
    __m128i l2 = _mm_loadl_epi64((const __m128i *)(src + 2 * rstride));
    __m128i l3 = _mm_loadl_epi64((const __m128i *)(src + 3 * rstride));
    __m128i l4 = (SNOVA_r > 4) ? _mm_loadl_epi64((const __m128i *)(src + 4 * rstride))
                               : _mm_setzero_si128();
    __m128i l5 = (SNOVA_r > 5) ? _mm_loadl_epi64((const __m128i *)(src + 5 * rstride))
                               : _mm_setzero_si128();
    __m128i l6 = (SNOVA_r > 6) ? _mm_loadl_epi64((const __m128i *)(src + 6 * rstride))
                               : _mm_setzero_si128();
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
#define RCT_QV_PA0 _mm256_setr_epi8( \
    0,-1,0,-1,0,-1,0,-1, 1,-1,1,-1,1,-1,1,-1, 2,-1,2,-1,2,-1,2,-1, 3,-1,3,-1,3,-1,3,-1)
#define RCT_QV_PA1 _mm256_setr_epi8( \
    4,-1,4,-1,4,-1,4,-1, 5,-1,5,-1,5,-1,5,-1, 6,-1,6,-1,6,-1,6,-1, 7,-1,7,-1,7,-1,7,-1)
#define RCT_QV_PB _mm256_setr_epi8( \
    0,-1,1,-1,2,-1,3,-1, 0,-1,1,-1,2,-1,3,-1, 0,-1,1,-1,2,-1,3,-1, 0,-1,1,-1,2,-1,3,-1)
static inline void rct_qv_mm_rx4(__m256i *o0, __m256i *o1, const uint8_t *AT, const uint8_t *B) {
    const __m256i pa0 = RCT_QV_PA0, pa1 = RCT_QV_PA1, pb = RCT_QV_PB;
    __m256i a0 = _mm256_setzero_si256(), a1 = _mm256_setzero_si256();
    for (int k = 0; k < SNOVA_r; k++) {
        long long aq;
        memcpy(&aq, AT + k * 8, 8);
        __m256i va = _mm256_set1_epi64x(aq);
        int32_t bd;
        memcpy(&bd, B + k * SNOVA_l, 4);
        __m256i vb = _mm256_shuffle_epi8(_mm256_set1_epi32(bd), pb);
        a0 = _mm256_add_epi16(a0, _mm256_mullo_epi16(_mm256_shuffle_epi8(va, pa0), vb));
        a1 = _mm256_add_epi16(a1, _mm256_mullo_epi16(_mm256_shuffle_epi8(va, pa1), vb));
    }
    *o0 = a0;
    *o1 = a1;
}
#endif
#endif
#else
#define RCT_Q_SIMD 0
#define RCT_Q_HAVE_MAGIC 0
#define RCT_Q_MADD 0
#endif
#ifndef RCT_Q_EMM
#define RCT_Q_EMM 0
#endif

#endif
