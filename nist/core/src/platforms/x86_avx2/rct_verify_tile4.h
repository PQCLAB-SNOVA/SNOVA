#ifndef RCT_VERIFY_TILE4_H
#define RCT_VERIFY_TILE4_H

#if RCT_VF_TILE4

#ifdef RCT_T4_SELFTEST
#include <stdio.h>
#endif

#define RCT_T4_NB ((RCT_JOG_NL + 3) / 4)
#define RCT_T4_RP4 (RCT_T4_NB * 4)
#define RCT_T4_NQ ((RCT_JOG_NL + 15) / 16)
#define RCT_T4_CP (RCT_T4_NQ * 16)
#define RCT_T4_NCB (RCT_T4_NQ * 4)
#define RCT_T4_NP ((SNOVA_lr + 1) / 2)
#define RCT_T4_LH SNOVA_lr16
_Static_assert(RCT_T4_CP <= RCT_JOG_L32, "rct_vf_tile4: column pad must sit inside jogress zero pad");
_Static_assert(SNOVA_lr <= 64, "rct_vf_tile4: lane axis assumes lr <= 64 (RCT_T4_LH <= 2 lr32 halves)");
_Static_assert(RCT_T4_NP <= 32, "rct_vf_tile4: Left dispatch + dual-LOFF-block addressing cover <=32 pairs");
_Static_assert(!RCT_USE_SIMD, "rct_vf_tile4: whip_flat copies the l!=4 in-kernel whip branch");

static const _Alignas(32) uint8_t rct_t4_zrow[RCT_JOG_L32] = {0};

static __m256i rct_t4_thv[4];
static __m256i rct_t4_aqko;
static __m256i rct_t4_aqlo;
static __m256i rct_t4_madp;
static void rct_build_t4(void) {
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
        rct_t4_thv[t] = _mm256_set1_epi64x((long long)qw);
    }
    _Alignas(32) uint8_t ko[32], lo[32];
    static const uint8_t ord8[8] = {2, 0, 6, 4, 3, 1, 7, 5};
    for (int b = 0; b < 32; ++b) {
        ko[b] = (uint8_t)((b & 8) | ord8[b & 7]);
        lo[b] = (uint8_t)((b & 16) | ((b & 15) ^ 1));
    }
    rct_t4_aqko = _mm256_load_si256((const __m256i *)ko);
    rct_t4_aqlo = _mm256_load_si256((const __m256i *)lo);
    rct_t4_madp = _mm256_set1_epi16(0x1001);
}
static inline void rct_t4_tree(__m256i keys, __m256i out[4]) {
    __m256i R0 = _mm256_gf2p8affine_epi64_epi8(keys, rct_t4_thv[0], 0);
    __m256i R1 = _mm256_gf2p8affine_epi64_epi8(keys, rct_t4_thv[1], 0);
    __m256i R2 = _mm256_gf2p8affine_epi64_epi8(keys, rct_t4_thv[2], 0);
    __m256i R3 = _mm256_gf2p8affine_epi64_epi8(keys, rct_t4_thv[3], 0);
    __m256i u32lo = _mm256_unpacklo_epi8(R3, R2);
    __m256i u10lo = _mm256_unpacklo_epi8(R1, R0);
    __m256i u32hi = _mm256_unpackhi_epi8(R3, R2);
    __m256i u10hi = _mm256_unpackhi_epi8(R1, R0);
    out[0] = _mm256_unpacklo_epi16(u32lo, u10lo);
    out[1] = _mm256_unpackhi_epi16(u32lo, u10lo);
    out[2] = _mm256_unpacklo_epi16(u32hi, u10hi);
    out[3] = _mm256_unpackhi_epi16(u32hi, u10hi);
}
static inline void rct_t4_quad(const uint8_t *const rp[4], int coff, uint8_t *dst) {
    __m128i r0 = _mm_loadu_si128((const __m128i *)(rp[0] + coff));
    __m128i r1 = _mm_loadu_si128((const __m128i *)(rp[1] + coff));
    __m128i r2 = _mm_loadu_si128((const __m128i *)(rp[2] + coff));
    __m128i r3 = _mm_loadu_si128((const __m128i *)(rp[3] + coff));
    __m128i t0 = _mm_unpacklo_epi32(r0, r1), t1 = _mm_unpackhi_epi32(r0, r1);
    __m128i t2 = _mm_unpacklo_epi32(r2, r3), t3 = _mm_unpackhi_epi32(r2, r3);
    __m256i c01 = _mm256_set_m128i(_mm_unpackhi_epi64(t0, t2), _mm_unpacklo_epi64(t0, t2));
    __m256i c23 = _mm256_set_m128i(_mm_unpackhi_epi64(t1, t3), _mm_unpacklo_epi64(t1, t3));
    __m256i k01 = _mm256_maddubs_epi16(c01, rct_t4_madp);
    __m256i k23 = _mm256_maddubs_epi16(c23, rct_t4_madp);
    __m256i keys = _mm256_shuffle_epi8(_mm256_packus_epi16(k01, k23), rct_t4_aqko);
    __m256i q[4];
    rct_t4_tree(keys, q);
    _mm256_store_si256((__m256i *)(dst + 0), q[0]);
    _mm256_store_si256((__m256i *)(dst + 32), q[1]);
    _mm256_store_si256((__m256i *)(dst + 64), q[2]);
    _mm256_store_si256((__m256i *)(dst + 96), q[3]);
}
#define RCT_T4_LOFF(p) (32 * (((p) & 7) >> 1) + 16 * ((p) >> 3) + 8 * ((p) & 1))
#define RCT_T4_LQOFF(p) (((p) >> 4) * 128 + RCT_T4_LOFF((p) & 15))
static inline __m256i rct_t4_bq(const uint8_t *base, int off) {
    int64_t w;
    memcpy(&w, base + off, 8);
    return _mm256_set1_epi64x(w);
}

static inline __attribute__((always_inline)) void rct_t4_left(
        const uint8_t *st0p_t, const uint8_t *t4_lqt, gf_t *st1_mi,
        const int g0, const int gw) {
#if RCT_T4_LH == 1
    __m256i acc[5];
    for (int pi = 0; pi < gw; ++pi) acc[pi] = _mm256_setzero_si256();
    for (int tb = 0; tb < RCT_T4_NB; ++tb)
        for (int kp = 0; kp < 2; ++kp) {
            const __m256i s0 =
                _mm256_load_si256((const __m256i *)&st0p_t[((size_t)tb * 2 + kp) * 32]);
            const uint8_t *lb = &t4_lqt[((size_t)tb * 2 + kp) * 128];
            for (int pi = 0; pi < gw; ++pi)
                acc[pi] = _mm256_xor_si256(acc[pi], _mm256_gf2p8affine_epi64_epi8(
                              s0, rct_t4_bq(lb, RCT_T4_LOFF(g0 + pi)), 0));
        }
    for (int pi = 0; pi < gw; ++pi) {
        _Alignas(32) uint8_t rowlo[32], rowhi[32];
        _mm256_store_si256((__m256i *)rowlo, rct_nib_lo(acc[pi]));
        _mm256_store_si256((__m256i *)rowhi, rct_nib_hi(acc[pi]));
#else
    __m256i acc[5][RCT_T4_LH];
    for (int pi = 0; pi < gw; ++pi)
        for (int h = 0; h < RCT_T4_LH; ++h) acc[pi][h] = _mm256_setzero_si256();
    for (int tb = 0; tb < RCT_T4_NB; ++tb)
        for (int kp = 0; kp < 2; ++kp) {
            __m256i s0[RCT_T4_LH];
            for (int h = 0; h < RCT_T4_LH; ++h)
                s0[h] = _mm256_load_si256((const __m256i *)
                            &st0p_t[(((size_t)tb * 2 + kp) * RCT_T4_LH + h) * 32]);
            const uint8_t *lb = &t4_lqt[((size_t)tb * 2 + kp) * (RCT_T4_LH * 128)];
            for (int pi = 0; pi < gw; ++pi) {
                const __m256i mq = rct_t4_bq(lb, RCT_T4_LQOFF(g0 + pi));
                for (int h = 0; h < RCT_T4_LH; ++h)
                    acc[pi][h] = _mm256_xor_si256(acc[pi][h],
                                     _mm256_gf2p8affine_epi64_epi8(s0[h], mq, 0));
            }
        }
    for (int pi = 0; pi < gw; ++pi) {
        _Alignas(32) uint8_t rowlo[SNOVA_lr32], rowhi[SNOVA_lr32];
        for (int h = 0; h < RCT_T4_LH; ++h) {
            _mm256_store_si256((__m256i *)(rowlo + h * 32), rct_nib_lo(acc[pi][h]));
            _mm256_store_si256((__m256i *)(rowhi + h * 32), rct_nib_hi(acc[pi][h]));
        }
#endif
        for (int half = 0; half < 2; ++half) {
            const int jrow = 2 * (g0 + pi) + half;
            if (jrow >= SNOVA_lr) break;
            const int a1 = jrow / SNOVA_r, i1 = jrow % SNOVA_r;
            const uint8_t *row = half ? rowhi : rowlo;
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int j1 = 0; j1 < SNOVA_r; ++j1)
                    st1_mi[((size_t)a1 * SNOVA_l + b1) * SNOVA_r2 + (size_t)i1 * SNOVA_r + j1]
                        = row[(size_t)b1 * SNOVA_r + j1];
        }
    }
}

#ifdef RCT_T4_SELFTEST
static uint8_t rct_t4_affb(uint64_t qm, uint8_t src) {
    uint8_t out = 0;
    for (int i = 0; i < 8; ++i) {
        const uint8_t row = (uint8_t)(qm >> (8 * (7 - i)));
        out |= (uint8_t)((__builtin_parity((unsigned)(row & src)) & 1) << i);
    }
    return out;
}
static int rct_t4_ckmat(uint64_t qm, gf_t a, gf_t b, gf_t cq, gf_t d) {
    for (int s = 0; s < 256; ++s) {
        const gf_t in0 = (gf_t)(s & 15), in1 = (gf_t)(s >> 4);
        const uint8_t want = (uint8_t)((gf_mult(a, in0) ^ gf_mult(b, in1)) |
                                       ((gf_mult(cq, in0) ^ gf_mult(d, in1)) << 4));
        if (rct_t4_affb(qm, (uint8_t)s) != want) return 1;
    }
    return 0;
}
#endif

#if RCT_TILE4_PREBAKE
static _Alignas(32) uint8_t rct_t4_qall[(size_t)SNOVA_m1 * RCT_T4_NB * RCT_T4_NQ * 128];

static void rct_t4_prebake(const rct_pk_t *pkx) {
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        const uint8_t *PJ = &pkx->P[(size_t)mi * RCT_JOG_NL * RCT_JOG_L32];
        for (int tb = 0; tb < RCT_T4_NB; ++tb) {
            const uint8_t *rp[4];
            for (int k = 0; k < 4; ++k) {
                const int p = 4 * tb + k;
                rp[k] = (p < RCT_JOG_NL) ? (PJ + (size_t)p * RCT_JOG_L32) : rct_t4_zrow;
            }
            uint8_t *dst = rct_t4_qall + ((size_t)mi * RCT_T4_NB + tb) * RCT_T4_NQ * 128;
            for (int q4 = 0; q4 < RCT_T4_NQ; ++q4)
                rct_t4_quad(rp, 16 * q4, dst + (size_t)q4 * 128);
        }
    }
}
#endif

static void rct_vf_tile4(rct_vf_ctx *c) {
    RCT_SCRATCH _Alignas(32) uint8_t t4_whip[RCT_T4_CP * SNOVA_lr32];
    memset(t4_whip, 0, sizeof(t4_whip));
    for (int ab = 0; ab < SNOVA_l; ++ab)
        for (int idx = 0; idx < SNOVA_n; ++idx) {
            const gf_t *sig = &c->sig_gf[(size_t)idx * SNOVA_lr];
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                const gf_t *Srow = &rct_S[ab * SNOVA_l2 + i1 * SNOVA_l];
                __m128i acc = _mm_setzero_si128();
                for (int k1 = 0; k1 < SNOVA_l; ++k1)
                    acc = _mm_xor_si128(acc,
                        rct_jog_sv128(Srow[k1],
                                      _mm_loadu_si128((const __m128i *)&sig[k1 * SNOVA_r])));
                _Alignas(16) uint8_t wb[16];
                _mm_store_si128((__m128i *)wb, acc);
                memcpy(&t4_whip[((size_t)idx * SNOVA_l + i1) * SNOVA_lr32
                                + (size_t)ab * SNOVA_r],
                       wb, SNOVA_r);
            }
        }

    RCT_SCRATCH _Alignas(32) uint8_t t4_wpk[RCT_T4_NCB * 2 * SNOVA_lr32];
    RCT_SCRATCH _Alignas(32) uint8_t t4_lqt[RCT_T4_NCB * 2 * (RCT_T4_LH * 128)];
#if RCT_T4_LH == 1
    for (int tb = 0; tb < RCT_T4_NCB; ++tb) {
        __m256i r0 = _mm256_load_si256((const __m256i *)&t4_whip[((size_t)tb * 4 + 0) * SNOVA_lr32]);
        __m256i r1 = _mm256_load_si256((const __m256i *)&t4_whip[((size_t)tb * 4 + 1) * SNOVA_lr32]);
        __m256i r2 = _mm256_load_si256((const __m256i *)&t4_whip[((size_t)tb * 4 + 2) * SNOVA_lr32]);
        __m256i r3 = _mm256_load_si256((const __m256i *)&t4_whip[((size_t)tb * 4 + 3) * SNOVA_lr32]);
        __m256i p0 = _mm256_or_si256(r0, _mm256_slli_epi16(r1, 4));
        __m256i p1 = _mm256_or_si256(r2, _mm256_slli_epi16(r3, 4));
        _mm256_store_si256((__m256i *)&t4_wpk[((size_t)tb * 2 + 0) * 32], p0);
        _mm256_store_si256((__m256i *)&t4_wpk[((size_t)tb * 2 + 1) * 32], p1);
        __m256i q[4];
        rct_t4_tree(_mm256_shuffle_epi8(p0, rct_t4_aqlo), q);
        _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 0) * 128 + 0], q[0]);
        _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 0) * 128 + 32], q[1]);
        _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 0) * 128 + 64], q[2]);
        _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 0) * 128 + 96], q[3]);
        rct_t4_tree(_mm256_shuffle_epi8(p1, rct_t4_aqlo), q);
        _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 1) * 128 + 0], q[0]);
        _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 1) * 128 + 32], q[1]);
        _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 1) * 128 + 64], q[2]);
        _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 1) * 128 + 96], q[3]);
    }
#else
    for (int tb = 0; tb < RCT_T4_NCB; ++tb)
        for (int h = 0; h < RCT_T4_LH; ++h) {
            __m256i r0 = _mm256_load_si256((const __m256i *)&t4_whip[((size_t)tb * 4 + 0) * SNOVA_lr32 + (size_t)h * 32]);
            __m256i r1 = _mm256_load_si256((const __m256i *)&t4_whip[((size_t)tb * 4 + 1) * SNOVA_lr32 + (size_t)h * 32]);
            __m256i r2 = _mm256_load_si256((const __m256i *)&t4_whip[((size_t)tb * 4 + 2) * SNOVA_lr32 + (size_t)h * 32]);
            __m256i r3 = _mm256_load_si256((const __m256i *)&t4_whip[((size_t)tb * 4 + 3) * SNOVA_lr32 + (size_t)h * 32]);
            __m256i p0 = _mm256_or_si256(r0, _mm256_slli_epi16(r1, 4));
            __m256i p1 = _mm256_or_si256(r2, _mm256_slli_epi16(r3, 4));
            _mm256_store_si256((__m256i *)&t4_wpk[(((size_t)tb * 2 + 0) * RCT_T4_LH + h) * 32], p0);
            _mm256_store_si256((__m256i *)&t4_wpk[(((size_t)tb * 2 + 1) * RCT_T4_LH + h) * 32], p1);
            __m256i q[4];
            rct_t4_tree(_mm256_shuffle_epi8(p0, rct_t4_aqlo), q);
            _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 0) * (RCT_T4_LH * 128) + (size_t)h * 128 + 0], q[0]);
            _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 0) * (RCT_T4_LH * 128) + (size_t)h * 128 + 32], q[1]);
            _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 0) * (RCT_T4_LH * 128) + (size_t)h * 128 + 64], q[2]);
            _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 0) * (RCT_T4_LH * 128) + (size_t)h * 128 + 96], q[3]);
            rct_t4_tree(_mm256_shuffle_epi8(p1, rct_t4_aqlo), q);
            _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 1) * (RCT_T4_LH * 128) + (size_t)h * 128 + 0], q[0]);
            _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 1) * (RCT_T4_LH * 128) + (size_t)h * 128 + 32], q[1]);
            _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 1) * (RCT_T4_LH * 128) + (size_t)h * 128 + 64], q[2]);
            _mm256_store_si256((__m256i *)&t4_lqt[((size_t)tb * 2 + 1) * (RCT_T4_LH * 128) + (size_t)h * 128 + 96], q[3]);
        }
#endif

#ifdef RCT_T4_SELFTEST
    {
        static int t4_ut_done = 0;
        if (!t4_ut_done) {
            t4_ut_done = 1;
            size_t bad_wpk = 0, bad_pad = 0, bad_lqt = 0, bad_quad = 0;
            for (int tb = 0; tb < RCT_T4_NCB; ++tb)
                for (int kp = 0; kp < 2; ++kp) {
                    const int c0 = 4 * tb + 2 * kp;
                    const int rowpad = (c0 >= RCT_JOG_NL);
                    for (int j = 0; j < SNOVA_lr32; ++j) {
                        const uint8_t got = t4_wpk[((size_t)tb * 2 + kp) * SNOVA_lr32 + j];
                        const uint8_t want = (uint8_t)(t4_whip[(size_t)c0 * SNOVA_lr32 + j] |
                                                       (t4_whip[((size_t)c0 + 1) * SNOVA_lr32 + j] << 4));
                        if (got != want) ++bad_wpk;
                        if ((rowpad || j >= SNOVA_lr) && got != 0) ++bad_pad;
                    }
                    if (rowpad) {
                        for (int j = 0; j < RCT_T4_LH * 128; ++j)
                            if (t4_lqt[((size_t)tb * 2 + kp) * (RCT_T4_LH * 128) + j] != 0) ++bad_pad;
                    }
                    for (int p = 0; p < RCT_T4_NP; ++p) {
                        uint64_t qm;
                        memcpy(&qm, &t4_lqt[((size_t)tb * 2 + kp) * (RCT_T4_LH * 128) + RCT_T4_LQOFF(p)], 8);
                        bad_lqt += (size_t)rct_t4_ckmat(qm,
                            t4_whip[(size_t)c0 * SNOVA_lr32 + 2 * p],
                            t4_whip[((size_t)c0 + 1) * SNOVA_lr32 + 2 * p],
                            t4_whip[(size_t)c0 * SNOVA_lr32 + 2 * p + 1],
                            t4_whip[((size_t)c0 + 1) * SNOVA_lr32 + 2 * p + 1]);
                    }
                }
            {
                _Alignas(32) uint8_t qd[128];
                const int mis[2] = {0, SNOVA_m1 - 1};
                const int tbs[2] = {0, RCT_T4_NB - 1};
                const int qs[2] = {0, RCT_T4_NQ - 1};
                for (int t = 0; t < 2; ++t) {
                    const uint8_t *PJt = &c->pkx->P[(size_t)mis[t] * RCT_JOG_NL * RCT_JOG_L32];
                    const uint8_t *rp[4];
                    for (int k = 0; k < 4; ++k) {
                        const int p = 4 * tbs[t] + k;
                        rp[k] = (p < RCT_JOG_NL) ? (PJt + (size_t)p * RCT_JOG_L32) : rct_t4_zrow;
                    }
                    const int coff = 16 * qs[t];
                    rct_t4_quad(rp, coff, qd);
                    for (int cc = 0; cc < 4; ++cc)
                        for (int ip = 0; ip < 2; ++ip)
                            for (int kp = 0; kp < 2; ++kp) {
                                uint64_t qm;
                                memcpy(&qm, qd + (cc >> 1) * 64 + kp * 32 + (cc & 1) * 16 + ip * 8, 8);
                                bad_quad += (size_t)rct_t4_ckmat(qm,
                                    rp[2 * ip][coff + 4 * cc + 2 * kp],
                                    rp[2 * ip][coff + 4 * cc + 2 * kp + 1],
                                    rp[2 * ip + 1][coff + 4 * cc + 2 * kp],
                                    rp[2 * ip + 1][coff + 4 * cc + 2 * kp + 1]);
                            }
                }
            }
            printf("[RCT_T4_UNIT] wpk=%zu pad=%zu lqt=%zu quad=%zu (all 0 = PASS)\n",
                   bad_wpk, bad_pad, bad_lqt, bad_quad);
        }
    }
#endif

    gf_t *sum_t1 = c->sum_t1;
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        const uint8_t *PJ = &c->pkx->P[(size_t)mi * RCT_JOG_NL * RCT_JOG_L32];
        (void)PJ;

        RCT_SCRATCH _Alignas(32) uint8_t st0p_t[RCT_T4_NB * 2 * SNOVA_lr32];
        for (int tb = 0; tb < RCT_T4_NB; ++tb) {
#if !RCT_TILE4_PREBAKE
            const uint8_t *rp[4];
            for (int k = 0; k < 4; ++k) {
                const int p = 4 * tb + k;
                rp[k] = (p < RCT_JOG_NL) ? (PJ + (size_t)p * RCT_JOG_L32) : rct_t4_zrow;
            }
            _Alignas(32) uint8_t qrow[RCT_T4_NQ * 128];
            for (int q4 = 0; q4 < RCT_T4_NQ; ++q4)
                rct_t4_quad(rp, 16 * q4, qrow + q4 * 128);
#endif
#if RCT_T4_LH == 1
            __m256i a01a = _mm256_setzero_si256(), a23a = _mm256_setzero_si256();
            __m256i a01b = _mm256_setzero_si256(), a23b = _mm256_setzero_si256();
            for (int njq = 0; njq < RCT_T4_NCB; njq += 4) {
#if RCT_TILE4_PREBAKE
                const uint8_t *qb = rct_t4_qall
                    + (((size_t)mi * RCT_T4_NB + tb) * RCT_T4_NQ + (njq >> 2)) * 128;
#else
                const uint8_t *qb = qrow + (njq >> 2) * 128;
#endif
                const uint8_t *wb = &t4_wpk[(size_t)njq * 64];
                for (int cc = 0; cc < 4; cc += 2) {
                    const int cb = (cc >> 1) * 64;
                    __m256i wv0 = _mm256_load_si256((const __m256i *)(wb + cc * 64));
                    __m256i wv1 = _mm256_load_si256((const __m256i *)(wb + cc * 64 + 32));
                    __m256i xv0 = _mm256_load_si256((const __m256i *)(wb + cc * 64 + 64));
                    __m256i xv1 = _mm256_load_si256((const __m256i *)(wb + cc * 64 + 96));
                    a01a = _mm256_xor_si256(a01a, _mm256_gf2p8affine_epi64_epi8(
                                wv0, rct_t4_bq(qb, cb + 0), 0));
                    a23a = _mm256_xor_si256(a23a, _mm256_gf2p8affine_epi64_epi8(
                                wv0, rct_t4_bq(qb, cb + 8), 0));
                    a01a = _mm256_xor_si256(a01a, _mm256_gf2p8affine_epi64_epi8(
                                wv1, rct_t4_bq(qb, cb + 32), 0));
                    a23a = _mm256_xor_si256(a23a, _mm256_gf2p8affine_epi64_epi8(
                                wv1, rct_t4_bq(qb, cb + 40), 0));
                    a01b = _mm256_xor_si256(a01b, _mm256_gf2p8affine_epi64_epi8(
                                xv0, rct_t4_bq(qb, cb + 16), 0));
                    a23b = _mm256_xor_si256(a23b, _mm256_gf2p8affine_epi64_epi8(
                                xv0, rct_t4_bq(qb, cb + 24), 0));
                    a01b = _mm256_xor_si256(a01b, _mm256_gf2p8affine_epi64_epi8(
                                xv1, rct_t4_bq(qb, cb + 48), 0));
                    a23b = _mm256_xor_si256(a23b, _mm256_gf2p8affine_epi64_epi8(
                                xv1, rct_t4_bq(qb, cb + 56), 0));
                }
            }
            _mm256_store_si256((__m256i *)&st0p_t[((size_t)tb * 2 + 0) * 32],
                               _mm256_xor_si256(a01a, a01b));
            _mm256_store_si256((__m256i *)&st0p_t[((size_t)tb * 2 + 1) * 32],
                               _mm256_xor_si256(a23a, a23b));
#else
            __m256i a01a[RCT_T4_LH], a23a[RCT_T4_LH], a01b[RCT_T4_LH], a23b[RCT_T4_LH];
            for (int h = 0; h < RCT_T4_LH; ++h) {
                a01a[h] = _mm256_setzero_si256(); a23a[h] = _mm256_setzero_si256();
                a01b[h] = _mm256_setzero_si256(); a23b[h] = _mm256_setzero_si256();
            }
            for (int njq = 0; njq < RCT_T4_NCB; njq += 4) {
#if RCT_TILE4_PREBAKE
                const uint8_t *qb = rct_t4_qall
                    + (((size_t)mi * RCT_T4_NB + tb) * RCT_T4_NQ + (njq >> 2)) * 128;
#else
                const uint8_t *qb = qrow + (njq >> 2) * 128;
#endif
                const uint8_t *wb = &t4_wpk[(size_t)njq * (2 * RCT_T4_LH * 32)];
                for (int cc = 0; cc < 4; cc += 2) {
                    const int cb = (cc >> 1) * 64;
                    for (int h = 0; h < RCT_T4_LH; ++h) {
                        __m256i wv0 = _mm256_load_si256((const __m256i *)(wb + ((size_t)(cc * 2 + 0) * RCT_T4_LH + h) * 32));
                        __m256i wv1 = _mm256_load_si256((const __m256i *)(wb + ((size_t)(cc * 2 + 1) * RCT_T4_LH + h) * 32));
                        __m256i xv0 = _mm256_load_si256((const __m256i *)(wb + ((size_t)(cc * 2 + 2) * RCT_T4_LH + h) * 32));
                        __m256i xv1 = _mm256_load_si256((const __m256i *)(wb + ((size_t)(cc * 2 + 3) * RCT_T4_LH + h) * 32));
                        a01a[h] = _mm256_xor_si256(a01a[h], _mm256_gf2p8affine_epi64_epi8(
                                    wv0, rct_t4_bq(qb, cb + 0), 0));
                        a23a[h] = _mm256_xor_si256(a23a[h], _mm256_gf2p8affine_epi64_epi8(
                                    wv0, rct_t4_bq(qb, cb + 8), 0));
                        a01a[h] = _mm256_xor_si256(a01a[h], _mm256_gf2p8affine_epi64_epi8(
                                    wv1, rct_t4_bq(qb, cb + 32), 0));
                        a23a[h] = _mm256_xor_si256(a23a[h], _mm256_gf2p8affine_epi64_epi8(
                                    wv1, rct_t4_bq(qb, cb + 40), 0));
                        a01b[h] = _mm256_xor_si256(a01b[h], _mm256_gf2p8affine_epi64_epi8(
                                    xv0, rct_t4_bq(qb, cb + 16), 0));
                        a23b[h] = _mm256_xor_si256(a23b[h], _mm256_gf2p8affine_epi64_epi8(
                                    xv0, rct_t4_bq(qb, cb + 24), 0));
                        a01b[h] = _mm256_xor_si256(a01b[h], _mm256_gf2p8affine_epi64_epi8(
                                    xv1, rct_t4_bq(qb, cb + 48), 0));
                        a23b[h] = _mm256_xor_si256(a23b[h], _mm256_gf2p8affine_epi64_epi8(
                                    xv1, rct_t4_bq(qb, cb + 56), 0));
                    }
                }
            }
            for (int h = 0; h < RCT_T4_LH; ++h) {
                _mm256_store_si256((__m256i *)&st0p_t[(((size_t)tb * 2 + 0) * RCT_T4_LH + h) * 32],
                                   _mm256_xor_si256(a01a[h], a01b[h]));
                _mm256_store_si256((__m256i *)&st0p_t[(((size_t)tb * 2 + 1) * RCT_T4_LH + h) * 32],
                                   _mm256_xor_si256(a23a[h], a23b[h]));
            }
#endif
        }

        gf_t *st1_mi = &sum_t1[(size_t)mi * SNOVA_l2 * SNOVA_r2];
        rct_t4_left(st0p_t, t4_lqt, st1_mi, 0, (RCT_T4_NP < 5 ? RCT_T4_NP : 5));
#if RCT_T4_NP > 5
        rct_t4_left(st0p_t, t4_lqt, st1_mi, 5, (RCT_T4_NP - 5 < 5 ? RCT_T4_NP - 5 : 5));
#endif
#if RCT_T4_NP > 10
        rct_t4_left(st0p_t, t4_lqt, st1_mi, 10, (RCT_T4_NP - 10 < 5 ? RCT_T4_NP - 10 : 5));
#endif
#if RCT_T4_NP > 15
        rct_t4_left(st0p_t, t4_lqt, st1_mi, 15, (RCT_T4_NP - 15 < 5 ? RCT_T4_NP - 15 : 5));
#endif
#if RCT_T4_NP > 20
        rct_t4_left(st0p_t, t4_lqt, st1_mi, 20, (RCT_T4_NP - 20 < 5 ? RCT_T4_NP - 20 : 5));
#endif
#if RCT_T4_NP > 25
        rct_t4_left(st0p_t, t4_lqt, st1_mi, 25, (RCT_T4_NP - 25 < 5 ? RCT_T4_NP - 25 : 5));
#endif
#if RCT_T4_NP > 30
        rct_t4_left(st0p_t, t4_lqt, st1_mi, 30, (RCT_T4_NP - 30 < 5 ? RCT_T4_NP - 30 : 5));
#endif
    }

#ifdef RCT_T4_SELFTEST
    {
        RCT_SCRATCH gf_t t4_ref[SNOVA_m1 * SNOVA_l2 * SNOVA_r2];
        memcpy(t4_ref, c->sum_t1, sizeof(t4_ref));
        rct_vf_jog(c);
        size_t t4_mm = 0;
        for (size_t i = 0; i < sizeof(t4_ref); ++i)
            if (t4_ref[i] != c->sum_t1[i]) ++t4_mm;
        printf("[RCT_T4_SELFTEST] sum_t1 mismatches=%zu / %zu\n", t4_mm, sizeof(t4_ref));
    }
#else
    (void)rct_vf_jog;
#endif
}

#endif

#endif
