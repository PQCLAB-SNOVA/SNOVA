#ifndef RCT_SIGN_CM_H
#define RCT_SIGN_CM_H

#if defined(RCT_CM_ACTIVE) || defined(RCT_CMS3_ONLY)
_Static_assert(SNOVA_R >= SNOVA_L, "rct_sign_cm: requires r >= l (Bm scatter grid)");
static void rct_sign_wF_gauss_cm(rct_sign_ctx *c) {
    const gf_t *F21 = c->F21, *F12 = c->F12;
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm, *Q1 = c->Q1, *Q2 = c->Q2;
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
    uint16_t *rct_cm_whip = c->cm->cm_whip;
#if RCT_USE_GFNI
    uint8_t *rct_cm_whipb = c->cm->cm_whipb;
    (void)rct_cm_whipb;
#endif
    const uint16_t *rct_cm_Amx = c->cm->Amx, *rct_cm_Bmx = c->cm->Bmx,
                   *rct_cm_Q1x = c->cm->Q1x, *rct_cm_Q2x = c->cm->Q2x,
                   *rct_cm_q1x = c->cm->q1x, *rct_cm_q2x = c->cm->q2x;
    (void)F21; (void)F12; (void)q1; (void)q2; (void)Am; (void)Bm; (void)Q1; (void)Q2;
    (void)rct_cm_whip; (void)rct_cm_Amx; (void)rct_cm_Bmx; (void)rct_cm_Q1x;
    (void)rct_cm_Q2x; (void)rct_cm_q1x; (void)rct_cm_q2x;
        {
            RCT_SCRATCH _Alignas(32) uint16_t wF21w[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr32];
            RCT_SCRATCH _Alignas(32) uint16_t wF12w[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr32];
            for (int mi = 0; mi < SNOVA_m1; ++mi)
                for (int idx = 0; idx < SNOVA_o; ++idx)
                    for (int i1 = 0; i1 < SNOVA_l; ++i1) {
#if RCT_USE_GFNI
                        __m128i a21[RCT_CMW16], a12[RCT_CMW16];
                        for (int c = 0; c < RCT_CMW16; ++c) {
                            a21[c] = _mm_setzero_si128();
                            a12[c] = _mm_setzero_si128();
                        }
                        for (int nj = 0; nj < SNOVA_v; ++nj)
                            for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                const uint8_t *wb = &rct_cm_whipb[(k1 * SNOVA_v + nj) * RCT_CMW];
                                __m128i f21 = _mm_set1_epi8((char)
                                    F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
                                __m128i f12 = _mm_set1_epi8((char)
                                    F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1]);
                                for (int c = 0; c < RCT_CMW16; ++c) {
                                    __m128i wv = _mm_loadu_si128((const __m128i *)(wb + c * 16));
                                    a21[c] = _mm_xor_si128(a21[c], RCT_GFMUL128(f21, wv));
                                    a12[c] = _mm_xor_si128(a12[c], RCT_GFMUL128(f12, wv));
                                }
                            }
                        uint16_t *d21 = &wF21w[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32];
                        uint16_t *d12 = &wF12w[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32];
                        for (int c = 0; c < RCT_CMW16; ++c) {
                            _mm256_store_si256((__m256i *)(d21 + c * 16),
                                               _mm256_cvtepu8_epi16(rct_gfni_cleanup128(a21[c])));
                            _mm256_store_si256((__m256i *)(d12 + c * 16),
                                               _mm256_cvtepu8_epi16(rct_gfni_cleanup128(a12[c])));
                        }
#else
                        __m256i a21[RCT_CMW16], a12[RCT_CMW16];
                        for (int c = 0; c < RCT_CMW16; ++c) {
                            a21[c] = _mm256_setzero_si256();
                            a12[c] = _mm256_setzero_si256();
                        }
                        for (int nj = 0; nj < SNOVA_v; ++nj)
                            for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                const __m256i *w = (const __m256i *)&rct_cm_whip[(k1 * SNOVA_v + nj) * RCT_CMW];
                                __m256i f21 = _mm256_set1_epi16((short)
                                    F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
                                __m256i f12 = _mm256_set1_epi16((short)
                                    F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1]);
                                for (int c = 0; c < RCT_CMW16; ++c) {
                                    a21[c] = _mm256_xor_si256(a21[c], _mm256_mullo_epi16(f21, w[c]));
                                    a12[c] = _mm256_xor_si256(a12[c], _mm256_mullo_epi16(f12, w[c]));
                                }
                            }
                        uint16_t *d21 = &wF21w[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32];
                        uint16_t *d12 = &wF12w[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32];
                        for (int c = 0; c < RCT_CMW16; ++c) {
                            _mm256_store_si256((__m256i *)(d21 + c * 16), gf16_compress_u16x16(a21[c]));
                            _mm256_store_si256((__m256i *)(d12 + c * 16), gf16_compress_u16x16(a12[c]));
                        }
#endif
                    }
#if RCT_USE_GFNI
            RCT_SCRATCH _Alignas(32) uint8_t wF21sb[SNOVA_m1 * SNOVA_l * RCT_CM_OLR32];
            RCT_SCRATCH _Alignas(32) uint8_t wF12sb[SNOVA_m1 * SNOVA_l * RCT_CM_OLR32];
            memset(wF21sb, 0, sizeof(wF21sb));
            memset(wF12sb, 0, sizeof(wF12sb));
            for (int mi = 0; mi < SNOVA_m1; ++mi) {
                for (int i1 = 0; i1 < SNOVA_l; ++i1)
                    for (int idx = 0; idx < SNOVA_o; ++idx)
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            for (int j1 = 0; j1 < SNOVA_r; ++j1) {
                                wF21sb[(mi * SNOVA_l + b1) * RCT_CM_OLR32 + idx * SNOVA_lr + i1 * SNOVA_r + j1] =
                                    (uint8_t)wF21w[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32 + b1 * SNOVA_r + j1];
                                wF12sb[(mi * SNOVA_l + b1) * RCT_CM_OLR32 + idx * SNOVA_lr + i1 * SNOVA_r + j1] =
                                    (uint8_t)wF12w[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32 + b1 * SNOVA_r + j1];
                            }
            }
            RCT_SCRATCH _Alignas(32) uint8_t gausstmp1b[SNOVA_o * SNOVA_r2 * RCT_CM_OLR32];
            RCT_SCRATCH _Alignas(32) uint8_t gausstmp2b[SNOVA_o * SNOVA_r2 * RCT_CM_OLR32];
            memset(gausstmp1b, 0, sizeof(gausstmp1b));
            memset(gausstmp2b, 0, sizeof(gausstmp2b));
#if (SNOVA_l == 4) && (SNOVA_r == 4)
            const __m256i RCT_SAY[4] = {
                _mm256_broadcastsi128_si256(_mm_setr_epi8(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12)),
                _mm256_broadcastsi128_si256(_mm_setr_epi8(1,1,1,1,5,5,5,5,9,9,9,9,13,13,13,13)),
                _mm256_broadcastsi128_si256(_mm_setr_epi8(2,2,2,2,6,6,6,6,10,10,10,10,14,14,14,14)),
                _mm256_broadcastsi128_si256(_mm_setr_epi8(3,3,3,3,7,7,7,7,11,11,11,11,15,15,15,15)),
            };
            const __m256i RCT_SBY[4] = {
                _mm256_broadcastsi128_si256(_mm_setr_epi8(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3)),
                _mm256_broadcastsi128_si256(_mm_setr_epi8(4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7)),
                _mm256_broadcastsi128_si256(_mm_setr_epi8(8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11)),
                _mm256_broadcastsi128_si256(_mm_setr_epi8(12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15)),
            };
            const __m256i RCT_STY = _mm256_broadcastsi128_si256(
                _mm_setr_epi8(0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15));
#endif
            for (int mi = 0; mi < SNOVA_o; ++mi)
                for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
                    int mp = i_prime(mi, alpha);
                    const gf_t *q1raw = &q1[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                    const gf_t *q2raw = &q2[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                    const gf_t *Amraw = &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2];
                    const gf_t *Bmraw = &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr];
#if (SNOVA_l == 4) && (SNOVA_r == 4)
                    const __m256i Bm256 = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)Bmraw));
                    const __m256i Am256 = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)Amraw));
                    const __m256i Q1_256 = _mm256_broadcastsi128_si256(
                        _mm_loadu_si128((const __m128i *)&Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2]));
                    const __m256i Q2_256 = _mm256_broadcastsi128_si256(
                        _mm_loadu_si128((const __m128i *)&Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2]));
#else
                    const uint16_t *Amx = &rct_cm_Amx[(mi * SNOVA_alpha + alpha) * SNOVA_r2];
                    const uint16_t *Bmx = &rct_cm_Bmx[(mi * SNOVA_alpha + alpha) * SNOVA_lr];
                    const uint16_t *Q1x = &rct_cm_Q1x[(mi * SNOVA_alpha + alpha) * SNOVA_l2];
                    const uint16_t *Q2x = &rct_cm_Q2x[(mi * SNOVA_alpha + alpha) * SNOVA_l2];
#endif
                    {
#if (SNOVA_l == 4) && (SNOVA_r == 4)
                        _Alignas(32) uint8_t t2b[RCT_CM_OLR32];
                        {
                            __m256i acc[RCT_CM_OLR32N];
                            for (int c = 0; c < RCT_CM_OLR32N; ++c) acc[c] = _mm256_setzero_si256();
                            for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                                __m256i qv = _mm256_set1_epi8((char)q2raw[b1]);
                                const uint8_t *wf = &wF21sb[(mp * SNOVA_l + b1) * RCT_CM_OLR32];
                                for (int c = 0; c < RCT_CM_OLR32N; ++c)
                                    acc[c] = _mm256_xor_si256(acc[c], RCT_GFMUL256(qv,
                                        _mm256_loadu_si256((const __m256i *)(wf + c * 32))));
                            }
                            for (int c = 0; c < RCT_CM_OLR32N; ++c) {
                                __m256i T0 = rct_gfni_cleanup256(acc[c]);
                                __m256i t1c = rct_cm_cellmm256(T0, Bm256, RCT_SAY, RCT_SBY);
                                __m256i t2c = rct_cm_cellmm256(Q1_256, t1c, RCT_SAY, RCT_SBY);
                                _mm256_store_si256((__m256i *)&t2b[c * 32], t2c);
                            }
                        }
                        for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2) {
                                __m256i av = _mm256_set1_epi8((char)Amraw[ti1 * SNOVA_r + tj2]);
                                uint8_t *g = &gausstmp1b[((mi * SNOVA_r + ti1) * SNOVA_r + tj2) * RCT_CM_OLR32];
                                for (int c = 0; c < RCT_CM_OLR32N; ++c)
                                    _mm256_storeu_si256((__m256i *)(g + c * 32),
                                        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)(g + c * 32)),
                                            RCT_GFMUL256(av, _mm256_loadu_si256((const __m256i *)(t2b + c * 32)))));
                            }
#else
                        _Alignas(32) uint16_t t0[RCT_CM_OLR32];
                        {
                            __m256i acc[RCT_CM_OLR32N];
                            for (int c = 0; c < RCT_CM_OLR32N; ++c) acc[c] = _mm256_setzero_si256();
                            for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                                __m256i qv = _mm256_set1_epi8((char)q2raw[b1]);
                                const uint8_t *wf = &wF21sb[(mp * SNOVA_l + b1) * RCT_CM_OLR32];
                                for (int c = 0; c < RCT_CM_OLR32N; ++c)
                                    acc[c] = _mm256_xor_si256(acc[c], RCT_GFMUL256(qv,
                                        _mm256_loadu_si256((const __m256i *)(wf + c * 32))));
                            }
                            for (int c = 0; c < RCT_CM_OLR32N; ++c) {
                                __m256i r = rct_gfni_cleanup256(acc[c]);
                                _mm256_store_si256((__m256i *)&t0[c * 32],
                                    _mm256_cvtepu8_epi16(_mm256_castsi256_si128(r)));
                                _mm256_store_si256((__m256i *)&t0[c * 32 + 16],
                                    _mm256_cvtepu8_epi16(_mm256_extracti128_si256(r, 1)));
                            }
                        }
                        uint16_t t1[SNOVA_o * SNOVA_l2] = {0};
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int i1 = 0; i1 < SNOVA_l; ++i1)
                                for (int k1 = 0; k1 < SNOVA_r; ++k1) {
                                    uint16_t tv = t0[idx * SNOVA_lr + i1 * SNOVA_r + k1];
                                    for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                        t1[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
                                            (uint16_t)(tv * Bmx[k1 * SNOVA_l + j1]);
                                }
                        for (int i = 0; i < SNOVA_o * SNOVA_l2; ++i) t1[i] = rct_cm_cmp(t1[i]);
                        _Alignas(32) uint16_t t2[RCT_CM_OLR32] = {0};
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int i1 = 0; i1 < SNOVA_l; ++i1)
                                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                    uint16_t qv = Q1x[i1 * SNOVA_l + k1];
                                    for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                        t2[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
                                            (uint16_t)(qv * t1[idx * SNOVA_l2 + k1 * SNOVA_l + j1]);
                                }
                        for (int i = 0; i < SNOVA_o * SNOVA_l2; ++i) t2[i] = rct_cm_cmp(t2[i]);
                        _Alignas(32) uint8_t t2b[RCT_CM_OLR32];
                        for (int c = 0; c < RCT_CM_OLR32N; ++c) {
                            _mm_store_si128((__m128i *)&t2b[c * 32],
                                gf16_pack_u16_to_bytes(_mm256_load_si256((const __m256i *)&t2[c * 32])));
                            _mm_store_si128((__m128i *)&t2b[c * 32 + 16],
                                gf16_pack_u16_to_bytes(_mm256_load_si256((const __m256i *)&t2[c * 32 + 16])));
                        }
                        for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2) {
                                __m256i av = _mm256_set1_epi8((char)Amraw[ti1 * SNOVA_r + tj2]);
                                uint8_t *g = &gausstmp1b[((mi * SNOVA_r + ti1) * SNOVA_r + tj2) * RCT_CM_OLR32];
                                for (int c = 0; c < RCT_CM_OLR32N; ++c)
                                    _mm256_storeu_si256((__m256i *)(g + c * 32),
                                        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)(g + c * 32)),
                                            RCT_GFMUL256(av, _mm256_loadu_si256((const __m256i *)(t2b + c * 32)))));
                            }
#endif
                    }
                    {
#if (SNOVA_l == 4) && (SNOVA_r == 4)
                        _Alignas(32) uint8_t t2b[RCT_CM_OLR32];
                        {
                            __m256i acc[RCT_CM_OLR32N];
                            for (int c = 0; c < RCT_CM_OLR32N; ++c) acc[c] = _mm256_setzero_si256();
                            for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                                __m256i qv = _mm256_set1_epi8((char)q1raw[b1]);
                                const uint8_t *wf = &wF12sb[(mp * SNOVA_l + b1) * RCT_CM_OLR32];
                                for (int c = 0; c < RCT_CM_OLR32N; ++c)
                                    acc[c] = _mm256_xor_si256(acc[c], RCT_GFMUL256(qv,
                                        _mm256_loadu_si256((const __m256i *)(wf + c * 32))));
                            }
                            for (int c = 0; c < RCT_CM_OLR32N; ++c) {
                                __m256i T0 = rct_gfni_cleanup256(acc[c]);
                                __m256i T0t = _mm256_shuffle_epi8(T0, RCT_STY);
                                __m256i t1c = rct_cm_cellmm256(Am256, T0t, RCT_SAY, RCT_SBY);
                                __m256i t2c = rct_cm_cellmm256(t1c, Q2_256, RCT_SAY, RCT_SBY);
                                _mm256_store_si256((__m256i *)&t2b[c * 32], t2c);
                            }
                        }
                        for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2) {
                                __m256i bv = _mm256_set1_epi8((char)Bmraw[tj2 * SNOVA_l + ti2]);
                                uint8_t *g = &gausstmp2b[((mi * SNOVA_r + ti2) * SNOVA_r + tj2) * RCT_CM_OLR32];
                                for (int c = 0; c < RCT_CM_OLR32N; ++c)
                                    _mm256_storeu_si256((__m256i *)(g + c * 32),
                                        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)(g + c * 32)),
                                            RCT_GFMUL256(bv, _mm256_loadu_si256((const __m256i *)(t2b + c * 32)))));
                            }
#else
                        _Alignas(32) uint16_t t0[RCT_CM_OLR32];
                        {
                            __m256i acc[RCT_CM_OLR32N];
                            for (int c = 0; c < RCT_CM_OLR32N; ++c) acc[c] = _mm256_setzero_si256();
                            for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                                __m256i qv = _mm256_set1_epi8((char)q1raw[b1]);
                                const uint8_t *wf = &wF12sb[(mp * SNOVA_l + b1) * RCT_CM_OLR32];
                                for (int c = 0; c < RCT_CM_OLR32N; ++c)
                                    acc[c] = _mm256_xor_si256(acc[c], RCT_GFMUL256(qv,
                                        _mm256_loadu_si256((const __m256i *)(wf + c * 32))));
                            }
                            for (int c = 0; c < RCT_CM_OLR32N; ++c) {
                                __m256i r = rct_gfni_cleanup256(acc[c]);
                                _mm256_store_si256((__m256i *)&t0[c * 32],
                                    _mm256_cvtepu8_epi16(_mm256_castsi256_si128(r)));
                                _mm256_store_si256((__m256i *)&t0[c * 32 + 16],
                                    _mm256_cvtepu8_epi16(_mm256_extracti128_si256(r, 1)));
                            }
                        }
                        _Alignas(32) uint16_t t1[RCT_CM_OLR32] = {0};
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int i1 = 0; i1 < SNOVA_r; ++i1)
                                for (int k1 = 0; k1 < SNOVA_r; ++k1) {
                                    uint16_t av = Amx[i1 * SNOVA_r + k1];
                                    for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                        t1[idx * SNOVA_lr + i1 * SNOVA_l + j1] ^=
                                            (uint16_t)(av * t0[idx * SNOVA_lr + j1 * SNOVA_r + k1]);
                                }
                        for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) t1[i] = rct_cm_cmp(t1[i]);
                        _Alignas(32) uint16_t t2[RCT_CM_OLR32] = {0};
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int i1 = 0; i1 < SNOVA_r; ++i1)
                                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                    uint16_t tv = t1[idx * SNOVA_lr + i1 * SNOVA_l + k1];
                                    for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                        t2[idx * SNOVA_lr + i1 * SNOVA_l + j1] ^=
                                            (uint16_t)(tv * Q2x[k1 * SNOVA_l + j1]);
                                }
                        for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) t2[i] = rct_cm_cmp(t2[i]);
                        _Alignas(32) uint8_t t2b[RCT_CM_OLR32];
                        for (int c = 0; c < RCT_CM_OLR32N; ++c) {
                            _mm_store_si128((__m128i *)&t2b[c * 32],
                                gf16_pack_u16_to_bytes(_mm256_load_si256((const __m256i *)&t2[c * 32])));
                            _mm_store_si128((__m128i *)&t2b[c * 32 + 16],
                                gf16_pack_u16_to_bytes(_mm256_load_si256((const __m256i *)&t2[c * 32 + 16])));
                        }
                        for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2) {
                                __m256i bv = _mm256_set1_epi8((char)Bmraw[tj2 * SNOVA_l + ti2]);
                                uint8_t *g = &gausstmp2b[((mi * SNOVA_r + ti2) * SNOVA_r + tj2) * RCT_CM_OLR32];
                                for (int c = 0; c < RCT_CM_OLR32N; ++c)
                                    _mm256_storeu_si256((__m256i *)(g + c * 32),
                                        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)(g + c * 32)),
                                            RCT_GFMUL256(bv, _mm256_loadu_si256((const __m256i *)(t2b + c * 32)))));
                            }
#endif
                    }
                }
            for (int i = 0; i < SNOVA_o * SNOVA_r2 * RCT_CM_OLR32; i += 32) {
                _mm256_store_si256((__m256i *)&gausstmp1b[i],
                    rct_gfni_cleanup256(_mm256_load_si256((const __m256i *)&gausstmp1b[i])));
                _mm256_store_si256((__m256i *)&gausstmp2b[i],
                    rct_gfni_cleanup256(_mm256_load_si256((const __m256i *)&gausstmp2b[i])));
            }
            for (int mi = 0; mi < SNOVA_o; ++mi)
                for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                    for (int ti2 = 0; ti2 < SNOVA_l; ++ti2) {
                        gf_t *grow = gauss[mi * SNOVA_lr + ti1 * SNOVA_l + ti2];
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int tj1 = 0; tj1 < SNOVA_l; ++tj1)
                                for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                                    grow[idx * SNOVA_lr + tj1 * SNOVA_r + tj2] ^= (gf_t)(
                                        gausstmp1b[((mi * SNOVA_r + ti1) * SNOVA_r + tj2) * RCT_CM_OLR32 +
                                                   idx * SNOVA_l2 + tj1 * SNOVA_l + ti2] ^
                                        gausstmp2b[((mi * SNOVA_r + ti2) * SNOVA_r + tj2) * RCT_CM_OLR32 +
                                                   idx * SNOVA_lr + ti1 * SNOVA_l + tj1]);
                    }
            SNOVA_CLEAR_OBJ(wF21sb);
            SNOVA_CLEAR_OBJ(wF12sb);
            SNOVA_CLEAR_OBJ(gausstmp1b);
            SNOVA_CLEAR_OBJ(gausstmp2b);
#else
            RCT_SCRATCH _Alignas(32) uint16_t wF21s[SNOVA_m1 * SNOVA_l * RCT_CM_OLR];
            RCT_SCRATCH _Alignas(32) uint16_t wF12s[SNOVA_m1 * SNOVA_l * RCT_CM_OLR];
            memset(wF21s, 0, sizeof(wF21s));
            memset(wF12s, 0, sizeof(wF12s));
            for (int mi = 0; mi < SNOVA_m1; ++mi)
                for (int i1 = 0; i1 < SNOVA_l; ++i1)
                    for (int idx = 0; idx < SNOVA_o; ++idx) {
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            for (int j1 = 0; j1 < SNOVA_r; ++j1) {
                                wF21s[(mi * SNOVA_l + b1) * RCT_CM_OLR + idx * SNOVA_lr + i1 * SNOVA_r + j1] =
                                    wF21w[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32 + b1 * SNOVA_r + j1];
                                wF12s[(mi * SNOVA_l + b1) * RCT_CM_OLR + idx * SNOVA_lr + i1 * SNOVA_r + j1] =
                                    wF12w[((mi * SNOVA_l + i1) * SNOVA_o + idx) * SNOVA_lr32 + b1 * SNOVA_r + j1];
                            }
            }
            RCT_SCRATCH _Alignas(32) uint16_t gausstmp1[SNOVA_o * SNOVA_r2 * RCT_CM_OLR];
            RCT_SCRATCH _Alignas(32) uint16_t gausstmp2[SNOVA_o * SNOVA_r2 * RCT_CM_OLR];
            memset(gausstmp1, 0, sizeof(gausstmp1));
            memset(gausstmp2, 0, sizeof(gausstmp2));
            for (int mi = 0; mi < SNOVA_o; ++mi)
                for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
                    int mp = i_prime(mi, alpha);
                    const uint16_t *q1r = &rct_cm_q1x[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                    const uint16_t *q2r = &rct_cm_q2x[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                    const uint16_t *Amr = &rct_cm_Amx[(mi * SNOVA_alpha + alpha) * SNOVA_r2];
                    const uint16_t *Bmr = &rct_cm_Bmx[(mi * SNOVA_alpha + alpha) * SNOVA_lr];
                    const uint16_t *Q1r = &rct_cm_Q1x[(mi * SNOVA_alpha + alpha) * SNOVA_l2];
                    const uint16_t *Q2r = &rct_cm_Q2x[(mi * SNOVA_alpha + alpha) * SNOVA_l2];
                    {
                        _Alignas(32) uint16_t t0[RCT_CM_OLR];
                        {
                            __m256i acc[RCT_CM_OLR16];
                            for (int c = 0; c < RCT_CM_OLR16; ++c) acc[c] = _mm256_setzero_si256();
                            for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                                __m256i qv = _mm256_set1_epi16((short)q2r[b1]);
                                const __m256i *wf = (const __m256i *)&wF21s[(mp * SNOVA_l + b1) * RCT_CM_OLR];
                                for (int c = 0; c < RCT_CM_OLR16; ++c)
                                    acc[c] = _mm256_xor_si256(acc[c], _mm256_mullo_epi16(qv, wf[c]));
                            }
                            for (int c = 0; c < RCT_CM_OLR16; ++c)
                                _mm256_store_si256((__m256i *)&t0[c * 16], gf16_compress_u16x16(acc[c]));
                        }
                        uint16_t t1[SNOVA_o * SNOVA_l2] = {0};
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int i1 = 0; i1 < SNOVA_l; ++i1)
                                for (int k1 = 0; k1 < SNOVA_r; ++k1) {
                                    uint16_t tv = t0[idx * SNOVA_lr + i1 * SNOVA_r + k1];
                                    for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                        t1[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
                                            (uint16_t)(tv * Bmr[k1 * SNOVA_l + j1]);
                                }
                        for (int i = 0; i < SNOVA_o * SNOVA_l2; ++i) t1[i] = rct_cm_cmp(t1[i]);
                        _Alignas(32) uint16_t t2[RCT_CM_OLR] = {0};
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int i1 = 0; i1 < SNOVA_l; ++i1)
                                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                    uint16_t qv = Q1r[i1 * SNOVA_l + k1];
                                    for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                        t2[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ^=
                                            (uint16_t)(qv * t1[idx * SNOVA_l2 + k1 * SNOVA_l + j1]);
                                }
                        for (int i = 0; i < SNOVA_o * SNOVA_l2; ++i) t2[i] = rct_cm_cmp(t2[i]);
                        for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2) {
                                __m256i av = _mm256_set1_epi16((short)Amr[ti1 * SNOVA_r + tj2]);
                                __m256i *g = (__m256i *)&gausstmp1[((mi * SNOVA_r + ti1) * SNOVA_r + tj2) * RCT_CM_OLR];
                                const __m256i *tv = (const __m256i *)t2;
                                for (int c = 0; c < RCT_CM_OLR16; ++c)
                                    g[c] = _mm256_xor_si256(g[c], _mm256_mullo_epi16(av, tv[c]));
                            }
                    }
                    {
                        _Alignas(32) uint16_t t0[RCT_CM_OLR];
                        {
                            __m256i acc[RCT_CM_OLR16];
                            for (int c = 0; c < RCT_CM_OLR16; ++c) acc[c] = _mm256_setzero_si256();
                            for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                                __m256i qv = _mm256_set1_epi16((short)q1r[b1]);
                                const __m256i *wf = (const __m256i *)&wF12s[(mp * SNOVA_l + b1) * RCT_CM_OLR];
                                for (int c = 0; c < RCT_CM_OLR16; ++c)
                                    acc[c] = _mm256_xor_si256(acc[c], _mm256_mullo_epi16(qv, wf[c]));
                            }
                            for (int c = 0; c < RCT_CM_OLR16; ++c)
                                _mm256_store_si256((__m256i *)&t0[c * 16], gf16_compress_u16x16(acc[c]));
                        }
                        _Alignas(32) uint16_t t1[RCT_CM_OLR] = {0};
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int i1 = 0; i1 < SNOVA_r; ++i1)
                                for (int k1 = 0; k1 < SNOVA_r; ++k1) {
                                    uint16_t av = Amr[i1 * SNOVA_r + k1];
                                    for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                        t1[idx * SNOVA_lr + i1 * SNOVA_l + j1] ^=
                                            (uint16_t)(av * t0[idx * SNOVA_lr + j1 * SNOVA_r + k1]);
                                }
                        for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) t1[i] = rct_cm_cmp(t1[i]);
                        _Alignas(32) uint16_t t2[RCT_CM_OLR] = {0};
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int i1 = 0; i1 < SNOVA_r; ++i1)
                                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                    uint16_t tv = t1[idx * SNOVA_lr + i1 * SNOVA_l + k1];
                                    for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                        t2[idx * SNOVA_lr + i1 * SNOVA_l + j1] ^=
                                            (uint16_t)(tv * Q2r[k1 * SNOVA_l + j1]);
                                }
                        for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) t2[i] = rct_cm_cmp(t2[i]);
                        for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2) {
                                __m256i bv = _mm256_set1_epi16((short)Bmr[tj2 * SNOVA_l + ti2]);
                                __m256i *g = (__m256i *)&gausstmp2[((mi * SNOVA_r + ti2) * SNOVA_r + tj2) * RCT_CM_OLR];
                                const __m256i *tv = (const __m256i *)t2;
                                for (int c = 0; c < RCT_CM_OLR16; ++c)
                                    g[c] = _mm256_xor_si256(g[c], _mm256_mullo_epi16(bv, tv[c]));
                            }
                    }
                }
            for (int i = 0; i < SNOVA_o * SNOVA_r2 * RCT_CM_OLR; i += 16) {
                _mm256_store_si256((__m256i *)&gausstmp1[i],
                    gf16_compress_u16x16(_mm256_load_si256((const __m256i *)&gausstmp1[i])));
                _mm256_store_si256((__m256i *)&gausstmp2[i],
                    gf16_compress_u16x16(_mm256_load_si256((const __m256i *)&gausstmp2[i])));
            }
            for (int mi = 0; mi < SNOVA_o; ++mi)
                for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                    for (int ti2 = 0; ti2 < SNOVA_l; ++ti2) {
                        gf_t *grow = gauss[mi * SNOVA_lr + ti1 * SNOVA_l + ti2];
                        for (int idx = 0; idx < SNOVA_o; ++idx)
                            for (int tj1 = 0; tj1 < SNOVA_l; ++tj1)
                                for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                                    grow[idx * SNOVA_lr + tj1 * SNOVA_r + tj2] ^= (gf_t)(
                                        gausstmp1[((mi * SNOVA_r + ti1) * SNOVA_r + tj2) * RCT_CM_OLR +
                                                  idx * SNOVA_l2 + tj1 * SNOVA_l + ti2] ^
                                        gausstmp2[((mi * SNOVA_r + ti2) * SNOVA_r + tj2) * RCT_CM_OLR +
                                                  idx * SNOVA_lr + ti1 * SNOVA_l + tj1]);
                    }
            SNOVA_CLEAR_OBJ(wF21s);
            SNOVA_CLEAR_OBJ(wF12s);
            SNOVA_CLEAR_OBJ(gausstmp1);
            SNOVA_CLEAR_OBJ(gausstmp2);
#endif
            SNOVA_CLEAR_OBJ(wF21w);
            SNOVA_CLEAR_OBJ(wF12w);
        }
}
#endif

#endif
