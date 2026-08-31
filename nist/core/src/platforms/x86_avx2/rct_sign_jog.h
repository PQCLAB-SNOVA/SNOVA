#ifndef RCT_SIGN_JOG_H
#define RCT_SIGN_JOG_H

#if RCT_SIGN_JOG

#ifndef RCT_SIGN_P2_B1LANE
#define RCT_SIGN_P2_B1LANE 1
#endif
#if RCT_SIGN_P2_B1LANE && (RCT_SJ_M4 || (RCT_SJ_A4 && SNOVA_r != 8))
#define RCT_SJ_B1L 1
#else
#define RCT_SJ_B1L 0
#endif
#if RCT_SJ_M4 && RCT_SJ_B1L
enum { RCT_B1L_WN = SNOVA_l * SNOVA_v * SNOVA_lr,
       RCT_B1L_WCH = (RCT_B1L_WN + 15) / 16 * 16,
       RCT_B1L_FN = SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2,
       RCT_B1L_FCH = (RCT_B1L_FN + 15) / 16 * 16 };
#endif

static void rct_sign_whipbuild_jog(rct_sign_ctx *c) {
    rct_sj_ensure();
    const gf_t *sig = c->signature_in_GF;
    gf_t *whipped_sig = c->whipped_sig;
    for (int ab = 0; ab < SNOVA_l; ++ab)
        for (int ni = 0; ni < SNOVA_v; ++ni)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                __m128i acc = _mm_setzero_si128();
                const gf_t *Srow = &rct_S[ab * SNOVA_l2 + i1 * SNOVA_l];
                for (int k1 = 0; k1 < SNOVA_l; ++k1)
                    acc = _mm_xor_si128(acc, rct_sj_sv128(rct_sj_bc128_pub(Srow[k1]),
                        _mm_loadu_si128((const __m128i *)&sig[ni * SNOVA_lr + k1 * SNOVA_r])));
#if RCT_SJ_B1L
                rct_sj_store_r(&whipped_sig[(ni * SNOVA_l + i1) * SNOVA_lr + ab * SNOVA_r], acc);
#else
                rct_sj_store_r(&whipped_sig[(ab * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r], acc);
#endif
            }
}

static void rct_sign_sumt_jog(rct_sign_ctx *c) {
    rct_sj_ensure();
    const gf_t *P11 = c->P11, *whipped_sig = c->whipped_sig;
    gf_t *sum_t1 = c->sum_t1;
    RCT_SCRATCH _Alignas(32) gf_t sum_t0[SNOVA_m1 * SNOVA_l * SNOVA_v * SNOVA_lr + 16];
#if (RCT_SJ_A4 || RCT_SJ_M4) && !RCT_SJ_B1L
    RCT_SCRATCH _Alignas(32) gf_t rct_a4_wsigT[SNOVA_l * SNOVA_v * SNOVA_lr + 16];
#endif
    memset(sum_t0 + SNOVA_m1 * SNOVA_l * SNOVA_v * SNOVA_lr, 0, 16);
#if (RCT_SJ_A4 || RCT_SJ_M4) && !RCT_SJ_B1L
    memset(rct_a4_wsigT + SNOVA_l * SNOVA_v * SNOVA_lr, 0, 16);
#endif
#if RCT_SJ_B1L
#if RCT_SJ_M4
    RCT_SCRATCH _Alignas(32) uint16_t whipe16[RCT_B1L_WCH + 16];
    RCT_SCRATCH _Alignas(32) uint16_t st0r16[RCT_B1L_WCH + 16];
    memset(&whipe16[RCT_B1L_WCH], 0, 16 * sizeof(uint16_t));
    memset(&st0r16[RCT_B1L_WCH], 0, 16 * sizeof(uint16_t));
    rct_m4_expand_buf(whipe16, whipped_sig, RCT_B1L_WN);
#endif
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        for (int ni = 0; ni < SNOVA_v; ++ni)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                const gf_t *pr = &P11[((mi * SNOVA_v + ni) * SNOVA_v) * SNOVA_l2 + i1 * SNOVA_l];
#if SNOVA_r == 8
                __m256i ae = _mm256_setzero_si256(), ao = _mm256_setzero_si256();
                __m128i te = _mm_setzero_si128(), to = _mm_setzero_si128();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const gf_t *wb = &whipped_sig[(nj * SNOVA_l) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m128i tb = rct_sj_bc128_pub(pr[nj * SNOVA_l2 + k1]);
                        __m256i pd = rct_sj_sv256(_mm256_broadcastsi128_si256(tb),
                            _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr]));
                        __m128i pt = rct_sj_sv128(tb,
                            _mm_loadu_si128((const __m128i *)&wb[k1 * SNOVA_lr + 32]));
                        if (k1 & 1) { ao = _mm256_xor_si256(ao, pd); to = _mm_xor_si128(to, pt); }
                        else        { ae = _mm256_xor_si256(ae, pd); te = _mm_xor_si128(te, pt); }
                    }
                }
                _Alignas(32) uint8_t tp[48];
                _mm256_store_si256((__m256i *)tp, RCT_SJ_CLEAN256(_mm256_xor_si256(ae, ao)));
                _mm_store_si128((__m128i *)&tp[32], RCT_SJ_CLEAN(_mm_xor_si128(te, to)));
                memcpy(&sum_t0[((mi * SNOVA_v + ni) * SNOVA_l + i1) * SNOVA_lr], tp, SNOVA_lr);
#else
                __m256i ae = _mm256_setzero_si256(), ao = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const gf_t *wb = &whipped_sig[(nj * SNOVA_l) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i pd = rct_sj_sv256(rct_sj_bc256_pub(pr[nj * SNOVA_l2 + k1]),
                            _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr]));
                        if (k1 & 1) ao = _mm256_xor_si256(ao, pd);
                        else        ae = _mm256_xor_si256(ae, pd);
                    }
                }
                _Alignas(32) uint8_t tp[32];
                _mm256_store_si256((__m256i *)tp, RCT_SJ_CLEAN256(_mm256_xor_si256(ae, ao)));
                memcpy(&sum_t0[((mi * SNOVA_v + ni) * SNOVA_l + i1) * SNOVA_lr], tp, SNOVA_lr);
#endif
            }
#if RCT_SJ_M4
        rct_m4_raw_buf(st0r16, &sum_t0[mi * SNOVA_v * SNOVA_l * SNOVA_lr], RCT_B1L_WN);
#endif
        for (int a1 = 0; a1 < SNOVA_l; ++a1)
            for (int i1 = 0; i1 < SNOVA_r; ++i1) {
#if RCT_SJ_M4
#if SNOVA_r == 8
                __m256i ae0 = _mm256_setzero_si256(), ae1 = _mm256_setzero_si256(), ae2 = _mm256_setzero_si256();
                __m256i ao0 = _mm256_setzero_si256(), ao1 = _mm256_setzero_si256(), ao2 = _mm256_setzero_si256();
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i s = _mm256_set1_epi16((short)whipe16[(ni * SNOVA_l + k1) * SNOVA_lr + a1 * SNOVA_r + i1]);
                        const uint16_t *sb = &st0r16[(ni * SNOVA_l + k1) * SNOVA_lr];
                        __m256i p0 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)sb));
                        __m256i p1 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&sb[16]));
                        __m256i p2 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&sb[32]));
                        if (k1 & 1) { ao0 = _mm256_xor_si256(ao0, p0); ao1 = _mm256_xor_si256(ao1, p1); ao2 = _mm256_xor_si256(ao2, p2); }
                        else        { ae0 = _mm256_xor_si256(ae0, p0); ae1 = _mm256_xor_si256(ae1, p1); ae2 = _mm256_xor_si256(ae2, p2); }
                    }
                _Alignas(32) uint8_t tp[48];
                _mm_store_si128((__m128i *)tp,
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae0, ao0))));
                _mm_store_si128((__m128i *)&tp[16],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae1, ao1))));
                _mm_store_si128((__m128i *)&tp[32],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae2, ao2))));
#else
                __m256i ae0 = _mm256_setzero_si256(), ae1 = _mm256_setzero_si256();
                __m256i ao0 = _mm256_setzero_si256(), ao1 = _mm256_setzero_si256();
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i s = _mm256_set1_epi16((short)whipe16[(ni * SNOVA_l + k1) * SNOVA_lr + a1 * SNOVA_r + i1]);
                        const uint16_t *sb = &st0r16[(ni * SNOVA_l + k1) * SNOVA_lr];
                        __m256i p0 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)sb));
                        __m256i p1 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&sb[16]));
                        if (k1 & 1) { ao0 = _mm256_xor_si256(ao0, p0); ao1 = _mm256_xor_si256(ao1, p1); }
                        else        { ae0 = _mm256_xor_si256(ae0, p0); ae1 = _mm256_xor_si256(ae1, p1); }
                    }
                _Alignas(32) uint8_t tp[32];
                _mm_store_si128((__m128i *)tp,
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae0, ao0))));
                _mm_store_si128((__m128i *)&tp[16],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae1, ao1))));
#endif
#else
                __m256i ae = _mm256_setzero_si256(), ao = _mm256_setzero_si256();
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i bc = _mm256_set1_epi8((char)whipped_sig[(ni * SNOVA_l + k1) * SNOVA_lr + a1 * SNOVA_r + i1]);
                        __m256i pd = _mm256_gf2p8mul_epi8(bc,
                            _mm256_loadu_si256((const __m256i *)&sum_t0[((mi * SNOVA_v + ni) * SNOVA_l + k1) * SNOVA_lr]));
                        if (k1 & 1) ao = _mm256_xor_si256(ao, pd);
                        else        ae = _mm256_xor_si256(ae, pd);
                    }
                _Alignas(32) uint8_t tp[32];
                _mm256_store_si256((__m256i *)tp, RCT_SJ_CLEAN256(_mm256_xor_si256(ae, ao)));
#endif
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    memcpy(&sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r],
                           &tp[b1 * SNOVA_r], SNOVA_r);
            }
    }
#if RCT_SJ_M4
    SNOVA_CLEAR_OBJ(whipe16);
    SNOVA_CLEAR_OBJ(st0r16);
#endif
#else
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
#if RCT_SJ_A4
        for (int ni = 0; ni < SNOVA_v; ++ni)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                rct_a4_acc_t acc[SNOVA_l];
                for (int b1 = 0; b1 < SNOVA_l; ++b1) acc[b1] = rct_a4_zero();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    rct_a4_bc_t bv = rct_a4_bc(&P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l]);
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        rct_a4_mac(&acc[b1], bv, &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr]);
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    rct_sj_store_r(&sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r], rct_a4_fold(acc[b1]));
            }
#else
        for (int b1 = 0; b1 < SNOVA_l; ++b1)
            for (int ni = 0; ni < SNOVA_v; ++ni)
                for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                    __m128i acc = _mm_setzero_si128();
                    for (int nj = 0; nj < SNOVA_v; ++nj) {
                        const gf_t *pc = &P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l];
                        const gf_t *w = &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr];
                        for (int k1 = 0; k1 < SNOVA_l; ++k1)
                            acc = _mm_xor_si128(acc, rct_sj_sv128(rct_sj_bc128_pub(pc[k1]),
                                _mm_loadu_si128((const __m128i *)&w[k1 * SNOVA_r])));
                    }
                    rct_sj_store_r(&sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r], acc);
                }
#endif
#if RCT_SJ_A4
        if (mi == 0) {
            for (int a1 = 0; a1 < SNOVA_l; ++a1)
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int i1 = 0; i1 < SNOVA_r; ++i1)
                        for (int k1 = 0; k1 < SNOVA_l; ++k1)
                            rct_a4_wsigT[(a1 * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_l + k1] =
                                whipped_sig[(a1 * SNOVA_v + ni) * SNOVA_lr + k1 * SNOVA_r + i1];
        }
        for (int a1 = 0; a1 < SNOVA_l; ++a1)
            for (int i1 = 0; i1 < SNOVA_r; ++i1) {
                rct_a4_acc_t acc[SNOVA_l];
                for (int b1 = 0; b1 < SNOVA_l; ++b1) acc[b1] = rct_a4_zero();
                for (int ni = 0; ni < SNOVA_v; ++ni) {
                    rct_a4_bc_t bv = rct_a4_bc(&rct_a4_wsigT[(a1 * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_l]);
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        rct_a4_mac(&acc[b1], bv, &sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr]);
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    rct_sj_store_r(&sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r], rct_a4_fold(acc[b1]));
            }
#elif RCT_SJ_M4
        if (mi == 0) {
            for (int a1 = 0; a1 < SNOVA_l; ++a1)
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int i1 = 0; i1 < SNOVA_r; ++i1)
                        for (int k1 = 0; k1 < SNOVA_l; ++k1)
                            rct_a4_wsigT[(a1 * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_l + k1] =
                                whipped_sig[(a1 * SNOVA_v + ni) * SNOVA_lr + k1 * SNOVA_r + i1];
        }
        for (int a1 = 0; a1 < SNOVA_l; ++a1)
            for (int i1 = 0; i1 < SNOVA_r; ++i1) {
                rct_m4_acc_t acc[SNOVA_l];
                for (int b1 = 0; b1 < SNOVA_l; ++b1) acc[b1] = rct_m4_zero();
                for (int ni = 0; ni < SNOVA_v; ++ni) {
                    rct_m4_bc_t bv = rct_m4_bc(&rct_a4_wsigT[(a1 * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_l]);
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        rct_m4_mac(&acc[b1], bv, &sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr]);
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    rct_sj_store_r(&sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r], rct_m4_fold(acc[b1]));
            }
#else
        for (int a1 = 0; a1 < SNOVA_l; ++a1)
            for (int i1 = 0; i1 < SNOVA_r; ++i1) {
                __m128i acc[SNOVA_l];
                for (int b1 = 0; b1 < SNOVA_l; ++b1) acc[b1] = _mm_setzero_si128();
                for (int ni = 0; ni < SNOVA_v; ++ni) {
                    const gf_t *w = &whipped_sig[(a1 * SNOVA_v + ni) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m128i bc = rct_sj_bc128(w[k1 * SNOVA_r + i1]);
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            acc[b1] = _mm_xor_si128(acc[b1], rct_sj_sv128(bc,
                                _mm_loadu_si128((const __m128i *)&sum_t0[(((mi * SNOVA_l + b1) * SNOVA_v + ni)) * SNOVA_lr + k1 * SNOVA_r])));
                    }
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    rct_sj_store_r(&sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r], acc[b1]);
            }
#endif
    }
#endif
    SNOVA_CLEAR_OBJ(sum_t0);
#if (RCT_SJ_A4 || RCT_SJ_M4) && !RCT_SJ_B1L
    SNOVA_CLEAR_OBJ(rct_a4_wsigT);
#endif
}

static void rct_sign_wf_jog(rct_sign_ctx *c) {
    rct_sj_ensure();
    const gf_t *F21 = c->F21, *F12 = c->F12, *whipped_sig = c->whipped_sig;
    gf_t *whipped_F21 = c->whipped_F21, *whipped_F12 = c->whipped_F12;
#if RCT_SJ_B1L
#if RCT_SJ_M4
    RCT_SCRATCH _Alignas(32) uint16_t whipr16[RCT_B1L_WCH + 16];
    RCT_SCRATCH _Alignas(32) uint16_t f21e16[RCT_B1L_FCH + 16];
    RCT_SCRATCH _Alignas(32) uint16_t f12e16[RCT_B1L_FCH + 16];
    memset(&whipr16[RCT_B1L_WCH], 0, 16 * sizeof(uint16_t));
    rct_m4_raw_buf(whipr16, whipped_sig, RCT_B1L_WN);
    rct_m4_expand_buf(f21e16, F21, RCT_B1L_FN);
    rct_m4_expand_buf(f12e16, F12, RCT_B1L_FN);
#endif
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        for (int idx = 0; idx < SNOVA_o; ++idx)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
#if RCT_SJ_M4
#if SNOVA_r == 8
                __m256i ae0 = _mm256_setzero_si256(), ae1 = _mm256_setzero_si256(), ae2 = _mm256_setzero_si256();
                __m256i ao0 = _mm256_setzero_si256(), ao1 = _mm256_setzero_si256(), ao2 = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const uint16_t *fr = &f21e16[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l];
                    const uint16_t *wb = &whipr16[(nj * SNOVA_l) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i s = _mm256_set1_epi16((short)fr[k1]);
                        __m256i p0 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr]));
                        __m256i p1 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr + 16]));
                        __m256i p2 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr + 32]));
                        if (k1 & 1) { ao0 = _mm256_xor_si256(ao0, p0); ao1 = _mm256_xor_si256(ao1, p1); ao2 = _mm256_xor_si256(ao2, p2); }
                        else        { ae0 = _mm256_xor_si256(ae0, p0); ae1 = _mm256_xor_si256(ae1, p1); ae2 = _mm256_xor_si256(ae2, p2); }
                    }
                }
                _Alignas(32) uint8_t tp[48];
                _mm_store_si128((__m128i *)tp,
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae0, ao0))));
                _mm_store_si128((__m128i *)&tp[16],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae1, ao1))));
                _mm_store_si128((__m128i *)&tp[32],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae2, ao2))));
#else
                __m256i ae0 = _mm256_setzero_si256(), ae1 = _mm256_setzero_si256();
                __m256i ao0 = _mm256_setzero_si256(), ao1 = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const uint16_t *fr = &f21e16[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l];
                    const uint16_t *wb = &whipr16[(nj * SNOVA_l) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i s = _mm256_set1_epi16((short)fr[k1]);
                        __m256i p0 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr]));
                        __m256i p1 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr + 16]));
                        if (k1 & 1) { ao0 = _mm256_xor_si256(ao0, p0); ao1 = _mm256_xor_si256(ao1, p1); }
                        else        { ae0 = _mm256_xor_si256(ae0, p0); ae1 = _mm256_xor_si256(ae1, p1); }
                    }
                }
                _Alignas(32) uint8_t tp[32];
                _mm_store_si128((__m128i *)tp,
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae0, ao0))));
                _mm_store_si128((__m128i *)&tp[16],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae1, ao1))));
#endif
#else
                __m256i ae = _mm256_setzero_si256(), ao = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const gf_t *fr = &F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l];
                    const gf_t *wb = &whipped_sig[(nj * SNOVA_l) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i pd = _mm256_gf2p8mul_epi8(_mm256_set1_epi8((char)fr[k1]),
                            _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr]));
                        if (k1 & 1) ao = _mm256_xor_si256(ao, pd);
                        else        ae = _mm256_xor_si256(ae, pd);
                    }
                }
                _Alignas(32) uint8_t tp[32];
                _mm256_store_si256((__m256i *)tp, RCT_SJ_CLEAN256(_mm256_xor_si256(ae, ao)));
#endif
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    memcpy(&whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r],
                           &tp[b1 * SNOVA_r], SNOVA_r);
            }
        for (int idx = 0; idx < SNOVA_o; ++idx)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
#if RCT_SJ_M4
#if SNOVA_r == 8
                __m256i ae0 = _mm256_setzero_si256(), ae1 = _mm256_setzero_si256(), ae2 = _mm256_setzero_si256();
                __m256i ao0 = _mm256_setzero_si256(), ao1 = _mm256_setzero_si256(), ao2 = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const uint16_t *fr = &f12e16[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + i1];
                    const uint16_t *wb = &whipr16[(nj * SNOVA_l) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i s = _mm256_set1_epi16((short)fr[k1 * SNOVA_l]);
                        __m256i p0 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr]));
                        __m256i p1 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr + 16]));
                        __m256i p2 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr + 32]));
                        if (k1 & 1) { ao0 = _mm256_xor_si256(ao0, p0); ao1 = _mm256_xor_si256(ao1, p1); ao2 = _mm256_xor_si256(ao2, p2); }
                        else        { ae0 = _mm256_xor_si256(ae0, p0); ae1 = _mm256_xor_si256(ae1, p1); ae2 = _mm256_xor_si256(ae2, p2); }
                    }
                }
                _Alignas(32) uint8_t tp[48];
                _mm_store_si128((__m128i *)tp,
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae0, ao0))));
                _mm_store_si128((__m128i *)&tp[16],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae1, ao1))));
                _mm_store_si128((__m128i *)&tp[32],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae2, ao2))));
#else
                __m256i ae0 = _mm256_setzero_si256(), ae1 = _mm256_setzero_si256();
                __m256i ao0 = _mm256_setzero_si256(), ao1 = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const uint16_t *fr = &f12e16[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + i1];
                    const uint16_t *wb = &whipr16[(nj * SNOVA_l) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i s = _mm256_set1_epi16((short)fr[k1 * SNOVA_l]);
                        __m256i p0 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr]));
                        __m256i p1 = _mm256_mullo_epi16(s, _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr + 16]));
                        if (k1 & 1) { ao0 = _mm256_xor_si256(ao0, p0); ao1 = _mm256_xor_si256(ao1, p1); }
                        else        { ae0 = _mm256_xor_si256(ae0, p0); ae1 = _mm256_xor_si256(ae1, p1); }
                    }
                }
                _Alignas(32) uint8_t tp[32];
                _mm_store_si128((__m128i *)tp,
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae0, ao0))));
                _mm_store_si128((__m128i *)&tp[16],
                    cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(_mm256_xor_si256(ae1, ao1))));
#endif
#else
                __m256i ae = _mm256_setzero_si256(), ao = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const gf_t *fr = &F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + i1];
                    const gf_t *wb = &whipped_sig[(nj * SNOVA_l) * SNOVA_lr];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m256i pd = _mm256_gf2p8mul_epi8(_mm256_set1_epi8((char)fr[k1 * SNOVA_l]),
                            _mm256_loadu_si256((const __m256i *)&wb[k1 * SNOVA_lr]));
                        if (k1 & 1) ao = _mm256_xor_si256(ao, pd);
                        else        ae = _mm256_xor_si256(ae, pd);
                    }
                }
                _Alignas(32) uint8_t tp[32];
                _mm256_store_si256((__m256i *)tp, RCT_SJ_CLEAN256(_mm256_xor_si256(ae, ao)));
#endif
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    memcpy(&whipped_F12[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r],
                           &tp[b1 * SNOVA_r], SNOVA_r);
            }
    }
#if RCT_SJ_M4
    SNOVA_CLEAR_OBJ(whipr16);
    SNOVA_CLEAR_OBJ(f21e16);
    SNOVA_CLEAR_OBJ(f12e16);
#endif
#else
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
#if RCT_SJ_A4
        for (int idx = 0; idx < SNOVA_o; ++idx)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                rct_a4_acc_t acc[SNOVA_l];
                for (int b1 = 0; b1 < SNOVA_l; ++b1) acc[b1] = rct_a4_zero();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    rct_a4_bc_t bv = rct_a4_bc(&F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l]);
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        rct_a4_mac(&acc[b1], bv, &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr]);
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    rct_sj_store_r(&whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r], rct_a4_fold(acc[b1]));
            }
#elif RCT_SJ_M4
        for (int idx = 0; idx < SNOVA_o; ++idx)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                rct_m4_acc_t acc[SNOVA_l];
                for (int b1 = 0; b1 < SNOVA_l; ++b1) acc[b1] = rct_m4_zero();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    rct_m4_bc_t bv = rct_m4_bc(&F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l]);
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        rct_m4_mac(&acc[b1], bv, &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr]);
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    rct_sj_store_r(&whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r], rct_m4_fold(acc[b1]));
            }
#else
        for (int idx = 0; idx < SNOVA_o; ++idx)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                __m128i acc[SNOVA_l];
                for (int b1 = 0; b1 < SNOVA_l; ++b1) acc[b1] = _mm_setzero_si128();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    const gf_t *fr = &F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l];
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m128i bc = rct_sj_bc128(fr[k1]);
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            acc[b1] = _mm_xor_si128(acc[b1], rct_sj_sv128(bc,
                                _mm_loadu_si128((const __m128i *)&whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr + k1 * SNOVA_r])));
                    }
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    rct_sj_store_r(&whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r], acc[b1]);
            }
#endif
        for (int idx = 0; idx < SNOVA_o; ++idx)
            for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                __m128i acc[SNOVA_l];
                for (int b1 = 0; b1 < SNOVA_l; ++b1) acc[b1] = _mm_setzero_si128();
                for (int nj = 0; nj < SNOVA_v; ++nj) {
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        __m128i bc = rct_sj_bc128(F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1]);
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            acc[b1] = _mm_xor_si128(acc[b1], rct_sj_sv128(bc,
                                _mm_loadu_si128((const __m128i *)&whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr + k1 * SNOVA_r])));
                    }
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    rct_sj_store_r(&whipped_F12[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r], acc[b1]);
            }
    }
#endif
}

static void rct_sign_fvv_jog(rct_sign_ctx *c) {
    rct_sj_ensure();
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm, *sum_t1 = c->sum_t1;
    gf_t *Fvv = c->Fvv;
    for (int mi = 0; mi < SNOVA_o; ++mi)
        for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
            int mi_prime = i_prime(mi, alpha);
            gf_t temp1[SNOVA_r2 + 16] = {0};
            gf_t temp2[SNOVA_lr + 16] = {0};
            const gf_t *q1r = &q1[(mi * SNOVA_alpha + alpha) * SNOVA_l];
            const gf_t *q2r = &q2[(mi * SNOVA_alpha + alpha) * SNOVA_l];
            for (int cc = 0; cc < SNOVA_r2; cc += 16) {
                __m128i t1 = _mm_setzero_si128();
                for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                    __m128i t0 = _mm_setzero_si128();
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        t0 = _mm_xor_si128(t0, rct_sj_sv128(rct_sj_bc128_pub(q2r[b1]),
                            _mm_loadu_si128((const __m128i *)
                                &sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + cc])));
                    t0 = RCT_SJ_CLEAN(t0);
                    t1 = _mm_xor_si128(t1, rct_sj_sv128(rct_sj_bc128_pub(q1r[a1]), t0)    );
                }
                _mm_storeu_si128((__m128i *)&temp1[cc], RCT_SJ_CLEAN(t1));
            }
            rct_sj_mm_add(temp2, temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_r, SNOVA_r, SNOVA_l);
            rct_sj_mm_add_pub(&Fvv[mi * SNOVA_lr], &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], temp2, SNOVA_r, SNOVA_r, SNOVA_l);
        }
}

#if SNOVA_Q == 16 && RCT_GFMUL_ANY && !defined(RCT_GAUSS_SCALAR)
static void rct_sign_apply_t12_jog(rct_sign_ctx *c, const gf_t *solpad) {
    rct_sj_ensure();
    gf_t *signature_in_GF = c->signature_in_GF;
    const gf_t *T12 = c->T12;
    for (int index = 0; index < SNOVA_v; ++index)
        for (int mi = 0; mi < SNOVA_o; ++mi)
            rct_sj_mm_add(&signature_in_GF[index * SNOVA_lr],
                          &T12[(index * SNOVA_o + mi) * SNOVA_l2],
                          &solpad[mi * SNOVA_lr], SNOVA_l, SNOVA_l, SNOVA_r);
}
#endif

static void rct_skx_fold_F_jog(gf_t *F21, gf_t *F12, const gf_t *T12, const gf_t *P11) {
    rct_sj_ensure();
#if RCT_F5_A4
    for (int k1 = 0; k1 < SNOVA_o; ++k1) {
        rct_a4f_bc_t bct[SNOVA_v][5];
        for (int j2 = 0; j2 < SNOVA_v; ++j2)
            for (int i = 0; i < 5; ++i)
                bct[j2][i] = rct_a4f_bc(&T12[(j2 * SNOVA_o + k1) * SNOVA_l2 + i * 5]);
        for (int i1 = 0; i1 < SNOVA_m1; ++i1)
            for (int j1 = 0; j1 < SNOVA_v; ++j1) {
                gf_t *C = &F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2];
                rct_a4f_acc_t a[5];
                for (int i = 0; i < 5; ++i) a[i] = rct_a4f_zero();
                for (int j2 = 0; j2 < SNOVA_v; ++j2) {
                    const gf_t *w = &P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2];
                    for (int i = 0; i < 5; ++i) rct_a4f_mac(&a[i], bct[j2][i], w);
                }
                for (int i = 0; i < 5; ++i) rct_a4f_xor5(&C[i * 5], rct_a4f_fold(a[i]));
            }
        SNOVA_CLEAR_OBJ(bct);
    }
#elif RCT_F5_M4 && !RCT_FOLD_WIDE
    for (int k1 = 0; k1 < SNOVA_o; ++k1) {
        rct_m4f_bc_t bct[SNOVA_v][5];
        for (int j2 = 0; j2 < SNOVA_v; ++j2) {
            _Alignas(32) uint16_t ev[32];
            rct_m4_expand_buf(ev, &T12[(j2 * SNOVA_o + k1) * SNOVA_l2], SNOVA_l2);
            for (int i = 0; i < 5; ++i) bct[j2][i] = rct_m4f2_bc_mirror(&ev[i * 5]);
        }
        for (int i1 = 0; i1 < SNOVA_m1; ++i1)
            for (int j1 = 0; j1 < SNOVA_v; ++j1) {
                gf_t *C = &F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2];
                __m256i a0[5], a1[5];
                for (int i = 0; i < 5; ++i) { a0[i] = _mm256_setzero_si256(); a1[i] = _mm256_setzero_si256(); }
                for (int j2 = 0; j2 < SNOVA_v; ++j2)
                    rct_m4f2_mac5(a0, a1, bct[j2],
                        &P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2]);
                rct_m4f2_fold5_xor(C, a0, a1);
            }
        SNOVA_CLEAR_OBJ(bct);
    }
#elif RCT_F5_WIDE
    {
        enum { RCT_MVL = RCT_WIDE_MVL };
        RCT_WIDE_SCRATCH_DECL;
        memset(Fw, 0, RCT_WIDE_FW_LEN * sizeof(uint16_t));
        for (int nk = 0; nk < SNOVA_v; ++nk)
            for (int k1 = 0; k1 < SNOVA_l; ++k1)
                for (int mi = 0; mi < SNOVA_m1; ++mi)
                    for (int nj = 0; nj < SNOVA_v; ++nj)
                        for (int j1 = 0; j1 < SNOVA_l; ++j1)
                            P11aw[(nk * SNOVA_l + k1) * RCT_MVL + (mi * SNOVA_v + nj) * SNOVA_l + j1] =
                                P11[((mi * SNOVA_v + nk) * SNOVA_v + nj) * SNOVA_l2 + k1 * SNOVA_l + j1];
        for (int ni = 0; ni < SNOVA_o; ++ni)
            for (int i1 = 0; i1 < SNOVA_l; ++i1)
                for (int nk = 0; nk < SNOVA_v; ++nk)
                    for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                        uint16_t s = cl_expand_scalar16(T12[(nk * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + k1]);
                        uint16_t *Fr = &Fw[(ni * SNOVA_l + i1) * RCT_MVL];
                        const uint16_t *Pr = &P11aw[(nk * SNOVA_l + k1) * RCT_MVL];
                        for (int mi = 0; mi < RCT_MVL; ++mi) Fr[mi] ^= (uint16_t)(s * Pr[mi]);
                    }
        for (int i = 0; i < RCT_MVL * SNOVA_o * SNOVA_l; i += 16)
            _mm256_storeu_si256((__m256i *)&Fw[i], cl_gf16_compress_u16x16(_mm256_loadu_si256((const __m256i *)&Fw[i])));
        for (int mi = 0; mi < SNOVA_m1; ++mi)
            for (int nj = 0; nj < SNOVA_v; ++nj)
                for (int ni = 0; ni < SNOVA_o; ++ni)
                    for (int i1 = 0; i1 < SNOVA_l; ++i1)
                        for (int j1 = 0; j1 < SNOVA_l; ++j1)
                            F21[((mi * SNOVA_o + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] =
                                (gf_t)Fw[(ni * SNOVA_l + i1) * RCT_MVL + (mi * SNOVA_v + nj) * SNOVA_l + j1];
        SNOVA_CLEAR(Fw, RCT_WIDE_FW_LEN * sizeof(uint16_t));
        rct_wide_fold_F12(F12, P11, T12, P11aw, Fw);
    }
#elif RCT_SIGN_JOG && !RCT_HAVE_GFNI
    for (int i1 = 0; i1 < SNOVA_m1; ++i1)
        for (int j2 = 0; j2 < SNOVA_v; ++j2)
            for (int k1 = 0; k1 < SNOVA_o; ++k1) {
                __m128i bcx[SNOVA_l2];
                for (int t = 0; t < SNOVA_l2; ++t)
                    bcx[t] = rct_sj_bc128(T12[(j2 * SNOVA_o + k1) * SNOVA_l2 + t]);
                for (int j1 = 0; j1 < SNOVA_v; ++j1) {
                    gf_t *C = &F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2];
                    const gf_t *B = &P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2];
                    for (int i = 0; i < SNOVA_l; ++i) {
                        __m128i acc = _mm_setzero_si128();
                        for (int k = 0; k < SNOVA_l; ++k)
                            acc = _mm_xor_si128(acc, rct_sj_sv128(bcx[i * SNOVA_l + k],
                                _mm_loadu_si128((const __m128i *)&B[k * SNOVA_l])));
                        _Alignas(16) uint8_t pb[16];
                        _mm_store_si128((__m128i *)pb, RCT_SJ_CLEAN(acc));
                        for (int j = 0; j < SNOVA_l; ++j)
                            C[i * SNOVA_l + j] = (gf_t)(C[i * SNOVA_l + j] ^ pb[j]);
                    }
                }
            }
#else
    for (int i1 = 0; i1 < SNOVA_m1; ++i1)
        for (int j1 = 0; j1 < SNOVA_v; ++j1)
            for (int j2 = 0; j2 < SNOVA_v; ++j2)
                for (int k1 = 0; k1 < SNOVA_o; ++k1)
                    rct_sj_mm_add(&F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2],
                                  &T12[(j2 * SNOVA_o + k1) * SNOVA_l2],
                                  &P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2],
                                  SNOVA_l, SNOVA_l, SNOVA_l);
#endif
#if RCT_F5_WIDE
#elif RCT_F5_A4 || RCT_F5_M4
    rct_f5_fold_F12(F12, P11, T12);
#else
    for (int i1 = 0; i1 < SNOVA_m1; ++i1)
        for (int j1 = 0; j1 < SNOVA_v; ++j1)
            for (int j2 = 0; j2 < SNOVA_v; ++j2)
                for (int k1 = 0; k1 < SNOVA_o; ++k1)
                    rct_sj_mm_add_pub(&F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
                                  &P11[((i1 * SNOVA_v + j1) * SNOVA_v + j2) * SNOVA_l2],
                                  &T12[(j2 * SNOVA_o + k1) * SNOVA_l2],
                                  SNOVA_l, SNOVA_l, SNOVA_l);
#endif
}

#endif

#endif
