#ifndef RCT_VERIFY_AVX2_H
#define RCT_VERIFY_AVX2_H

#if RCT_USE_SIMD

static void rct_vf_whip_q16(rct_vf_ctx *c) {
    const gf_t *signature_in_GF = c->sig_gf;
    uint8_t *whipped_sig2 = c->whipped_sig2;
#if RCT_USE_GFNI && RCT_VF_MTK2
    for (int idx = 0; idx < SNOVA_n; ++idx) {
        const gf_t *sg = &signature_in_GF[idx * SNOVA_lr];
        __m256i b0 = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)sg));
        __m256i b1 = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(sg + 2 * SNOVA_r)));
        __m256i R[SNOVA_l];
        for (int k1 = 0; k1 < SNOVA_l; k1++)
            R[k1] = _mm256_shuffle_epi8((k1 < 2) ? b0 : b1,
                                        _mm256_load_si256((const __m256i *)rct_whipR[k1]));
        for (int i1 = 0; i1 < SNOVA_l; i1++) {
            __m256i acc = _mm256_setzero_si256();
            for (int k1 = 0; k1 < SNOVA_l; k1++)
                acc = _mm256_xor_si256(acc, _mm256_gf2p8mul_epi8(
                    _mm256_load_si256((const __m256i *)rct_whipM[i1 * SNOVA_l + k1]), R[k1]));
            _mm256_store_si256((__m256i *)&whipped_sig2[(idx * SNOVA_l + i1) * SNOVA_lr32],
                               rct_gfni_cleanup256(acc));
        }
    }
#elif SNOVA_lr <= 32
    for (int idx = 0; idx < SNOVA_n; ++idx)
        for (int i1 = 0; i1 < SNOVA_l; i1++)
            for (int ab = 0; ab < SNOVA_l; ++ab) {
                __m256i acc = _mm256_setzero_si256();
                for (int k1 = 0; k1 < SNOVA_l; k1++) {
                    __m256i sv = _mm256_loadu_si256((const __m256i *)&signature_in_GF[idx * SNOVA_lr + k1 * SNOVA_r]);
                    __m256i sc = RCT_BC(rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1]);
                    acc = _mm256_xor_si256(acc, RCT_SV(sc, sv));
                }
                _Alignas(32) uint8_t tmp[32];
                _mm256_store_si256((__m256i *)tmp, rct_gfni_cleanup256(acc));
                for (int j1 = 0; j1 < SNOVA_r; j1++)
                    whipped_sig2[idx * SNOVA_l * SNOVA_lr32 + i1 * SNOVA_lr32 + ab * SNOVA_r + j1] = tmp[j1];
            }
#else
    for (int ab = 0; ab < SNOVA_l; ++ab)
        for (int idx = 0; idx < SNOVA_n; ++idx)
            for (int i1 = 0; i1 < SNOVA_l; i1++)
                for (int j1 = 0; j1 < SNOVA_r; j1++)
                    for (int k1 = 0; k1 < SNOVA_l; k1++)
                        gf_set_add(&whipped_sig2[idx * SNOVA_l * SNOVA_lr32 + i1 * SNOVA_lr32 + ab * SNOVA_r + j1],
                                   gf_mult(rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1],
                                           signature_in_GF[idx * SNOVA_lr + k1 * SNOVA_r + j1]));
#endif
}

#if !RCT_VF_JOG
static void rct_vf_reindex_q16(rct_vf_ctx *c) {
    uint8_t *sum_t1p = c->sum_t1p;
#if RCT_VF_EMM
    gf_t *sum_t1q = c->sum_t1q;
#else
    gf_t *sum_t1 = c->sum_t1;
#endif
#if !RCT_VF_MTK2
    for (int i = 0; i < SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32; ++i) sum_t1p[i] = rct_gfni_cleanup(sum_t1p[i]);
#endif
#if RCT_VF_EMM
    for (int mi = 0; mi < SNOVA_m1; ++mi)
        for (int a1 = 0; a1 < SNOVA_l; ++a1)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                rct_vf_tr8(&sum_t1q[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * 64],
                           &sum_t1p[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + b1 * SNOVA_r],
                           SNOVA_lr32);
#else
    for (int mi = 0; mi < SNOVA_m1; ++mi)
        for (int a1 = 0; a1 < SNOVA_l; ++a1)
            for (int i1 = 0; i1 < SNOVA_r; i1++)
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    for (int j1 = 0; j1 < SNOVA_r; j1++)
                        sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1] =
                            sum_t1p[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32 + b1 * SNOVA_r + j1];
#endif
}

#endif

#if !RCT_VF_AQ && !RCT_VF_JOG
static void rct_vf_contract_nonaq(rct_vf_ctx *c) {
    const rct_pk_t *pkx = c->pkx;
    uint8_t *whipped_sig2 = c->whipped_sig2;
    uint8_t *sum_t1p = c->sum_t1p;
#if RCT_VF_MTK2
        static uint16_t rct_vf_ppair[SNOVA_m1 * SNOVA_n * SNOVA_n * 8];
        RCT_SCRATCH _Alignas(32) uint16_t rct_vf_wpair[SNOVA_n * SNOVA_l * 16];
        {
            const __m256i ppl = _mm256_load_si256((const __m256i *)RCT_VF_PPL);
            const __m256i pph = _mm256_load_si256((const __m256i *)RCT_VF_PPH);
            const int ncell = SNOVA_m1 * SNOVA_n * SNOVA_n;
            int c = 0;
            for (; c + 2 <= ncell; c += 2)
                _mm256_storeu_si256(
                    (__m256i *)&rct_vf_ppair[c * 8],
                    _mm256_slli_epi16(
                        _mm256_cvtepu8_epi16(rct_vf_pack_pair32(
                            _mm256_loadu_si256((const __m256i *)&pkx->P[c * SNOVA_l2]), ppl, pph)),
                        4));
            if (c < ncell)
                _mm_storeu_si128(
                    (__m128i *)&rct_vf_ppair[c * 8],
                    _mm256_castsi256_si128(_mm256_slli_epi16(
                        _mm256_cvtepu8_epi16(rct_vf_pack_pair32(
                            _mm256_zextsi128_si256(
                                _mm_loadu_si128((const __m128i *)&pkx->P[c * SNOVA_l2])),
                            ppl, pph)),
                        4)));
            const __m256i wpl = _mm256_load_si256((const __m256i *)RCT_VF_WPL);
            const __m256i wph = _mm256_load_si256((const __m256i *)RCT_VF_WPH);
            for (int rw = 0; rw < SNOVA_n * SNOVA_l; ++rw)
                _mm256_store_si256(
                    (__m256i *)&rct_vf_wpair[rw * 16],
                    _mm256_slli_epi16(
                        _mm256_cvtepu8_epi16(rct_vf_pack_pair32(
                            _mm256_load_si256((const __m256i *)&whipped_sig2[rw * SNOVA_lr32]),
                            wpl, wph)),
                        4));
        }
        for (int mi = 0; mi < SNOVA_m1; ++mi) {
            RCT_SCRATCH _Alignas(32) uint8_t sum_t0[SNOVA_n * SNOVA_l * SNOVA_lr32];
            for (int ni = 0; ni < SNOVA_n; ++ni) {
                __m256i acc01 = _mm256_setzero_si256(), acc23 = _mm256_setzero_si256();
                const uint16_t *pp = &rct_vf_ppair[(mi * SNOVA_n + ni) * SNOVA_n * 8];
                for (int nj = 0; nj < SNOVA_n; ++nj)
                    for (int k1 = 0; k1 < SNOVA_l; k1++) {
                        __m256i wp = _mm256_load_si256(
                            (const __m256i *)&whipped_sig2[(nj * SNOVA_l + k1) * SNOVA_lr32]);
                        uint32_t ii;
                        memcpy(&ii, &pp[nj * 8 + k1 * 2], 4);
                        acc01 = _mm256_xor_si256(
                            acc01, _mm256_shuffle_epi8(rct_mtk2t16((uint16_t)ii), wp));
                        acc23 = _mm256_xor_si256(
                            acc23, _mm256_shuffle_epi8(rct_mtk2t16((uint16_t)(ii >> 16)), wp));
                    }
                uint8_t *s0 = &sum_t0[ni * SNOVA_l * SNOVA_lr32];
                _mm256_store_si256((__m256i *)(s0 + 0 * SNOVA_lr32), rct_nib_lo(acc01));
                _mm256_store_si256((__m256i *)(s0 + 1 * SNOVA_lr32), rct_nib_hi(acc01));
                _mm256_store_si256((__m256i *)(s0 + 2 * SNOVA_lr32), rct_nib_lo(acc23));
                _mm256_store_si256((__m256i *)(s0 + 3 * SNOVA_lr32), rct_nib_hi(acc23));
            }
            for (int h = 0; h < 2; ++h) {
                __m256i acc[SNOVA_r];
                for (int p = 0; p < SNOVA_r; p++) acc[p] = _mm256_setzero_si256();
                for (int ni = 0; ni < SNOVA_n; ++ni)
                    for (int k1 = 0; k1 < SNOVA_l; k1++) {
                        __m256i s0 = _mm256_load_si256(
                            (const __m256i *)&sum_t0[(ni * SNOVA_l + k1) * SNOVA_lr32]);
                        const uint16_t *wq = &rct_vf_wpair[(ni * SNOVA_l + k1) * 16 + h * SNOVA_r];
                        uint64_t w0, w1;
                        memcpy(&w0, wq, 8);
                        memcpy(&w1, wq + 4, 8);
                        for (int p = 0; p < SNOVA_r; p++) {
                            uint16_t idx = (uint16_t)((p < 4 ? (w0 >> (16 * p))
                                                             : (w1 >> (16 * (p - 4)))));
                            acc[p] = _mm256_xor_si256(acc[p],
                                                      _mm256_shuffle_epi8(rct_mtk2t16(idx), s0));
                        }
                    }
                for (int p = 0; p < SNOVA_r; p++) {
                    int j0 = h * 2 * SNOVA_r + 2 * p;
                    _mm256_store_si256(
                        (__m256i *)&sum_t1p[(mi * SNOVA_l * SNOVA_r + j0) * SNOVA_lr32],
                        rct_nib_lo(acc[p]));
                    _mm256_store_si256(
                        (__m256i *)&sum_t1p[(mi * SNOVA_l * SNOVA_r + j0 + 1) * SNOVA_lr32],
                        rct_nib_hi(acc[p]));
                }
            }
        }
#elif SNOVA_lr16 == 1
        for (int mi = 0; mi < SNOVA_m1; ++mi) {
            RCT_SCRATCH _Alignas(32) uint8_t sum_t0[SNOVA_n * SNOVA_l * SNOVA_lr32];
            for (int ni = 0; ni < SNOVA_n; ++ni) {
                __m256i acc[SNOVA_l];
                for (int i1 = 0; i1 < SNOVA_l; i1++) acc[i1] = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_n; ++nj)
                    for (int k1 = 0; k1 < SNOVA_l; k1++) {
                        __m256i wp = _mm256_load_si256(
                            (const __m256i *)&whipped_sig2[(nj * SNOVA_l + k1) * SNOVA_lr32]);
                        const gf_t *prow =
                            &pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + k1];
                        for (int i1 = 0; i1 < SNOVA_l; i1++)
                            acc[i1] = _mm256_xor_si256(acc[i1], RCT_SV(RCT_BC(prow[i1 * SNOVA_l]), wp));
                    }
                for (int i1 = 0; i1 < SNOVA_l; i1++)
                    _mm256_store_si256((__m256i *)&sum_t0[(ni * SNOVA_l + i1) * SNOVA_lr32],
                                       rct_gfni_cleanup256(acc[i1]));
            }
            for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                __m256i acc[SNOVA_r];
                for (int i1 = 0; i1 < SNOVA_r; i1++) acc[i1] = _mm256_setzero_si256();
                for (int ni = 0; ni < SNOVA_n; ++ni)
                    for (int k1 = 0; k1 < SNOVA_l; k1++) {
                        __m256i s0 = _mm256_load_si256(
                            (const __m256i *)&sum_t0[(ni * SNOVA_l + k1) * SNOVA_lr32]);
                        const gf_t *wrow =
                            &whipped_sig2[(ni * SNOVA_l + k1) * SNOVA_lr32 + a1 * SNOVA_r];
                        for (int i1 = 0; i1 < SNOVA_r; i1++)
                            acc[i1] = _mm256_xor_si256(acc[i1], RCT_SV(RCT_BC(wrow[i1]), s0));
                    }
                for (int i1 = 0; i1 < SNOVA_r; i1++)
                    _mm256_store_si256(
                        (__m256i *)&sum_t1p[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32],
                        acc[i1]);
            }
        }
#else
        for (int mi = 0; mi < SNOVA_m1; ++mi) {
            RCT_SCRATCH _Alignas(32) uint8_t sum_t0[SNOVA_n * SNOVA_l * SNOVA_lr32];
            memset(sum_t0, 0, sizeof(sum_t0));
            for (int ni = 0; ni < SNOVA_n; ++ni)
                for (int nj = 0; nj < SNOVA_n; ++nj)
                    for (int i1 = 0; i1 < SNOVA_l; i1++)
                        for (int k1 = 0; k1 < SNOVA_l; k1++)
                            for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
                                __m256i *s0 = (__m256i *)&sum_t0[(ni * SNOVA_l + i1) * SNOVA_lr32];
                                __m256i q1v = RCT_BC(
                                    pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
                                __m256i *wp = (__m256i *)&whipped_sig2[nj * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32];
                                s0[b1] = _mm256_xor_si256(s0[b1], RCT_SV(q1v, wp[b1]));
                            }
            for (int i = 0; i < SNOVA_n * SNOVA_l * SNOVA_lr32; ++i) sum_t0[i] = rct_gfni_cleanup(sum_t0[i]);
            for (int ni = 0; ni < SNOVA_n; ++ni)
                for (int a1 = 0; a1 < SNOVA_l; ++a1)
                    for (int k1 = 0; k1 < SNOVA_l; k1++)
                        for (int i1 = 0; i1 < SNOVA_r; i1++)
                            for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
                                __m256i *s1 = (__m256i *)&sum_t1p[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32];
                                __m256i wp = RCT_BC(
                                    whipped_sig2[ni * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32 + a1 * SNOVA_r + i1]);
                                __m256i *s0 = (__m256i *)&sum_t0[(ni * SNOVA_l + k1) * SNOVA_lr32];
                                s1[b1] = _mm256_xor_si256(s1[b1], RCT_SV(wp, s0[b1]));
                            }
        }
#endif
}
#endif

#endif

#endif
