#ifndef RCT_VERIFY_JOG_H
#define RCT_VERIFY_JOG_H

#if RCT_VF_JOG

static inline __m256i rct_jog_sv256(gf_t s, __m256i v) {
    return _mm256_shuffle_epi8(rct_mtk2t((uint8_t)s), v);
}
static inline __m128i rct_jog_sv128(gf_t s, __m128i v) {
    return _mm_shuffle_epi8(_mm_load_si128((const __m128i *)rct_mtk2[(uint8_t)s]), v);
}

static void rct_vf_jog(rct_vf_ctx *c) {
#if RCT_USE_SIMD
    const uint8_t *W2 = c->whipped_sig2;
#define RCT_JOG_W(ab, ni, row, col) \
    (W2[(((size_t)(ni) * SNOVA_l + (row)) * SNOVA_lr32) + (size_t)(ab) * SNOVA_r + (col)])
#else
    RCT_SCRATCH _Alignas(32) gf_t jog_whip[SNOVA_l * SNOVA_n * SNOVA_lr];
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
                memcpy(&jog_whip[((size_t)ab * SNOVA_n + idx) * SNOVA_lr + i1 * SNOVA_r],
                       wb, SNOVA_r);
            }
        }
#define RCT_JOG_W(ab, ni, row, col) \
    (jog_whip[(((size_t)(ab) * SNOVA_n + (ni)) * SNOVA_lr) + (size_t)(row) * SNOVA_r + (col)])
#endif

    RCT_SCRATCH _Alignas(32) uint16_t jog_wpair[SNOVA_l][RCT_JOG_RP][RCT_JOG_NL];
    RCT_SCRATCH _Alignas(32) uint8_t jog_WALL[RCT_JOG_NL][SNOVA_lr32];
    memset(jog_WALL, 0, sizeof(jog_WALL));
    for (int a1 = 0; a1 < SNOVA_l; ++a1)
        for (int ni = 0; ni < SNOVA_n; ++ni)
            for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                const int p = ni * SNOVA_l + k1;
                for (int i1p = 0; i1p < RCT_JOG_RP; ++i1p) {
                    const unsigned lo = RCT_JOG_W(a1, ni, k1, 2 * i1p);
                    const unsigned hi =
                        (2 * i1p + 1 < SNOVA_r) ? RCT_JOG_W(a1, ni, k1, 2 * i1p + 1) : 0u;
                    jog_wpair[a1][i1p][p] = (uint16_t)(((hi << 4) | lo) << 4);
                }
            }
    for (int nj = 0; nj < SNOVA_n; ++nj)
        for (int ej = 0; ej < SNOVA_l; ++ej) {
            const int cc = nj * SNOVA_l + ej;
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int j1 = 0; j1 < SNOVA_r; ++j1)
                    jog_WALL[cc][b1 * SNOVA_r + j1] = (uint8_t)RCT_JOG_W(b1, nj, ej, j1);
        }

#if !RCT_JOG_PKXJOG
    RCT_SCRATCH _Alignas(32) uint8_t jog_PJ[RCT_JOG_NL][RCT_JOG_L32];
    memset(jog_PJ, 0, sizeof(jog_PJ));
#endif
    RCT_SCRATCH _Alignas(32) uint8_t jog_tall[SNOVA_l][RCT_JOG_RP][RCT_JOG_L32];
    memset(jog_tall, 0, sizeof(jog_tall));

    gf_t *sum_t1 = c->sum_t1;
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
#if RCT_JOG_PKXJOG
        const uint8_t *PJ = &c->pkx->P[(size_t)mi * RCT_JOG_NL * RCT_JOG_L32];
#define RCT_JOG_PROW(p) (PJ + (size_t)(p) * RCT_JOG_L32)
#else
        {
            const __m256i perm = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);
            for (int ni = 0; ni < SNOVA_n; ++ni) {
                const gf_t *crow =
                    &c->pkx->P[(((size_t)mi * SNOVA_n + ni) * SNOVA_n) * SNOVA_l2];
                uint8_t *prow = jog_PJ[ni * SNOVA_l];
                int nj = 0;
                for (; nj + 8 <= SNOVA_n; nj += 8) {
                    const uint8_t *src = crow + (size_t)nj * SNOVA_l2;
                    __m256i A = _mm256_loadu_si256((const __m256i *)(src + 0));
                    __m256i B = _mm256_loadu_si256((const __m256i *)(src + 32));
                    __m256i C = _mm256_loadu_si256((const __m256i *)(src + 64));
                    __m256i D = _mm256_loadu_si256((const __m256i *)(src + 96));
                    __m256i t0 = _mm256_unpacklo_epi32(A, B);
                    __m256i t1 = _mm256_unpackhi_epi32(A, B);
                    __m256i t2 = _mm256_unpacklo_epi32(C, D);
                    __m256i t3 = _mm256_unpackhi_epi32(C, D);
                    __m256i u0 = _mm256_unpacklo_epi64(t0, t2);
                    __m256i u1 = _mm256_unpackhi_epi64(t0, t2);
                    __m256i u2 = _mm256_unpacklo_epi64(t1, t3);
                    __m256i u3 = _mm256_unpackhi_epi64(t1, t3);
                    _mm256_storeu_si256((__m256i *)(prow + 0 * RCT_JOG_L32 + nj * 4),
                                        _mm256_permutevar8x32_epi32(u0, perm));
                    _mm256_storeu_si256((__m256i *)(prow + 1 * RCT_JOG_L32 + nj * 4),
                                        _mm256_permutevar8x32_epi32(u1, perm));
                    _mm256_storeu_si256((__m256i *)(prow + 2 * RCT_JOG_L32 + nj * 4),
                                        _mm256_permutevar8x32_epi32(u2, perm));
                    _mm256_storeu_si256((__m256i *)(prow + 3 * RCT_JOG_L32 + nj * 4),
                                        _mm256_permutevar8x32_epi32(u3, perm));
                }
                for (; nj < SNOVA_n; ++nj) {
                    const gf_t *cell = crow + (size_t)nj * SNOVA_l2;
                    for (int ei = 0; ei < SNOVA_l; ++ei)
                        memcpy(prow + (size_t)ei * RCT_JOG_L32 + (size_t)nj * SNOVA_l,
                               cell + (size_t)ei * SNOVA_l, SNOVA_l);
                }
            }
        }
#define RCT_JOG_PROW(p) (jog_PJ[(p)])
#endif
#define RCT_JOG_BG ((RCT_JOG_VTL <= 5) ? RCT_JOG_VTL : 4)
        for (int vb = 0; vb < RCT_JOG_VTL; vb += RCT_JOG_BG) {
            const int gw = (vb + RCT_JOG_BG <= RCT_JOG_VTL) ? RCT_JOG_BG : (RCT_JOG_VTL - vb);
            for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                int i1p = 0;
                for (; i1p + 2 <= RCT_JOG_RP; i1p += 2) {
                    __m256i acc0[RCT_JOG_BG], acc1[RCT_JOG_BG];
                    for (int g = 0; g < gw; ++g) {
                        acc0[g] = _mm256_setzero_si256();
                        acc1[g] = _mm256_setzero_si256();
                    }
                    const uint16_t *wp0 = jog_wpair[a1][i1p];
                    const uint16_t *wp1 = jog_wpair[a1][i1p + 1];
                    for (int p = 0; p < RCT_JOG_NL; ++p) {
                        const __m256i tb0 = rct_mtk2t16(wp0[p]);
                        const __m256i tb1 = rct_mtk2t16(wp1[p]);
                        const uint8_t *pr = RCT_JOG_PROW(p) + (size_t)vb * 32;
                        for (int g = 0; g < gw; ++g) {
                            const __m256i pv =
                                _mm256_loadu_si256((const __m256i *)(pr + g * 32));
                            acc0[g] = _mm256_xor_si256(acc0[g], _mm256_shuffle_epi8(tb0, pv));
                            acc1[g] = _mm256_xor_si256(acc1[g], _mm256_shuffle_epi8(tb1, pv));
                        }
                    }
                    for (int g = 0; g < gw; ++g) {
                        _mm256_store_si256(
                            (__m256i *)&jog_tall[a1][i1p][(vb + g) * 32], acc0[g]);
                        _mm256_store_si256(
                            (__m256i *)&jog_tall[a1][i1p + 1][(vb + g) * 32], acc1[g]);
                    }
                }
                for (; i1p < RCT_JOG_RP; ++i1p) {
                    __m256i acc[RCT_JOG_BG];
                    for (int g = 0; g < gw; ++g) acc[g] = _mm256_setzero_si256();
                    const uint16_t *wp = jog_wpair[a1][i1p];
                    for (int p = 0; p < RCT_JOG_NL; ++p) {
                        const __m256i tbl = rct_mtk2t16(wp[p]);
                        const uint8_t *pr = RCT_JOG_PROW(p) + (size_t)vb * 32;
                        for (int g = 0; g < gw; ++g)
                            acc[g] = _mm256_xor_si256(
                                acc[g], _mm256_shuffle_epi8(
                                            tbl, _mm256_loadu_si256((const __m256i *)(pr + g * 32))));
                    }
                    for (int g = 0; g < gw; ++g)
                        _mm256_store_si256((__m256i *)&jog_tall[a1][i1p][(vb + g) * 32], acc[g]);
                }
            }
        }
#undef RCT_JOG_BG
        for (int a1 = 0; a1 < SNOVA_l; ++a1) {
            __m256i accB[RCT_JOG_RP][SNOVA_lr16];
            for (int i1p = 0; i1p < RCT_JOG_RP; ++i1p)
                for (int w = 0; w < SNOVA_lr16; ++w) accB[i1p][w] = _mm256_setzero_si256();
            for (int cc = 0; cc < RCT_JOG_NL; ++cc) {
                __m256i wall[SNOVA_lr16];
                for (int w = 0; w < SNOVA_lr16; ++w)
                    wall[w] = _mm256_load_si256((const __m256i *)&jog_WALL[cc][w * 32]);
                for (int i1p = 0; i1p < RCT_JOG_RP; ++i1p) {
                    const __m256i tbl = rct_mtk2t(jog_tall[a1][i1p][cc]);
                    for (int w = 0; w < SNOVA_lr16; ++w)
                        accB[i1p][w] =
                            _mm256_xor_si256(accB[i1p][w], _mm256_shuffle_epi8(tbl, wall[w]));
                }
            }
            gf_t *st1_a = &sum_t1[((size_t)mi * SNOVA_l2 + (size_t)a1 * SNOVA_l) * SNOVA_r2];
            for (int i1p = 0; i1p < RCT_JOG_RP; ++i1p) {
                _Alignas(32) uint8_t rowlo[SNOVA_lr32], rowhi[SNOVA_lr32];
                for (int w = 0; w < SNOVA_lr16; ++w) {
                    _mm256_store_si256((__m256i *)&rowlo[w * 32], rct_nib_lo(accB[i1p][w]));
                    _mm256_store_si256((__m256i *)&rowhi[w * 32], rct_nib_hi(accB[i1p][w]));
                }
                for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                    gf_t *st1 = st1_a + (size_t)b1 * SNOVA_r2;
                    for (int j1 = 0; j1 < SNOVA_r; ++j1) {
                        st1[(2 * i1p) * SNOVA_r + j1] = rowlo[b1 * SNOVA_r + j1];
                        if (2 * i1p + 1 < SNOVA_r)
                            st1[(2 * i1p + 1) * SNOVA_r + j1] = rowhi[b1 * SNOVA_r + j1];
                    }
                }
            }
        }
#undef RCT_JOG_PROW
    }
#undef RCT_JOG_W
}

_Static_assert(SNOVA_r2 <= 64, "rct_vf_emat_jog: temp1[64]/2-ymm q1q2 fold assumes r2<=64");
static void rct_vf_emat_jog(rct_vf_ctx *c) {
    const rct_pk_t *pkx = c->pkx;
    const gf_t *sum_t1 = c->sum_t1;
    gf_t *hash_in_GF = c->hash_gf;
    for (int mi = 0; mi < SNOVA_o; ++mi) {
        for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
            const int mi_prime = i_prime(mi, alpha);
            const gf_t *q1 = &pkx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l];
            const gf_t *q2 = &pkx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l];
            _Alignas(32) uint8_t temp1[64];
            __m256i t1lo = _mm256_setzero_si256(), t1hi = _mm256_setzero_si256();
            for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                __m256i slo = _mm256_setzero_si256(), shi = _mm256_setzero_si256();
                for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                    const gf_t *base =
                        &sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2];
                    slo = _mm256_xor_si256(slo,
                        rct_jog_sv256(q2[b1], _mm256_loadu_si256((const __m256i *)base)));
                    shi = _mm256_xor_si256(shi,
                        rct_jog_sv256(q2[b1], _mm256_loadu_si256((const __m256i *)(base + 32))));
                }
                t1lo = _mm256_xor_si256(t1lo, rct_jog_sv256(q1[a1], slo));
                t1hi = _mm256_xor_si256(t1hi, rct_jog_sv256(q1[a1], shi));
            }
            _mm256_store_si256((__m256i *)temp1, t1lo);
            _mm256_store_si256((__m256i *)(temp1 + 32), t1hi);
            const gf_t *Bm = &pkx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr];
            const gf_t *Am = &pkx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2];
            _Alignas(16) uint8_t temp2[SNOVA_r * 16];
            for (int i1 = 0; i1 < SNOVA_r; ++i1) {
                __m128i acc = _mm_setzero_si128();
                for (int k1 = 0; k1 < SNOVA_r; ++k1)
                    acc = _mm_xor_si128(acc,
                        rct_jog_sv128(temp1[i1 * SNOVA_r + k1],
                                      _mm_loadu_si128((const __m128i *)&Bm[k1 * SNOVA_l])));
                _mm_store_si128((__m128i *)&temp2[i1 * 16], acc);
            }
            gf_t *hrow = &hash_in_GF[mi * SNOVA_lr];
            for (int i1 = 0; i1 < SNOVA_r; ++i1) {
                _Alignas(16) uint8_t hb[16];
                __m128i acc = _mm_setzero_si128();
                for (int k1 = 0; k1 < SNOVA_r; ++k1)
                    acc = _mm_xor_si128(acc,
                        rct_jog_sv128(Am[i1 * SNOVA_r + k1],
                                      _mm_load_si128((const __m128i *)&temp2[k1 * 16])));
                _mm_store_si128((__m128i *)hb, acc);
                for (int j = 0; j < SNOVA_l; ++j)
                    hrow[i1 * SNOVA_l + j] = (gf_t)(hrow[i1 * SNOVA_l + j] ^ hb[j]);
            }
        }
    }
}

#endif

#endif
