#ifndef RCT_VERIFY_ODDQ_H
#define RCT_VERIFY_ODDQ_H

#if RCT_Q_SIMD

#if !RCT_Q_MADD
static _Alignas(32) uint16_t whipw[SNOVA_l * SNOVA_n * RCT_Q_LRP];
#endif
static _Alignas(32) uint16_t sum_t1w[SNOVA_m1 * SNOVA_l * SNOVA_r * RCT_Q_LRP];

static void rct_vf_whip_oddq(rct_vf_ctx *c) {
    const gf_t *signature_in_GF = c->sig_gf;
#if RCT_Q_MADD
        rct_qv_build();
        for (int idx = 0; idx < SNOVA_n; ++idx) {
            __m256i se[SNOVA_l][RCT_Q_LR16];
            for (int k1 = 0; k1 < SNOVA_l; k1++) {
                long long sq;
                memcpy(&sq, &signature_in_GF[idx * SNOVA_lr + k1 * SNOVA_r], 8);
                __m256i bq = _mm256_set1_epi64x(sq);
                for (int g = 0; g < RCT_Q_LR16; g++)
                    se[k1][g] = _mm256_shuffle_epi8(bq, _mm256_load_si256((const __m256i *)rct_qv_sigpat[g]));
            }
#if RCT_Q_LR16 == 1
            __m128i row[SNOVA_l];
            for (int i1 = 0; i1 < SNOVA_l; i1++) {
                __m256i acc = _mm256_mullo_epi16(
                    _mm256_load_si256((const __m256i *)rct_qv_sseg[i1 * SNOVA_l]), se[0][0]);
                for (int k1 = 1; k1 < SNOVA_l; k1++)
                    acc = _mm256_add_epi16(acc, _mm256_mullo_epi16(
                        _mm256_load_si256((const __m256i *)rct_qv_sseg[i1 * SNOVA_l + k1]), se[k1][0]));
                row[i1] = rct_qv_pack16(rct_q_barrett16(acc));
            }
            uint8_t *wd = &rct_qv_wpair[idx * 2 * 2 * RCT_Q_LRP];
            _mm_store_si128((__m128i *)wd, _mm_unpacklo_epi8(row[0], row[1]));
            _mm_store_si128((__m128i *)(wd + 16), _mm_unpackhi_epi8(row[0], row[1]));
            _mm_store_si128((__m128i *)(wd + 32), _mm_unpacklo_epi8(row[2], row[3]));
            _mm_store_si128((__m128i *)(wd + 48), _mm_unpackhi_epi8(row[2], row[3]));
#if RCT_Q_VNNI
            {
                __m256i p0 = _mm256_load_si256((const __m256i *)wd);
                __m256i p1 = _mm256_load_si256((const __m256i *)(wd + 32));
                __m256i lo = _mm256_unpacklo_epi16(p0, p1), hi = _mm256_unpackhi_epi16(p0, p1);
                uint8_t *qd = &rct_qv_wquad[idx * 4 * RCT_Q_LRP];
                _mm256_store_si256((__m256i *)qd, _mm256_permute2x128_si256(lo, hi, 0x20));
                _mm256_store_si256((__m256i *)(qd + 32), _mm256_permute2x128_si256(lo, hi, 0x31));
            }
#endif
#else
            __m256i row[SNOVA_l];
            for (int i1 = 0; i1 < SNOVA_l; i1++) {
                __m256i a0 = _mm256_mullo_epi16(
                    _mm256_load_si256((const __m256i *)&rct_qv_sseg[i1 * SNOVA_l][0]), se[0][0]);
                __m256i a1 = _mm256_mullo_epi16(
                    _mm256_load_si256((const __m256i *)&rct_qv_sseg[i1 * SNOVA_l][16]), se[0][1]);
                for (int k1 = 1; k1 < SNOVA_l; k1++) {
                    a0 = _mm256_add_epi16(a0, _mm256_mullo_epi16(
                        _mm256_load_si256((const __m256i *)&rct_qv_sseg[i1 * SNOVA_l + k1][0]), se[k1][0]));
                    a1 = _mm256_add_epi16(a1, _mm256_mullo_epi16(
                        _mm256_load_si256((const __m256i *)&rct_qv_sseg[i1 * SNOVA_l + k1][16]), se[k1][1]));
                }
                row[i1] = rct_qv_pack32(rct_q_barrett16(a0), rct_q_barrett16(a1));
            }
            uint8_t *wd = &rct_qv_wpair[idx * 2 * 2 * RCT_Q_LRP];
            rct_qv_ilv32(wd, row[0], row[1]);
            rct_qv_ilv32(wd + 2 * RCT_Q_LRP, row[2], row[3]);
#if RCT_Q_VNNI
            {
                uint8_t *qd = &rct_qv_wquad[idx * 4 * RCT_Q_LRP];
                for (int blk = 0; blk < 2; ++blk) {
                    __m256i p0 = _mm256_load_si256((const __m256i *)(wd + blk * 32));
                    __m256i p1 = _mm256_load_si256((const __m256i *)(wd + 2 * RCT_Q_LRP + blk * 32));
                    __m256i lo = _mm256_unpacklo_epi16(p0, p1), hi = _mm256_unpackhi_epi16(p0, p1);
                    _mm256_store_si256((__m256i *)(qd + blk * 64),
                                       _mm256_permute2x128_si256(lo, hi, 0x20));
                    _mm256_store_si256((__m256i *)(qd + blk * 64 + 32),
                                       _mm256_permute2x128_si256(lo, hi, 0x31));
                }
            }
#endif
#endif
        }
#else
        memset(whipw, 0, sizeof(whipw));
        for (int ab = 0; ab < SNOVA_l; ++ab)
            for (int idx = 0; idx < SNOVA_n; ++idx)
                for (int i1 = 0; i1 < SNOVA_l; i1++)
                    for (int j1 = 0; j1 < SNOVA_r; j1++)
                        for (int k1 = 0; k1 < SNOVA_l; k1++)
                            whipw[idx * SNOVA_l * RCT_Q_LRP + i1 * RCT_Q_LRP + ab * SNOVA_r + j1] +=
                                (uint16_t)rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1] *
                                (uint16_t)signature_in_GF[idx * SNOVA_lr + k1 * SNOVA_r + j1];
        for (int i1 = 0; i1 < SNOVA_l * SNOVA_n * RCT_Q_LRP; i1++) whipw[i1] %= SNOVA_q;
#endif
}

static void rct_vf_rl_oddq(rct_vf_ctx *c) {
    const rct_pk_t *pkx = c->pkx;

#if !RCT_Q_MADD
        memset(sum_t1w, 0, sizeof(sum_t1w));
#endif
#if RCT_VERIFY_STREAM
        {
            _Alignas(32) static uint16_t rct_vf_qst0[SNOVA_m1 * SNOVA_n * SNOVA_l * RCT_Q_LRP];
            memset(rct_vf_qst0, 0, sizeof(rct_vf_qst0));
            snova_pgen_t pg;
            snova_pgen_init(&pg, pkx->pk_seed);
            snova_prow_t rw;
            while (snova_pgen_next_row(&pg, &rw)) {
                const int mi = rw.mi;
                const int ni = (rw.block == 2) ? SNOVA_v + rw.ni : rw.ni;
                const int col0 = (rw.block == 1) ? SNOVA_v : 0;
                uint16_t *am = &rct_vf_qst0[((size_t)mi * SNOVA_n + ni) * SNOVA_l * RCT_Q_LRP];
                rct_vf_qs_rowseg(am, rw.cells, col0, rw.ncols);
                if (rw.block == 2) {
                    rct_vf_qs_rowseg(am, &pkx->P22[((size_t)mi * SNOVA_o + rw.ni) * SNOVA_o * SNOVA_l2],
                                     SNOVA_v, SNOVA_o);
                    if (rw.ni == SNOVA_o - 1) {
                        for (int nn = 0; nn < SNOVA_n; ++nn) {
                            const uint16_t *am2 =
                                &rct_vf_qst0[((size_t)mi * SNOVA_n + nn) * SNOVA_l * RCT_Q_LRP];
                            uint8_t *const s0p[2] = { &rct_qv_s0all[nn * 4 * RCT_Q_LRP],
                                                      &rct_qv_s0all[nn * 4 * RCT_Q_LRP + 2 * RCT_Q_LRP] };
#if RCT_Q_LR16 == 1
                            __m128i r0 = rct_qv_pack16(rct_q_barrett16(_mm256_load_si256((const __m256i *)&am2[0])));
                            __m128i r1 = rct_qv_pack16(rct_q_barrett16(_mm256_load_si256((const __m256i *)&am2[16])));
                            __m128i r2 = rct_qv_pack16(rct_q_barrett16(_mm256_load_si256((const __m256i *)&am2[32])));
                            __m128i r3 = rct_qv_pack16(rct_q_barrett16(_mm256_load_si256((const __m256i *)&am2[48])));
                            _mm_store_si128((__m128i *)&s0p[0][0], _mm_unpacklo_epi8(r0, r1));
                            _mm_store_si128((__m128i *)&s0p[0][16], _mm_unpackhi_epi8(r0, r1));
                            _mm_store_si128((__m128i *)&s0p[1][0], _mm_unpacklo_epi8(r2, r3));
                            _mm_store_si128((__m128i *)&s0p[1][16], _mm_unpackhi_epi8(r2, r3));
#else
                            for (int hh = 0; hh < 2; ++hh) {
                                __m256i A = rct_qv_pack32(
                                    rct_q_barrett16(_mm256_load_si256(
                                        (const __m256i *)&am2[(2 * hh) * 2 * 16])),
                                    rct_q_barrett16(_mm256_load_si256(
                                        (const __m256i *)&am2[(2 * hh) * 2 * 16 + 16])));
                                __m256i B = rct_qv_pack32(
                                    rct_q_barrett16(_mm256_load_si256(
                                        (const __m256i *)&am2[(2 * hh + 1) * 2 * 16])),
                                    rct_q_barrett16(_mm256_load_si256(
                                        (const __m256i *)&am2[(2 * hh + 1) * 2 * 16 + 16])));
                                rct_qv_ilv32(s0p[hh], A, B);
                            }
#endif
                        }
                        for (int col = 0; col < SNOVA_l * SNOVA_r; ++col) {
                            __m256i v0 = _mm256_setzero_si256();
#if RCT_Q_LR16 == 2
                            __m256i v1 = _mm256_setzero_si256();
#endif
                            for (int nn = 0; nn < SNOVA_n; ++nn) {
                                const uint8_t *sb = &rct_qv_s0all[nn * 4 * RCT_Q_LRP];
                                const uint16_t *wc = (const uint16_t *)&rct_qv_wpair[nn * 4 * RCT_Q_LRP];
                                __m256i b0 = _mm256_set1_epi16((short)wc[col]);
                                __m256i b1 = _mm256_set1_epi16((short)wc[RCT_Q_LRP + col]);
                                v0 = _mm256_add_epi16(v0, _mm256_maddubs_epi16(
                                    _mm256_load_si256((const __m256i *)sb), b0));
                                v0 = _mm256_add_epi16(v0, _mm256_maddubs_epi16(
                                    _mm256_load_si256((const __m256i *)(sb + 2 * RCT_Q_LRP)), b1));
#if RCT_Q_LR16 == 2
                                v1 = _mm256_add_epi16(v1, _mm256_maddubs_epi16(
                                    _mm256_load_si256((const __m256i *)(sb + 32)), b0));
                                v1 = _mm256_add_epi16(v1, _mm256_maddubs_epi16(
                                    _mm256_load_si256((const __m256i *)(sb + 2 * RCT_Q_LRP + 32)), b1));
#endif
                            }
                            __m256i *s1 = (__m256i *)&sum_t1w[(mi * SNOVA_l * SNOVA_r + col) * RCT_Q_LRP];
                            _mm256_store_si256(s1, v0);
#if RCT_Q_LR16 == 2
                            _mm256_store_si256(s1 + 1, v1);
#endif
                        }
                    }
                }
            }
        }
#elif RCT_Q_MADD
        for (int mi = 0; mi < SNOVA_m1; ++mi) {
            for (int ni = 0; ni < SNOVA_n; ++ni) {
                const gf_t *pc = &pkx->P[(mi * SNOVA_n + ni) * SNOVA_n * SNOVA_l2];
                uint8_t *const s0p[2] = { &rct_qv_s0all[ni * 4 * RCT_Q_LRP],
                                          &rct_qv_s0all[ni * 4 * RCT_Q_LRP + 2 * RCT_Q_LRP] };
#if RCT_Q_LR16 == 1
#if RCT_Q_VNNI
                __m256i qa[SNOVA_l][2];
                for (int i1 = 0; i1 < SNOVA_l; i1++) {
                    qa[i1][0] = _mm256_setzero_si256();
                    qa[i1][1] = _mm256_setzero_si256();
                }
                for (int nj = 0; nj < SNOVA_n; ++nj) {
                    const uint8_t *wq = &rct_qv_wquad[nj * 4 * RCT_Q_LRP];
                    __m256i w0 = _mm256_load_si256((const __m256i *)wq);
                    __m256i w1 = _mm256_load_si256((const __m256i *)(wq + 32));
                    const gf_t *pcell = pc + nj * SNOVA_l2;
                    for (int i1 = 0; i1 < SNOVA_l; i1++) {
                        int32_t pd;
                        memcpy(&pd, pcell + i1 * SNOVA_l, 4);
                        __m256i pq = _mm256_set1_epi32(pd);
                        qa[i1][0] = _mm256_dpbusd_avx_epi32(qa[i1][0], w0, pq);
                        qa[i1][1] = _mm256_dpbusd_avx_epi32(qa[i1][1], w1, pq);
                    }
                }
                {
                    __m128i r0 = rct_qv_pack16(rct_q_barrett16(rct_qv_low16(qa[0][0], qa[0][1])));
                    __m128i r1 = rct_qv_pack16(rct_q_barrett16(rct_qv_low16(qa[1][0], qa[1][1])));
                    __m128i r2 = rct_qv_pack16(rct_q_barrett16(rct_qv_low16(qa[2][0], qa[2][1])));
                    __m128i r3 = rct_qv_pack16(rct_q_barrett16(rct_qv_low16(qa[3][0], qa[3][1])));
                    _mm_store_si128((__m128i *)&s0p[0][0], _mm_unpacklo_epi8(r0, r1));
                    _mm_store_si128((__m128i *)&s0p[0][16], _mm_unpackhi_epi8(r0, r1));
                    _mm_store_si128((__m128i *)&s0p[1][0], _mm_unpacklo_epi8(r2, r3));
                    _mm_store_si128((__m128i *)&s0p[1][16], _mm_unpackhi_epi8(r2, r3));
                }
#else
                __m256i acc0 = _mm256_setzero_si256(), acc1 = _mm256_setzero_si256();
                __m256i acc2 = _mm256_setzero_si256(), acc3 = _mm256_setzero_si256();
                for (int nj = 0; nj < SNOVA_n; ++nj) {
                    const uint8_t *wb = &rct_qv_wpair[nj * 2 * 2 * RCT_Q_LRP];
                    const gf_t *pcell = pc + nj * SNOVA_l2;
                    for (int h = 0; h < 2; ++h) {
                        __m256i w0 = _mm256_load_si256((const __m256i *)(wb + h * 2 * RCT_Q_LRP));
                        uint16_t pw;
                        memcpy(&pw, pcell + 0 * SNOVA_l + 2 * h, 2);
                        acc0 = _mm256_add_epi16(acc0, _mm256_maddubs_epi16(w0, _mm256_set1_epi16((short)pw)));
                        memcpy(&pw, pcell + 1 * SNOVA_l + 2 * h, 2);
                        acc1 = _mm256_add_epi16(acc1, _mm256_maddubs_epi16(w0, _mm256_set1_epi16((short)pw)));
                        memcpy(&pw, pcell + 2 * SNOVA_l + 2 * h, 2);
                        acc2 = _mm256_add_epi16(acc2, _mm256_maddubs_epi16(w0, _mm256_set1_epi16((short)pw)));
                        memcpy(&pw, pcell + 3 * SNOVA_l + 2 * h, 2);
                        acc3 = _mm256_add_epi16(acc3, _mm256_maddubs_epi16(w0, _mm256_set1_epi16((short)pw)));
                    }
                }
                {
                    __m128i r0 = rct_qv_pack16(rct_q_barrett16(acc0));
                    __m128i r1 = rct_qv_pack16(rct_q_barrett16(acc1));
                    __m128i r2 = rct_qv_pack16(rct_q_barrett16(acc2));
                    __m128i r3 = rct_qv_pack16(rct_q_barrett16(acc3));
                    _mm_store_si128((__m128i *)&s0p[0][0], _mm_unpacklo_epi8(r0, r1));
                    _mm_store_si128((__m128i *)&s0p[0][16], _mm_unpackhi_epi8(r0, r1));
                    _mm_store_si128((__m128i *)&s0p[1][0], _mm_unpacklo_epi8(r2, r3));
                    _mm_store_si128((__m128i *)&s0p[1][16], _mm_unpackhi_epi8(r2, r3));
                }
#endif
#else
#if RCT_Q_VNNI
                __m256i row[SNOVA_l];
                for (int pass = 0; pass < 2; ++pass) {
                    __m256i qa0[4], qa1[4];
                    for (int g = 0; g < 4; g++) {
                        qa0[g] = _mm256_setzero_si256();
                        qa1[g] = _mm256_setzero_si256();
                    }
                    for (int nj = 0; nj < SNOVA_n; ++nj) {
                        const uint8_t *wq = &rct_qv_wquad[nj * 4 * RCT_Q_LRP];
                        const gf_t *pcell = pc + nj * SNOVA_l2 + pass * 2 * SNOVA_l;
                        int32_t pd;
                        memcpy(&pd, pcell, 4);
                        __m256i pq = _mm256_set1_epi32(pd);
                        for (int g = 0; g < 4; g++)
                            qa0[g] = _mm256_dpbusd_avx_epi32(
                                qa0[g], _mm256_load_si256((const __m256i *)(wq + g * 32)), pq);
                        memcpy(&pd, pcell + SNOVA_l, 4);
                        pq = _mm256_set1_epi32(pd);
                        for (int g = 0; g < 4; g++)
                            qa1[g] = _mm256_dpbusd_avx_epi32(
                                qa1[g], _mm256_load_si256((const __m256i *)(wq + g * 32)), pq);
                    }
                    row[pass * 2 + 0] = rct_qv_pack32(
                        rct_q_barrett16(rct_qv_low16(qa0[0], qa0[1])),
                        rct_q_barrett16(rct_qv_low16(qa0[2], qa0[3])));
                    row[pass * 2 + 1] = rct_qv_pack32(
                        rct_q_barrett16(rct_qv_low16(qa1[0], qa1[1])),
                        rct_q_barrett16(rct_qv_low16(qa1[2], qa1[3])));
                }
                rct_qv_ilv32(s0p[0], row[0], row[1]);
                rct_qv_ilv32(s0p[1], row[2], row[3]);
#else
                __m256i acc[SNOVA_l][2];
                for (int i1 = 0; i1 < SNOVA_l; i1++) {
                    acc[i1][0] = _mm256_setzero_si256();
                    acc[i1][1] = _mm256_setzero_si256();
                }
                for (int nj = 0; nj < SNOVA_n; ++nj) {
                    const uint8_t *wb = &rct_qv_wpair[nj * 2 * 2 * RCT_Q_LRP];
                    const gf_t *pcell = pc + nj * SNOVA_l2;
                    for (int h = 0; h < 2; ++h) {
                        __m256i w0 = _mm256_load_si256((const __m256i *)(wb + h * 2 * RCT_Q_LRP));
                        __m256i w1 = _mm256_load_si256((const __m256i *)(wb + h * 2 * RCT_Q_LRP + 32));
                        for (int i1 = 0; i1 < SNOVA_l; i1++) {
                            uint16_t pw;
                            memcpy(&pw, pcell + i1 * SNOVA_l + 2 * h, 2);
                            __m256i pb = _mm256_set1_epi16((short)pw);
                            acc[i1][0] = _mm256_add_epi16(acc[i1][0], _mm256_maddubs_epi16(w0, pb));
                            acc[i1][1] = _mm256_add_epi16(acc[i1][1], _mm256_maddubs_epi16(w1, pb));
                        }
                    }
                }
                for (int hh = 0; hh < 2; ++hh) {
                    __m256i A = rct_qv_pack32(rct_q_barrett16(acc[2 * hh][0]), rct_q_barrett16(acc[2 * hh][1]));
                    __m256i B = rct_qv_pack32(rct_q_barrett16(acc[2 * hh + 1][0]), rct_q_barrett16(acc[2 * hh + 1][1]));
                    rct_qv_ilv32(s0p[hh], A, B);
                }
#endif
#endif
            }
            for (int col = 0; col < SNOVA_l * SNOVA_r; ++col) {
                __m256i v0 = _mm256_setzero_si256();
#if RCT_Q_LR16 == 2
                __m256i v1 = _mm256_setzero_si256();
#endif
                for (int ni = 0; ni < SNOVA_n; ++ni) {
                    const uint8_t *sb = &rct_qv_s0all[ni * 4 * RCT_Q_LRP];
                    const uint16_t *wc = (const uint16_t *)&rct_qv_wpair[ni * 4 * RCT_Q_LRP];
                    __m256i b0 = _mm256_set1_epi16((short)wc[col]);
                    __m256i b1 = _mm256_set1_epi16((short)wc[RCT_Q_LRP + col]);
                    v0 = _mm256_add_epi16(v0, _mm256_maddubs_epi16(
                        _mm256_load_si256((const __m256i *)sb), b0));
                    v0 = _mm256_add_epi16(v0, _mm256_maddubs_epi16(
                        _mm256_load_si256((const __m256i *)(sb + 2 * RCT_Q_LRP)), b1));
#if RCT_Q_LR16 == 2
                    v1 = _mm256_add_epi16(v1, _mm256_maddubs_epi16(
                        _mm256_load_si256((const __m256i *)(sb + 32)), b0));
                    v1 = _mm256_add_epi16(v1, _mm256_maddubs_epi16(
                        _mm256_load_si256((const __m256i *)(sb + 2 * RCT_Q_LRP + 32)), b1));
#endif
                }
                __m256i *s1 = (__m256i *)&sum_t1w[(mi * SNOVA_l * SNOVA_r + col) * RCT_Q_LRP];
                _mm256_store_si256(s1, v0);
#if RCT_Q_LR16 == 2
                _mm256_store_si256(s1 + 1, v1);
#endif
            }
        }
#else
        for (int mi = 0; mi < SNOVA_m1; ++mi)
            for (int ni = 0; ni < SNOVA_n; ++ni) {
                _Alignas(32) uint16_t Prow[SNOVA_n * SNOVA_l2];
                {
                    const gf_t *psrc = &pkx->P[(mi * SNOVA_n + ni) * SNOVA_n * SNOVA_l2];
                    for (int i = 0; i < SNOVA_n * SNOVA_l2; i += 16)
                        _mm256_store_si256((__m256i *)&Prow[i],
                            _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)&psrc[i])));
                }
                _Alignas(32) uint16_t sum_t0[SNOVA_l * RCT_Q_LRP] = {0};
                for (int nj = 0; nj < SNOVA_n; ++nj)
                    for (int i1 = 0; i1 < SNOVA_l; i1++) {
                        __m256i *s0 = (__m256i *)&sum_t0[i1 * RCT_Q_LRP];
                        const __m256i *wp = (const __m256i *)&whipw[(nj * SNOVA_l) * RCT_Q_LRP];
                        const uint16_t *prow =
                            &Prow[nj * SNOVA_l2 + i1 * SNOVA_l];
                        __m256i pv[SNOVA_l];
                        for (int k1 = 0; k1 < SNOVA_l; k1++) pv[k1] = _mm256_set1_epi16((short)prow[k1]);
                        for (int gk = 0; gk < RCT_Q_LR16; ++gk) {
                            __m256i acc = s0[gk];
                            for (int k1 = 0; k1 < SNOVA_l; k1++)
                                acc = _mm256_add_epi16(acc, _mm256_mullo_epi16(pv[k1], wp[k1 * RCT_Q_LR16 + gk]));
                            s0[gk] = acc;
                        }
                    }
                for (int b1 = 0; b1 < SNOVA_l * RCT_Q_LRP; ++b1) sum_t0[b1] %= SNOVA_q;
                for (int a1 = 0; a1 < SNOVA_l; ++a1)
                    for (int k1 = 0; k1 < SNOVA_l; k1++) {
                        const __m256i *s0 = (const __m256i *)&sum_t0[k1 * RCT_Q_LRP];
                        const uint16_t *wrow =
                            &whipw[ni * SNOVA_l * RCT_Q_LRP + k1 * RCT_Q_LRP + a1 * SNOVA_r];
                        for (int i1 = 0; i1 < SNOVA_r; i1++) {
                            __m256i wv = _mm256_set1_epi16((short)wrow[i1]);
                            __m256i *s1 = (__m256i *)&sum_t1w[((mi * SNOVA_l + a1) * SNOVA_r + i1) * RCT_Q_LRP];
                            for (int gk = 0; gk < RCT_Q_LR16; ++gk)
                                s1[gk] = _mm256_add_epi16(s1[gk], _mm256_mullo_epi16(wv, s0[gk]));
                        }
                    }
            }
#endif
}

static void rct_vf_reindex_oddq(rct_vf_ctx *c) {
#if !RCT_Q_EMM
    uint16_t *sum_t1s = c->sum_t1s;
#endif
        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_r * RCT_Q_LRP; i1++) sum_t1w[i1] %= SNOVA_q;
#if RCT_Q_EMM
        {
#define RCT_QV_IL(j2) ((j2) < 2 * SNOVA_r ? (char)((((j2) & 1) ? SNOVA_r : 0) + ((j2) >> 1)) : (char)-1)
            const __m128i ilp = _mm_setr_epi8(
                RCT_QV_IL(0), RCT_QV_IL(1), RCT_QV_IL(2), RCT_QV_IL(3),
                RCT_QV_IL(4), RCT_QV_IL(5), RCT_QV_IL(6), RCT_QV_IL(7),
                RCT_QV_IL(8), RCT_QV_IL(9), RCT_QV_IL(10), RCT_QV_IL(11),
                RCT_QV_IL(12), RCT_QV_IL(13), RCT_QV_IL(14), RCT_QV_IL(15));
#undef RCT_QV_IL
            for (int mi = 0; mi < SNOVA_m1; ++mi)
                for (int a1 = 0; a1 < SNOVA_l; ++a1)
                    for (int i1 = 0; i1 < SNOVA_r; i1++) {
                        const uint16_t *srow = &sum_t1w[((mi * SNOVA_l + a1) * SNOVA_r + i1) * RCT_Q_LRP];
                        __m256i pk = rct_qv_pack32(
                            _mm256_load_si256((const __m256i *)srow),
                            _mm256_load_si256((const __m256i *)(srow + 16)));
                        __m128i x0 = _mm256_castsi256_si128(pk);
                        __m128i x1 = _mm256_extracti128_si256(pk, 1);
                        uint8_t *d0 = &rct_qv_s1p8[(mi * 8 + a1 * 2) * 128 + i1 * 2 * SNOVA_r];
                        _mm_storeu_si128((__m128i *)d0, _mm_shuffle_epi8(x0, ilp));
                        _mm_storeu_si128((__m128i *)(d0 + 128),
                                         _mm_shuffle_epi8(_mm_alignr_epi8(x1, x0, 2 * SNOVA_r), ilp));
                    }
        }
#else
        for (int mi = 0; mi < SNOVA_m1; ++mi)
            for (int a1 = 0; a1 < SNOVA_l; ++a1)
                for (int i1 = 0; i1 < SNOVA_r; i1++)
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        for (int j1 = 0; j1 < SNOVA_r; j1++)
                            sum_t1s[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1] =
                                sum_t1w[((mi * SNOVA_l + a1) * SNOVA_r + i1) * RCT_Q_LRP + b1 * SNOVA_r + j1];
#endif
}

#endif

#endif
