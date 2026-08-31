#ifndef RCT_VERIFY_AQ_H
#define RCT_VERIFY_AQ_H

#if RCT_VF_AQ

static void rct_vf_contract_aq(rct_vf_ctx *c) {
    const rct_pk_t *pkx = c->pkx;
    uint8_t *whipped_sig2 = c->whipped_sig2;
    uint8_t *sum_t1p = c->sum_t1p;
#define RCT_VF_N4 (((SNOVA_n + 3) / 4) * 4)
        RCT_SCRATCH _Alignas(32) uint8_t rct_vf_wpk[RCT_VF_N4 * 2 * 32];
        memset(rct_vf_wpk, 0, sizeof(rct_vf_wpk));
        RCT_SCRATCH _Alignas(32) uint8_t rct_vf_lqt[SNOVA_n * 2 * 128];
        for (int nj = 0; nj < SNOVA_n; ++nj) {
            __m256i r0 = _mm256_load_si256((const __m256i *)&whipped_sig2[(nj * 4 + 0) * SNOVA_lr32]);
            __m256i r1 = _mm256_load_si256((const __m256i *)&whipped_sig2[(nj * 4 + 1) * SNOVA_lr32]);
            __m256i r2 = _mm256_load_si256((const __m256i *)&whipped_sig2[(nj * 4 + 2) * SNOVA_lr32]);
            __m256i r3 = _mm256_load_si256((const __m256i *)&whipped_sig2[(nj * 4 + 3) * SNOVA_lr32]);
            __m256i p0 = _mm256_or_si256(r0, _mm256_slli_epi16(r1, 4));
            __m256i p1 = _mm256_or_si256(r2, _mm256_slli_epi16(r3, 4));
            _mm256_store_si256((__m256i *)&rct_vf_wpk[(nj * 2 + 0) * 32], p0);
            _mm256_store_si256((__m256i *)&rct_vf_wpk[(nj * 2 + 1) * 32], p1);
            __m256i q[4];
            rct_vf_aq_tree(_mm256_shuffle_epi8(p0, rct_vf_aqlo), q);
            _mm256_store_si256((__m256i *)&rct_vf_lqt[(nj * 2 + 0) * 128 + 0], q[0]);
            _mm256_store_si256((__m256i *)&rct_vf_lqt[(nj * 2 + 0) * 128 + 32], q[1]);
            _mm256_store_si256((__m256i *)&rct_vf_lqt[(nj * 2 + 0) * 128 + 64], q[2]);
            _mm256_store_si256((__m256i *)&rct_vf_lqt[(nj * 2 + 0) * 128 + 96], q[3]);
            rct_vf_aq_tree(_mm256_shuffle_epi8(p1, rct_vf_aqlo), q);
            _mm256_store_si256((__m256i *)&rct_vf_lqt[(nj * 2 + 1) * 128 + 0], q[0]);
            _mm256_store_si256((__m256i *)&rct_vf_lqt[(nj * 2 + 1) * 128 + 32], q[1]);
            _mm256_store_si256((__m256i *)&rct_vf_lqt[(nj * 2 + 1) * 128 + 64], q[2]);
            _mm256_store_si256((__m256i *)&rct_vf_lqt[(nj * 2 + 1) * 128 + 96], q[3]);
        }
#if RCT_VERIFY_STREAM
        RCT_SCRATCH _Alignas(32) uint8_t vf_st0p_all[SNOVA_m1 * SNOVA_n * 2 * 32];
        memset(vf_st0p_all, 0, sizeof(vf_st0p_all));
        {
            snova_pgen_t pg;
            snova_pgen_init(&pg, pkx->pk_seed);
            snova_prow_t rw;
            _Alignas(32) uint8_t segbuf[RCT_VF_N4 * SNOVA_l2];
            _Alignas(32) uint8_t qrow[(RCT_VF_N4 / 4) * 128];
            int blk, bmi, bni, bnc;
            (void)bmi; (void)bni; (void)bnc;
            while (snova_pgen_peek(&pg, &blk, &bmi, &bni, &bnc)) {
                const int col0 = (blk == 1) ? SNOVA_v : 0;
                const int qa = col0 & ~3;
                const int pre = col0 - qa;
                (void)snova_pgen_next_row_into(&pg, &rw,
                    (gf_t *)(segbuf + (size_t)pre * SNOVA_l2));
                const int mi = rw.mi;
                const int ni = (rw.block == 2) ? SNOVA_v + rw.ni : rw.ni;
                int ncols = rw.ncols;
                if (rw.block == 2) {
                    memcpy(segbuf + (size_t)SNOVA_v * SNOVA_l2,
                           &pkx->P22[((size_t)mi * SNOVA_o + rw.ni) * SNOVA_o * SNOVA_l2],
                           (size_t)SNOVA_o * SNOVA_l2);
                    ncols = SNOVA_n;
                }
                const int nq = (pre + ncols + 3) >> 2;
                memset(segbuf, 0, (size_t)pre * SNOVA_l2);
                memset(segbuf + (size_t)(pre + ncols) * SNOVA_l2, 0,
                       (size_t)(nq * 4 - pre - ncols) * SNOVA_l2);
                for (int q4 = 0; q4 < nq; ++q4)
                    rct_vf_aq_quad(segbuf + q4 * 64, qrow + q4 * 128);
                __m256i a01a = _mm256_setzero_si256(), a23a = _mm256_setzero_si256();
                __m256i a01b = _mm256_setzero_si256(), a23b = _mm256_setzero_si256();
                for (int njq = 0; njq < nq * 4; njq += 4) {
                    const uint8_t *qb = qrow + (njq >> 2) * 128;
                    const uint8_t *wb = &rct_vf_wpk[(qa + njq) * 64];
                    for (int c = 0; c < 4; c += 2) {
                        const int cb = (c >> 1) * 64;
                        __m256i wv0 = _mm256_load_si256((const __m256i *)(wb + c * 64));
                        __m256i wv1 = _mm256_load_si256((const __m256i *)(wb + c * 64 + 32));
                        __m256i xv0 = _mm256_load_si256((const __m256i *)(wb + c * 64 + 64));
                        __m256i xv1 = _mm256_load_si256((const __m256i *)(wb + c * 64 + 96));
                        a01a = _mm256_xor_si256(a01a, _mm256_gf2p8affine_epi64_epi8(
                                    wv0, rct_vf_aq_bq(qb, cb + 0), 0));
                        a23a = _mm256_xor_si256(a23a, _mm256_gf2p8affine_epi64_epi8(
                                    wv0, rct_vf_aq_bq(qb, cb + 8), 0));
                        a01a = _mm256_xor_si256(a01a, _mm256_gf2p8affine_epi64_epi8(
                                    wv1, rct_vf_aq_bq(qb, cb + 32), 0));
                        a23a = _mm256_xor_si256(a23a, _mm256_gf2p8affine_epi64_epi8(
                                    wv1, rct_vf_aq_bq(qb, cb + 40), 0));
                        a01b = _mm256_xor_si256(a01b, _mm256_gf2p8affine_epi64_epi8(
                                    xv0, rct_vf_aq_bq(qb, cb + 16), 0));
                        a23b = _mm256_xor_si256(a23b, _mm256_gf2p8affine_epi64_epi8(
                                    xv0, rct_vf_aq_bq(qb, cb + 24), 0));
                        a01b = _mm256_xor_si256(a01b, _mm256_gf2p8affine_epi64_epi8(
                                    xv1, rct_vf_aq_bq(qb, cb + 48), 0));
                        a23b = _mm256_xor_si256(a23b, _mm256_gf2p8affine_epi64_epi8(
                                    xv1, rct_vf_aq_bq(qb, cb + 56), 0));
                    }
                }
                uint8_t *sp = &vf_st0p_all[((size_t)mi * SNOVA_n + ni) * 64];
                _mm256_store_si256((__m256i *)sp,
                    _mm256_xor_si256(_mm256_load_si256((const __m256i *)sp),
                                     _mm256_xor_si256(a01a, a01b)));
                _mm256_store_si256((__m256i *)(sp + 32),
                    _mm256_xor_si256(_mm256_load_si256((const __m256i *)(sp + 32)),
                                     _mm256_xor_si256(a23a, a23b)));
                if (rw.block == 2 && rw.ni == SNOVA_o - 1) {
                    for (int h = 0; h < 2; ++h) {
                        __m256i acc[SNOVA_r];
                        for (int p = 0; p < SNOVA_r; p++) acc[p] = _mm256_setzero_si256();
                        for (int nn = 0; nn < SNOVA_n; ++nn)
                            for (int kp = 0; kp < 2; ++kp) {
                                __m256i s0 = _mm256_load_si256((const __m256i *)
                                    &vf_st0p_all[(((size_t)mi * SNOVA_n + nn) * 2 + kp) * 32]);
                                const uint8_t *lb = &rct_vf_lqt[(nn * 2 + kp) * 128];
                                for (int p = 0; p < SNOVA_r; p++)
                                    acc[p] = _mm256_xor_si256(acc[p], _mm256_gf2p8affine_epi64_epi8(
                                                 s0, rct_vf_aq_bq(lb, RCT_VF_AQ_LOFF(h * SNOVA_r + p)), 0));
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
            }
        }
#else
        for (int mi = 0; mi < SNOVA_m1; ++mi) {
            RCT_SCRATCH _Alignas(32) uint8_t rct_vf_st0p[SNOVA_n * 2 * 32];
            _Alignas(32) uint8_t qrow[(RCT_VF_N4 / 4) * 128];
            for (int ni = 0; ni < SNOVA_n; ++ni) {
                const uint8_t *prow = &pkx->P[(size_t)(mi * SNOVA_n + ni) * SNOVA_n * SNOVA_l2];
                for (int q4 = 0; q4 < RCT_VF_N4 / 4; ++q4)
                    rct_vf_aq_quad(prow + q4 * 64, qrow + q4 * 128);
                __m256i a01a = _mm256_setzero_si256(), a23a = _mm256_setzero_si256();
                __m256i a01b = _mm256_setzero_si256(), a23b = _mm256_setzero_si256();
                for (int njq = 0; njq < RCT_VF_N4; njq += 4) {
                    const uint8_t *qb = qrow + (njq >> 2) * 128;
                    const uint8_t *wb = &rct_vf_wpk[njq * 64];
                    for (int c = 0; c < 4; c += 2) {
                        const int cb = (c >> 1) * 64;
                        __m256i wv0 = _mm256_load_si256((const __m256i *)(wb + c * 64));
                        __m256i wv1 = _mm256_load_si256((const __m256i *)(wb + c * 64 + 32));
                        __m256i xv0 = _mm256_load_si256((const __m256i *)(wb + c * 64 + 64));
                        __m256i xv1 = _mm256_load_si256((const __m256i *)(wb + c * 64 + 96));
                        a01a = _mm256_xor_si256(a01a, _mm256_gf2p8affine_epi64_epi8(
                                    wv0, rct_vf_aq_bq(qb, cb + 0), 0));
                        a23a = _mm256_xor_si256(a23a, _mm256_gf2p8affine_epi64_epi8(
                                    wv0, rct_vf_aq_bq(qb, cb + 8), 0));
                        a01a = _mm256_xor_si256(a01a, _mm256_gf2p8affine_epi64_epi8(
                                    wv1, rct_vf_aq_bq(qb, cb + 32), 0));
                        a23a = _mm256_xor_si256(a23a, _mm256_gf2p8affine_epi64_epi8(
                                    wv1, rct_vf_aq_bq(qb, cb + 40), 0));
                        a01b = _mm256_xor_si256(a01b, _mm256_gf2p8affine_epi64_epi8(
                                    xv0, rct_vf_aq_bq(qb, cb + 16), 0));
                        a23b = _mm256_xor_si256(a23b, _mm256_gf2p8affine_epi64_epi8(
                                    xv0, rct_vf_aq_bq(qb, cb + 24), 0));
                        a01b = _mm256_xor_si256(a01b, _mm256_gf2p8affine_epi64_epi8(
                                    xv1, rct_vf_aq_bq(qb, cb + 48), 0));
                        a23b = _mm256_xor_si256(a23b, _mm256_gf2p8affine_epi64_epi8(
                                    xv1, rct_vf_aq_bq(qb, cb + 56), 0));
                    }
                }
                _mm256_store_si256((__m256i *)&rct_vf_st0p[(ni * 2 + 0) * 32],
                                   _mm256_xor_si256(a01a, a01b));
                _mm256_store_si256((__m256i *)&rct_vf_st0p[(ni * 2 + 1) * 32],
                                   _mm256_xor_si256(a23a, a23b));
            }
            for (int h = 0; h < 2; ++h) {
                __m256i acc[SNOVA_r];
                for (int p = 0; p < SNOVA_r; p++) acc[p] = _mm256_setzero_si256();
                for (int ni = 0; ni < SNOVA_n; ++ni)
                    for (int kp = 0; kp < 2; ++kp) {
                        __m256i s0 = _mm256_load_si256(
                            (const __m256i *)&rct_vf_st0p[(ni * 2 + kp) * 32]);
                        const uint8_t *lb = &rct_vf_lqt[(ni * 2 + kp) * 128];
                        for (int p = 0; p < SNOVA_r; p++)
                            acc[p] = _mm256_xor_si256(acc[p], _mm256_gf2p8affine_epi64_epi8(
                                         s0, rct_vf_aq_bq(lb, RCT_VF_AQ_LOFF(h * SNOVA_r + p)), 0));
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
#endif
#if defined(RCT_AQ_SELFTEST) && !RCT_VERIFY_STREAM
        {
            static int aq_st_done = 0;
            if (!aq_st_done) {
                aq_st_done = 1;
                int bad = 0;
                RCT_SCRATCH gf_t st0r[SNOVA_n * SNOVA_l][SNOVA_lr];
                RCT_SCRATCH gf_t st1r[SNOVA_l * SNOVA_r][SNOVA_lr];
                for (int mi = 0; mi < SNOVA_m1; ++mi) {
                    memset(st0r, 0, sizeof(st0r));
                    memset(st1r, 0, sizeof(st1r));
                    for (int ni = 0; ni < SNOVA_n; ++ni)
                        for (int nj = 0; nj < SNOVA_n; ++nj)
                            for (int i1 = 0; i1 < SNOVA_l; ++i1)
                                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                    gf_t s = pkx->P[((size_t)(mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + k1];
                                    for (int b = 0; b < SNOVA_lr; ++b)
                                        st0r[ni * SNOVA_l + i1][b] ^=
                                            rct_multtab[s * SNOVA_q + whipped_sig2[(nj * SNOVA_l + k1) * SNOVA_lr32 + b]];
                                }
                    for (int ni = 0; ni < SNOVA_n; ++ni)
                        for (int k1 = 0; k1 < SNOVA_l; ++k1)
                            for (int j = 0; j < SNOVA_l * SNOVA_r; ++j) {
                                gf_t s = whipped_sig2[(ni * SNOVA_l + k1) * SNOVA_lr32 + j];
                                for (int b = 0; b < SNOVA_lr; ++b)
                                    st1r[j][b] ^= rct_multtab[s * SNOVA_q + st0r[ni * SNOVA_l + k1][b]];
                            }
                    for (int j = 0; j < SNOVA_l * SNOVA_r; ++j)
                        for (int b = 0; b < SNOVA_lr; ++b)
                            if (st1r[j][b] != sum_t1p[(mi * SNOVA_l * SNOVA_r + j) * SNOVA_lr32 + b]) bad++;
                }
                fprintf(stderr, "[AQ-SELFTEST] mismatches=%d\n", bad);
            }
        }
#endif
}

#endif

#endif
