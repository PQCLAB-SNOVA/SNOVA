#ifndef RCT_VERIFY_EMAT_H
#define RCT_VERIFY_EMAT_H

static void rct_vf_emat(rct_vf_ctx *c) {
#if RCT_VF_JOG
    rct_vf_emat_jog(c);
    return;
#endif
    RCT_SCRATCH __attribute__((unused)) _Alignas(64) uint8_t rct_vf_amt[SNOVA_o * SNOVA_alpha * 64];
    RCT_SCRATCH __attribute__((unused)) _Alignas(16) uint8_t rct_vf_q12[SNOVA_o * SNOVA_alpha * 16];
    const rct_pk_t *pkx = c->pkx;
    gf_t *hash_in_GF = c->hash_gf;
#if RCT_VF_EMM
    gf_t *sum_t1q = c->sum_t1q;
#elif !RCT_Q_SIMD
    gf_t *sum_t1 = c->sum_t1;
#endif
#if RCT_Q_SIMD && !RCT_Q_EMM
    uint16_t *sum_t1s = c->sum_t1s;
#endif
#if RCT_Q_SIMD
    uint16_t hash_u16[SNOVA_o * SNOVA_lr + 16] = {0};
#endif
#if RCT_Q_EMM
    for (int t_idx = 0; t_idx < SNOVA_o * SNOVA_alpha; t_idx++)
        rct_qv_tr8(&rct_qv_amt[t_idx * 64], &pkx->Am[t_idx * SNOVA_r2], SNOVA_r);
#endif
#if RCT_VF_EMM
    for (int t_idx = 0; t_idx < SNOVA_o * SNOVA_alpha; t_idx++)
        rct_vf_tr8(&rct_vf_amt[t_idx * 64], &pkx->Am[t_idx * SNOVA_r2], SNOVA_r);
    {
        const __m128i m0f_ = _mm_set1_epi8(0x0f);
        for (int t_idx = 0; t_idx < SNOVA_o * SNOVA_alpha; t_idx++) {
            const gf_t *q1p = &pkx->q1[t_idx * SNOVA_l];
            const gf_t *q2p = &pkx->q2[t_idx * SNOVA_l];
            uint16_t a01, a23;
            memcpy(&a01, q1p, 2);
            memcpy(&a23, q1p + 2, 2);
            unsigned i01 = (a01 & 0xFFu) | (unsigned)((a01 >> 8) << 4);
            unsigned i23 = (a23 & 0xFFu) | (unsigned)((a23 >> 8) << 4);
            int32_t q2w;
            memcpy(&q2w, q2p, 4);
            __m128i q2v = _mm_cvtsi32_si128(q2w);
            __m128i r01 = _mm_shuffle_epi8(_mm_load_si128((const __m128i *)rct_mtk2[i01]), q2v);
            __m128i r23 = _mm_shuffle_epi8(_mm_load_si128((const __m128i *)rct_mtk2[i23]), q2v);
            __m128i p01 = _mm_unpacklo_epi32(_mm_and_si128(r01, m0f_),
                                             _mm_and_si128(_mm_srli_epi16(r01, 4), m0f_));
            __m128i p23 = _mm_unpacklo_epi32(_mm_and_si128(r23, m0f_),
                                             _mm_and_si128(_mm_srli_epi16(r23, 4), m0f_));
            _mm_store_si128((__m128i *)&rct_vf_q12[t_idx * 16], _mm_unpacklo_epi64(p01, p23));
        }
    }
#endif
    for (int mi = 0; mi < SNOVA_o; ++mi) {
#if RCT_Q_EMM
        __m256i qv_h0 = _mm256_setzero_si256(), qv_h1 = _mm256_setzero_si256();
#endif
        for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
            int mi_prime = i_prime(mi, alpha);
            gf_t gf16m_temp1[SNOVA_r2] = {0};
            gf_t gf16m_temp2[SNOVA_lr + 16] = {0};
            (void)gf16m_temp2;
            (void)gf16m_temp1;
#if RCT_USE_SIMD && SNOVA_l == 4 && SNOVA_r2 <= 64
            _Alignas(32) uint8_t t1b[64];
#endif
#if RCT_VF_EMM
            {
                const uint8_t *cc = &rct_vf_q12[(mi * SNOVA_alpha + alpha) * 16];
                __m256i t1lo = _mm256_setzero_si256(), t1hi = _mm256_setzero_si256();
                for (int ab = 0; ab < SNOVA_l2; ab++) {
                    const gf_t *base = &sum_t1q[(mi_prime * SNOVA_l2 + ab) * 64];
                    __m256i cv = RCT_BC(cc[ab]);
                    t1lo = _mm256_xor_si256(t1lo, RCT_SV(cv, _mm256_load_si256((const __m256i *)base)));
                    t1hi = _mm256_xor_si256(t1hi, RCT_SV(cv, _mm256_load_si256((const __m256i *)(base + 32))));
                }
                _mm256_store_si256((__m256i *)t1b, rct_gfni_cleanup256(t1lo));
                _mm256_store_si256((__m256i *)(t1b + 32), rct_gfni_cleanup256(t1hi));
            }
#elif RCT_USE_SIMD && SNOVA_l == 4 && SNOVA_r2 <= 64
            {
#if SNOVA_r2 <= 32
                __m256i t1 = _mm256_setzero_si256();
                for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                    __m256i s = _mm256_setzero_si256();
                    for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                        const gf_t *base = &sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2];
                        __m256i qv = RCT_BC(pkx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
                        s = _mm256_xor_si256(s, RCT_SV(qv, _mm256_loadu_si256((const __m256i *)base)));
                    }
                    __m256i qa = RCT_BC(pkx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]);
                    t1 = _mm256_xor_si256(t1, RCT_SV(qa, rct_gfni_cleanup256(s)));
                }
                _mm256_store_si256((__m256i *)t1b, rct_gfni_cleanup256(t1));
                memcpy(gf16m_temp1, t1b, SNOVA_r2);
#else
                __m256i t1lo = _mm256_setzero_si256(), t1hi = _mm256_setzero_si256();
                for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                    __m256i slo = _mm256_setzero_si256(), shi = _mm256_setzero_si256();
                    for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                        const gf_t *base = &sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2];
                        __m256i qv = RCT_BC(pkx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
                        slo = _mm256_xor_si256(slo, RCT_SV(qv, _mm256_loadu_si256((const __m256i *)base)));
                        shi = _mm256_xor_si256(shi, RCT_SV(qv, _mm256_loadu_si256((const __m256i *)(base + 32))));
                    }
                    slo = rct_gfni_cleanup256(slo);
                    shi = rct_gfni_cleanup256(shi);
                    __m256i qa = RCT_BC(pkx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]);
                    t1lo = _mm256_xor_si256(t1lo, RCT_SV(qa, slo));
                    t1hi = _mm256_xor_si256(t1hi, RCT_SV(qa, shi));
                }
                _mm256_store_si256((__m256i *)t1b, rct_gfni_cleanup256(t1lo));
                _mm256_store_si256((__m256i *)(t1b + 32), rct_gfni_cleanup256(t1hi));
                memcpy(gf16m_temp1, t1b, SNOVA_r2);
#endif
            }
#elif RCT_Q_EMM
            if ((alpha & 1) == 0) {
                const int t0i = mi * SNOVA_alpha + alpha, t1i = t0i + 1;
                const int mp1 = i_prime(mi, alpha + 1);
                _Alignas(32) uint8_t q12b[32];
                {
                    int32_t q1d, q2d;
                    memcpy(&q1d, &pkx->q1[t0i * SNOVA_l], 4);
                    memcpy(&q2d, &pkx->q2[t0i * SNOVA_l], 4);
                    _mm_store_si128((__m128i *)q12b, rct_qv_pack16(rct_q_barrett16(_mm256_mullo_epi16(
                        _mm256_shuffle_epi8(_mm256_set1_epi32(q1d), RCT_QV_PA0),
                        _mm256_shuffle_epi8(_mm256_set1_epi32(q2d), RCT_QV_PB)))));
                    memcpy(&q1d, &pkx->q1[t1i * SNOVA_l], 4);
                    memcpy(&q2d, &pkx->q2[t1i * SNOVA_l], 4);
                    _mm_store_si128((__m128i *)(q12b + 16), rct_qv_pack16(rct_q_barrett16(_mm256_mullo_epi16(
                        _mm256_shuffle_epi8(_mm256_set1_epi32(q1d), RCT_QV_PA0),
                        _mm256_shuffle_epi8(_mm256_set1_epi32(q2d), RCT_QV_PB)))));
                }
                const uint16_t *qp0 = (const uint16_t *)q12b;
                const uint16_t *qp1 = (const uint16_t *)(q12b + 16);
                const uint8_t *sp0 = &rct_qv_s1p8[mi_prime * 8 * 128];
                const uint8_t *sp1 = &rct_qv_s1p8[mp1 * 8 * 128];
                _Alignas(64) uint8_t t1b8[2][64];
                {
                    __m256i ta = _mm256_setzero_si256(), tb = _mm256_setzero_si256();
                    __m256i tc = _mm256_setzero_si256(), td = _mm256_setzero_si256();
                    __m256i ua = _mm256_setzero_si256(), ub = _mm256_setzero_si256();
                    __m256i uc = _mm256_setzero_si256(), ud = _mm256_setzero_si256();
                    for (int p = 0; p < 8; p++) {
                        __m256i bq = _mm256_set1_epi16((short)qp0[p]);
                        ta = _mm256_add_epi16(ta, _mm256_maddubs_epi16(
                            _mm256_load_si256((const __m256i *)(sp0 + p * 128)), bq));
                        tb = _mm256_add_epi16(tb, _mm256_maddubs_epi16(
                            _mm256_load_si256((const __m256i *)(sp0 + p * 128 + 32)), bq));
                        tc = _mm256_add_epi16(tc, _mm256_maddubs_epi16(
                            _mm256_load_si256((const __m256i *)(sp0 + p * 128 + 64)), bq));
                        td = _mm256_add_epi16(td, _mm256_maddubs_epi16(
                            _mm256_load_si256((const __m256i *)(sp0 + p * 128 + 96)), bq));
                        bq = _mm256_set1_epi16((short)qp1[p]);
                        ua = _mm256_add_epi16(ua, _mm256_maddubs_epi16(
                            _mm256_load_si256((const __m256i *)(sp1 + p * 128)), bq));
                        ub = _mm256_add_epi16(ub, _mm256_maddubs_epi16(
                            _mm256_load_si256((const __m256i *)(sp1 + p * 128 + 32)), bq));
                        uc = _mm256_add_epi16(uc, _mm256_maddubs_epi16(
                            _mm256_load_si256((const __m256i *)(sp1 + p * 128 + 64)), bq));
                        ud = _mm256_add_epi16(ud, _mm256_maddubs_epi16(
                            _mm256_load_si256((const __m256i *)(sp1 + p * 128 + 96)), bq));
                    }
                    _mm256_store_si256((__m256i *)t1b8[0],
                                       rct_qv_pack32(rct_q_barrett16(ta), rct_q_barrett16(tb)));
                    _mm256_store_si256((__m256i *)(t1b8[0] + 32),
                                       rct_qv_pack32(rct_q_barrett16(tc), rct_q_barrett16(td)));
                    _mm256_store_si256((__m256i *)t1b8[1],
                                       rct_qv_pack32(rct_q_barrett16(ua), rct_q_barrett16(ub)));
                    _mm256_store_si256((__m256i *)(t1b8[1] + 32),
                                       rct_qv_pack32(rct_q_barrett16(uc), rct_q_barrett16(ud)));
                }
                _Alignas(64) uint8_t t1t8[2][64];
                rct_qv_tr8(t1t8[0], t1b8[0], SNOVA_r);
                rct_qv_tr8(t1t8[1], t1b8[1], SNOVA_r);
                __m256i m0, m1, m2, m3;
                rct_qv_mm_rx4(&m0, &m1, t1t8[0], &pkx->Bm[t0i * SNOVA_lr]);
                rct_qv_mm_rx4(&m2, &m3, t1t8[1], &pkx->Bm[t1i * SNOVA_lr]);
                _Alignas(32) uint8_t t2b8[2][32];
                _mm256_store_si256((__m256i *)t2b8[0],
                                   rct_qv_pack32(rct_q_barrett16(m0), rct_q_barrett16(m1)));
                _mm256_store_si256((__m256i *)t2b8[1],
                                   rct_qv_pack32(rct_q_barrett16(m2), rct_q_barrett16(m3)));
                rct_qv_mm_rx4(&m0, &m1, &rct_qv_amt[t0i * 64], t2b8[0]);
                rct_qv_mm_rx4(&m2, &m3, &rct_qv_amt[t1i * 64], t2b8[1]);
                qv_h0 = _mm256_add_epi16(qv_h0, _mm256_add_epi16(m0, m2));
                qv_h1 = _mm256_add_epi16(qv_h1, _mm256_add_epi16(m1, m3));
            }
#elif RCT_Q_SIMD && (SNOVA_r == SNOVA_l) && RCT_Q_HAVE_MAGIC
            {
                __m256i t1v = _mm256_setzero_si256();
                for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                    __m256i t0v = _mm256_setzero_si256();
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        t0v = _mm256_add_epi16(t0v, _mm256_mullo_epi16(
                            _mm256_set1_epi16((short)pkx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]),
                            _mm256_load_si256((const __m256i *)&sum_t1s[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2])));
                    t0v = rct_q_barrett16(t0v);
                    t1v = _mm256_add_epi16(t1v, _mm256_mullo_epi16(
                        _mm256_set1_epi16((short)pkx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]), t0v));
                }
                t1v = rct_q_barrett16(t1v);
                _Alignas(32) uint16_t t1buf[SNOVA_r2];
                _mm256_store_si256((__m256i *)t1buf, t1v);
                for (int k = 0; k < SNOVA_r2; ++k) gf16m_temp1[k] = (gf_t)t1buf[k];
            }
#elif RCT_Q_SIMD
            {
                uint16_t t0acc[SNOVA_r2];
                uint16_t t1[SNOVA_r2] = {0};
                for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                    for (int i1 = 0; i1 < SNOVA_r2; i1++) t0acc[i1] = 0;
                    for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                        uint16_t qb = pkx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1];
                        const uint16_t *st = &sum_t1s[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2];
                        for (int i1 = 0; i1 < SNOVA_r2; i1++) t0acc[i1] += qb * st[i1];
                    }
                    for (int i1 = 0; i1 < SNOVA_r2; i1++) t0acc[i1] %= SNOVA_q;
                    uint16_t qa = pkx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1];
                    for (int i1 = 0; i1 < SNOVA_r2; i1++) t1[i1] += t0acc[i1] * qa;
                }
                for (int i1 = 0; i1 < SNOVA_r2; i1++) gf16m_temp1[i1] = (gf_t)(t1[i1] % SNOVA_q);
            }
#else
            for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                gf_t sumb[SNOVA_r2] = {0};
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    for (int i1 = 0; i1 < SNOVA_r; i1++)
                        for (int j1 = 0; j1 < SNOVA_r; j1++)
                            gf_set_add(&sumb[i1 * SNOVA_r + j1],
                                       gf_mult(sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1],
                                               pkx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]));
                for (int i1 = 0; i1 < SNOVA_r; i1++)
                    for (int j1 = 0; j1 < SNOVA_r; j1++)
                        gf_set_add(&gf16m_temp1[i1 * SNOVA_r + j1],
                                   gf_mult(sumb[i1 * SNOVA_r + j1], pkx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]));
            }
#endif
#if RCT_VF_EMM
            {
                _Alignas(32) uint8_t t2b[32];
                __m256i t2 = rct_vf_mm_dw(t1b, &pkx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr]);
                _mm256_store_si256((__m256i *)t2b, rct_gfni_cleanup256(t2));
                __m256i h = rct_vf_mm_dw(&rct_vf_amt[(mi * SNOVA_alpha + alpha) * 64], t2b);
                __m256i cur = _mm256_loadu_si256((const __m256i *)&hash_in_GF[mi * SNOVA_lr]);
                _mm256_storeu_si256((__m256i *)&hash_in_GF[mi * SNOVA_lr], _mm256_xor_si256(cur, h));
            }
#elif RCT_USE_SIMD && SNOVA_l == 4
            rct_matmul_l4rows(gf16m_temp2, gf16m_temp1, &pkx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_r, SNOVA_r);
            rct_matmul_l4rows_add(&hash_in_GF[mi * SNOVA_lr], &pkx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2],
                                  gf16m_temp2, SNOVA_r, SNOVA_r);
#elif RCT_Q_EMM
#elif RCT_Q_SIMD && (SNOVA_r == SNOVA_l)
            {
                _Alignas(32) uint16_t t2acc[SNOVA_l2] = {0};
                rct_q_matmul4_add(t2acc, gf16m_temp1, &pkx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr]);
                gf_t t2u8[SNOVA_l2];
                for (int k = 0; k < SNOVA_l2; ++k) t2u8[k] = (gf_t)(t2acc[k] % SNOVA_q);
                rct_q_matmul4_add(&hash_u16[mi * SNOVA_lr], &pkx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], t2u8);
            }
#elif RCT_Q_SIMD
            {
                uint16_t t2[SNOVA_lr] = {0};
                for (int i1 = 0; i1 < SNOVA_r; i1++)
                    for (int j1 = 0; j1 < SNOVA_l; j1++)
                        for (int k1 = 0; k1 < SNOVA_r; k1++)
                            t2[i1 * SNOVA_l + j1] += gf16m_temp1[i1 * SNOVA_r + k1] *
                                pkx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + k1 * SNOVA_l + j1];
                for (int k = 0; k < SNOVA_lr; ++k) t2[k] %= SNOVA_q;
                for (int i1 = 0; i1 < SNOVA_r; i1++)
                    for (int j1 = 0; j1 < SNOVA_l; j1++)
                        for (int k1 = 0; k1 < SNOVA_r; k1++)
                            hash_u16[mi * SNOVA_lr + i1 * SNOVA_l + j1] +=
                                pkx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + k1] *
                                t2[k1 * SNOVA_l + j1];
            }
#else
            gf_mat_mul_add_lr(gf16m_temp2, gf16m_temp1, &pkx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_r, SNOVA_r, SNOVA_l);
            gf_mat_mul_add_lr(&hash_in_GF[mi * SNOVA_lr], &pkx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp2,
                              SNOVA_r, SNOVA_r, SNOVA_l);
#endif
        }
#if RCT_Q_EMM
        {
            uint16_t *hp = &hash_u16[mi * SNOVA_lr];
            _mm256_storeu_si256((__m256i *)hp,
                                _mm256_add_epi16(_mm256_loadu_si256((const __m256i *)hp), qv_h0));
            _mm256_storeu_si256((__m256i *)(hp + 16),
                                _mm256_add_epi16(_mm256_loadu_si256((const __m256i *)(hp + 16)), qv_h1));
        }
#endif
    }
#if RCT_Q_SIMD
    for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) hash_in_GF[i1] = (gf_t)(hash_u16[i1] % SNOVA_q);
#endif
#if RCT_VF_EMM
    for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1 += 32)
        _mm256_storeu_si256((__m256i *)&hash_in_GF[i1],
                            rct_gfni_cleanup256(_mm256_loadu_si256((const __m256i *)&hash_in_GF[i1])));
#endif
}

#endif
