#ifndef RCT_SIGN_ODDQ_H
#define RCT_SIGN_ODDQ_H

#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC
static void rct_sign_apply_t12_oddq(rct_sign_ctx *c, const uint16_t *sol16) {
    gf_t *signature_in_GF = c->signature_in_GF;
    const gf_t *T12 = c->T12;
    _Static_assert((uint32_t)SNOVA_o * SNOVA_l * (SNOVA_q - 1u) * (SNOVA_q - 1u)
                       + SNOVA_q < 65536u,
                   "odd-q SIMD T12-apply row accumulation overflows uint16");
    for (int index = 0; index < SNOVA_v; ++index) {
        gf_t *sigrow = &signature_in_GF[index * SNOVA_lr];
        for (int i1 = 0; i1 < SNOVA_l; ++i1) {
            __m128i acc = _mm_cvtepu8_epi16(
                _mm_loadl_epi64((const __m128i *)&sigrow[i1 * SNOVA_r]));
            const gf_t *tb = &T12[(index * SNOVA_o) * SNOVA_l2 + i1 * SNOVA_l];
            for (int mi = 0; mi < SNOVA_o; ++mi, tb += SNOVA_l2)
                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                    __m128i bv = _mm_set1_epi16((short)tb[k1]);
                    __m128i cv = _mm_loadu_si128(
                        (const __m128i *)&sol16[mi * SNOVA_lr + k1 * SNOVA_r]);
                    acc = _mm_add_epi16(acc, _mm_mullo_epi16(bv, cv));
                }
            __m128i quo = _mm_srli_epi16(
                _mm_mulhi_epu16(acc, _mm_set1_epi16((short)RCT_Q_MAGIC_M)), RCT_Q_MAGIC_S);
            acc = _mm_sub_epi16(acc, _mm_mullo_epi16(quo, _mm_set1_epi16((short)SNOVA_q)));
            uint8_t out8[8];
            _mm_storel_epi64((__m128i *)out8, _mm_packus_epi16(acc, acc));
            memcpy(&sigrow[i1 * SNOVA_r], out8, SNOVA_r);
        }
    }
}
#endif

#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC
static int rct_sign_gauss_oddq_magic(rct_sign_ctx *c, uint16_t (*gu)[RCT_GPAD],
                                     uint16_t *sol16, gf_t *solution) {
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
    int flag_redo;
    {
        enum { RCT_GNB = RCT_GPAD / 16,
               RCT_GW = (SNOVA_q > 23) ? 32 : 64 };
        _Static_assert(2u * ((uint32_t)RCT_GW * (SNOVA_q - 1u) * (SNOVA_q - 1u) + SNOVA_q) < 65536u,
                       "u16 accumulation overflow guard");
        for (int r = 0; r < RCT_GN; ++r) {
            for (int c = 0; c < RCT_GPAD; ++c) gu[r][c] = 0;
            for (int c = 0; c <= RCT_GN; ++c) gu[r][c] = gauss[r][c];
        }
        uint32_t redo_acc = 0;
        for (int i = 0; i < RCT_GN; ++i) {
            const int b0 = i / 16;
            __m256i rowi[RCT_GNB];
            for (int b = b0; b < RCT_GNB; ++b)
                rowi[b] = _mm256_load_si256((const __m256i *)&gu[i][b * 16]);
            uint16_t ii = gu[i][i];
            for (int j = i + 1; j < RCT_GN; ++j) {
                uint32_t need = (1u - ct_gf_nz((uint32_t)ii)) & ct_gf_nz((uint32_t)gu[j][i]);
                __m256i mv = _mm256_set1_epi16((short)(0u - need));
                for (int b = b0; b < RCT_GNB; ++b)
                    rowi[b] = _mm256_add_epi16(rowi[b],
                        _mm256_and_si256(_mm256_load_si256((const __m256i *)&gu[j][b * 16]), mv));
                ii = (uint16_t)(ii + (uint16_t)((0u - need) & (uint32_t)gu[j][i]));
            }
            redo_acc |= (1u - ct_gf_nz((uint32_t)ii));
            uint16_t t_inv = gf_inv_sec((gf_t)ii);
            __m256i tv = _mm256_set1_epi16((short)t_inv);
            for (int b = b0; b < RCT_GNB; ++b) {
                __m256i x = rct_q_barrett16(rowi[b]);
                x = rct_q_barrett16(_mm256_mullo_epi16(x, tv));
                rowi[b] = x;
                _mm256_store_si256((__m256i *)&gu[i][b * 16], x);
            }
            for (int j = i + 1; j < RCT_GN; ++j) {
                uint32_t nz = ct_gf_nz((uint32_t)gu[j][i]);
                uint16_t gji = (uint16_t)(((uint32_t)(SNOVA_q - gu[j][i])) & (0u - nz));
                __m256i gv = _mm256_set1_epi16((short)gji);
                for (int b = b0; b < RCT_GNB; ++b) {
                    __m256i d = _mm256_load_si256((const __m256i *)&gu[j][b * 16]);
                    _mm256_store_si256((__m256i *)&gu[j][b * 16],
                        _mm256_add_epi16(d, _mm256_mullo_epi16(rowi[b], gv)));
                }
            }
            if (!(i % RCT_GW)) {
                for (int j = i + 1; j < RCT_GN; ++j)
                    for (int b = b0; b < RCT_GNB; ++b)
                        _mm256_store_si256((__m256i *)&gu[j][b * 16],
                            rct_q_barrett16(_mm256_load_si256((const __m256i *)&gu[j][b * 16])));
            } else {
                for (int j = i + 1; j < RCT_GN; ++j)
                    gu[j][i + 1] = (uint16_t)(gu[j][i + 1] % SNOVA_q);
            }
        }
        flag_redo = (int)redo_acc;
        SNOVA_CT_DECLASSIFY(&flag_redo, sizeof flag_redo);
        if (!flag_redo) {
            _Static_assert((uint32_t)RCT_GN * (SNOVA_q - 1u) * (SNOVA_q - 1u) < 65536u,
                           "odd-q SIMD backsub horizontal sum overflows uint16");
            for (int c = 0; c < RCT_GPAD; ++c) sol16[c] = 0;
            for (int i = RCT_GN - 1; i >= 0; --i) {
                __m256i acc = _mm256_setzero_si256();
                for (int b = (i + 1) / 16; b < RCT_GNB; ++b)
                    acc = _mm256_add_epi16(acc,
                        _mm256_mullo_epi16(_mm256_load_si256((const __m256i *)&gu[i][b * 16]),
                                           _mm256_load_si256((const __m256i *)&sol16[b * 16])));
                __m128i s = _mm_add_epi16(_mm256_castsi256_si128(acc),
                                          _mm256_extracti128_si256(acc, 1));
                s = _mm_add_epi16(s, _mm_srli_si128(s, 8));
                s = _mm_add_epi16(s, _mm_srli_si128(s, 4));
                s = _mm_add_epi16(s, _mm_srli_si128(s, 2));
                uint16_t sum = (uint16_t)((uint32_t)(uint16_t)_mm_extract_epi16(s, 0) % SNOVA_q);
                sol16[i] = (uint16_t)(((uint32_t)gu[i][RCT_GN] + SNOVA_q - sum) % SNOVA_q);
            }
            for (int c = 0; c < RCT_GN; ++c) solution[c] = (gf_t)sol16[c];
        }
    }
    return flag_redo;
}
#endif

#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && !RCT_Q_HAVE_MAGIC
static int rct_sign_gauss_oddq_nomagic(rct_sign_ctx *c) {
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
    int flag_redo;
    {
        enum { RCT_OLR = SNOVA_o * SNOVA_lr, RCT_OLR_PAD = (RCT_OLR / 16 + 1) * 16 };
        _Alignas(32) static uint16_t gu[RCT_OLR][RCT_OLR_PAD];
        for (int r = 0; r < RCT_OLR; ++r) {
            for (int c = 0; c < RCT_OLR_PAD; ++c) gu[r][c] = 0;
            for (int c = 0; c <= RCT_OLR; ++c) gu[r][c] = gauss[r][c];
        }
        uint32_t redo_acc = 0;
        for (int i = 0; i < RCT_OLR; ++i) {
            for (int j = i + 1; j < RCT_OLR; ++j) {
                uint32_t need = (1u - ct_gf_nz((uint32_t)gu[i][i])) & ct_gf_nz((uint32_t)gu[j][i]);
                uint16_t m = (uint16_t)(0u - need);
                for (int k = i; k < RCT_OLR_PAD; ++k)
                    gu[i][k] = (uint16_t)((gu[i][k] + (gu[j][k] & m)) % SNOVA_q);
            }
            redo_acc |= (1u - ct_gf_nz((uint32_t)gu[i][i]));
            uint16_t t_inv = gf_inv_sec((gf_t)gu[i][i]);
            for (int k = i; k < RCT_OLR_PAD; ++k)
                gu[i][k] = (uint16_t)(((uint32_t)gu[i][k] * t_inv) % SNOVA_q);
            for (int j = i + 1; j < RCT_OLR; ++j) {
                uint32_t nz = ct_gf_nz((uint32_t)gu[j][i]);
                uint16_t gji = (uint16_t)(((uint32_t)(SNOVA_q - gu[j][i])) & (0u - nz));
                for (int k = i; k < RCT_OLR_PAD; ++k)
                    gu[j][k] = (uint16_t)((gu[j][k] + gu[i][k] * gji) % SNOVA_q);
            }
        }
        flag_redo = (int)redo_acc;
        SNOVA_CT_DECLASSIFY(&flag_redo, sizeof flag_redo);
        if (!flag_redo)
            for (int r = 0; r < RCT_OLR; ++r)
                for (int c = r; c <= RCT_OLR; ++c) gauss[r][c] = (gf_t)gu[r][c];
        SNOVA_CLEAR_OBJ(gu);
    }
    return flag_redo;
}
#endif

#if RCT_Q_SIMD && (SNOVA_r == SNOVA_l)
static void rct_sign_wf_oddq_sq(rct_sign_ctx *c) {
    const gf_t *F21 = c->F21, *F12 = c->F12, *whipped_sig = c->whipped_sig;
    gf_t *whipped_F21 = c->whipped_F21, *whipped_F12 = c->whipped_F12;
    _Static_assert((uint32_t)SNOVA_v * SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "u16 accumulation overflow guard");
    {
        static uint16_t wF21u[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr];
        static uint16_t wF12u[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr];
        memset(wF21u, 0, sizeof(wF21u));
        memset(wF12u, 0, sizeof(wF12u));
        for (int mi = 0; mi < SNOVA_m1; mi++)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int idx = 0; idx < SNOVA_o; idx++)
                    for (int nj = 0; nj < SNOVA_v; ++nj) {
                        rct_q_matmul4_add(&wF21u[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr],
                                          &F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2],
                                          &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr]);
                        rct_q_matmul4T_add(&wF12u[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr],
                                           &F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2],
                                           &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr]);
                    }
        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr; i1++)
            whipped_F21[i1] = (gf_t)(wF21u[i1] % SNOVA_q);
        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr; i1++)
            whipped_F12[i1] = (gf_t)(wF12u[i1] % SNOVA_q);
        SNOVA_CLEAR_OBJ(wF21u);
        SNOVA_CLEAR_OBJ(wF12u);
    }
}
#endif

#if RCT_Q_SIMD && (SNOVA_r == SNOVA_l)
static void rct_sign_gauss_scatter_oddq_sq(rct_sign_ctx *c) {
    const gf_t *whipped_F21 = c->whipped_F21, *whipped_F12 = c->whipped_F12;
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm, *Q1 = c->Q1, *Q2 = c->Q2;
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
    _Static_assert((uint32_t)2u * SNOVA_alpha * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "u16 accumulation overflow guard");
    {
        _Alignas(32) static uint16_t gauss16[SNOVA_o * SNOVA_lr][SNOVA_o * SNOVA_lr];
        memset(gauss16, 0, sizeof(gauss16));

        for (int mi = 0; mi < SNOVA_o; mi++)
            for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
                int mi_prime = i_prime(mi, alpha);
                _Alignas(32) uint16_t t0[SNOVA_o * SNOVA_l2] = {0};
                _Alignas(32) uint16_t t1[SNOVA_o * SNOVA_l2] = {0};
                _Alignas(32) uint16_t t2[SNOVA_o * SNOVA_l2] = {0};
                for (int idx = 0; idx < SNOVA_o; idx++)
                    for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                        __m256i qv = _mm256_set1_epi16((short)q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
                        __m256i wf = _mm256_cvtepu8_epi16(_mm_loadu_si128(
                            (const __m128i *)&whipped_F21[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr]));
                        __m256i *d = (__m256i *)&t0[idx * SNOVA_l2];
                        *d = _mm256_add_epi16(*d, _mm256_mullo_epi16(qv, wf));
                    }
                for (int i = 0; i < SNOVA_o * SNOVA_l2; i++) t0[i] %= SNOVA_q;
                for (int idx = 0; idx < SNOVA_o; idx++)
                    for (int i1 = 0; i1 < SNOVA_l; i1++)
                        for (int j1 = 0; j1 < SNOVA_l; j1++)
                            for (int k1 = 0; k1 < SNOVA_r; k1++)
                                t1[idx * SNOVA_l2 + i1 * SNOVA_l + j1] +=
                                    t0[idx * SNOVA_lr + i1 * SNOVA_r + k1] *
                                    Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + k1 * SNOVA_l + j1];
                for (int i = 0; i < SNOVA_o * SNOVA_l2; i++) t1[i] %= SNOVA_q;
                for (int idx = 0; idx < SNOVA_o; idx++)
                    for (int i1 = 0; i1 < SNOVA_l; i1++)
                        for (int j1 = 0; j1 < SNOVA_l; j1++)
                            for (int k1 = 0; k1 < SNOVA_l; k1++)
                                t2[idx * SNOVA_l2 + i1 * SNOVA_l + j1] +=
                                    Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + i1 * SNOVA_l + k1] *
                                    t1[idx * SNOVA_l2 + k1 * SNOVA_l + j1];
                for (int i = 0; i < SNOVA_o * SNOVA_l2; i++) t2[i] %= SNOVA_q;
                for (int idx = 0; idx < SNOVA_o; idx++) {
                    __m256i gfm[4], am[4];
                    uint16_t *gfm16 = (uint16_t *)gfm;
                    uint16_t *am16 = (uint16_t *)am;
                    for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
                        for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                            for (int tj2 = 0; tj2 < SNOVA_l; tj2++)
                                gfm16[ti2 * SNOVA_l2 + tj1 * SNOVA_l + tj2] =
                                    t2[idx * SNOVA_l2 + tj1 * SNOVA_l + ti2];
                    for (int ti1 = 0; ti1 < SNOVA_l; ti1++)
                        for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                            for (int tj2 = 0; tj2 < SNOVA_l; tj2++)
                                am16[ti1 * SNOVA_l2 + tj1 * SNOVA_l + tj2] =
                                    Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r + tj2];
                    for (int ti1 = 0; ti1 < SNOVA_l; ti1++)
                        for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                            __m256i *g = (__m256i *)&gauss16[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx * SNOVA_lr];
                            *g = _mm256_add_epi16(*g, _mm256_mullo_epi16(gfm[ti2], am[ti1]));
                        }
                }
            }

        for (int mi = 0; mi < SNOVA_o; mi++)
            for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
                int mi_prime = i_prime(mi, alpha);
                _Alignas(32) uint16_t t0[SNOVA_o * SNOVA_lr] = {0};
                _Alignas(32) uint16_t t1[SNOVA_o * SNOVA_lr] = {0};
                _Alignas(32) uint16_t t2[SNOVA_o * SNOVA_lr] = {0};
                for (int idx = 0; idx < SNOVA_o; idx++)
                    for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                        __m256i qv = _mm256_set1_epi16((short)q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
                        __m256i wf = _mm256_cvtepu8_epi16(_mm_loadu_si128(
                            (const __m128i *)&whipped_F12[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr]));
                        __m256i *d = (__m256i *)&t0[idx * SNOVA_l2];
                        *d = _mm256_add_epi16(*d, _mm256_mullo_epi16(qv, wf));
                    }
                for (int i = 0; i < SNOVA_o * SNOVA_lr; i++) t0[i] %= SNOVA_q;
                for (int idx = 0; idx < SNOVA_o; idx++)
                    for (int i1 = 0; i1 < SNOVA_r; i1++)
                        for (int j1 = 0; j1 < SNOVA_l; j1++)
                            for (int k1 = 0; k1 < SNOVA_r; k1++)
                                t1[idx * SNOVA_lr + i1 * SNOVA_l + j1] +=
                                    Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + k1] *
                                    t0[idx * SNOVA_lr + j1 * SNOVA_r + k1];
                for (int i = 0; i < SNOVA_o * SNOVA_lr; i++) t1[i] %= SNOVA_q;
                for (int idx = 0; idx < SNOVA_o; idx++)
                    for (int i1 = 0; i1 < SNOVA_r; i1++)
                        for (int j1 = 0; j1 < SNOVA_l; j1++)
                            for (int k1 = 0; k1 < SNOVA_l; k1++)
                                t2[idx * SNOVA_lr + i1 * SNOVA_l + j1] +=
                                    t1[idx * SNOVA_lr + i1 * SNOVA_l + k1] *
                                    Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + k1 * SNOVA_l + j1];
                for (int i = 0; i < SNOVA_o * SNOVA_lr; i++) t2[i] %= SNOVA_q;
                for (int idx = 0; idx < SNOVA_o; idx++) {
                    __m256i gfm[4], bm[4];
                    uint16_t *gfm16 = (uint16_t *)gfm;
                    uint16_t *bm16 = (uint16_t *)bm;
                    for (int ti1 = 0; ti1 < SNOVA_l; ti1++)
                        for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                            for (int tj2 = 0; tj2 < SNOVA_l; tj2++)
                                gfm16[ti1 * SNOVA_l2 + tj1 * SNOVA_l + tj2] =
                                    t2[idx * SNOVA_l2 + ti1 * SNOVA_l + tj1];
                    for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
                        for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                            for (int tj2 = 0; tj2 < SNOVA_l; tj2++)
                                bm16[ti2 * SNOVA_l2 + tj1 * SNOVA_l + tj2] =
                                    Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2];
                    for (int ti1 = 0; ti1 < SNOVA_l; ti1++)
                        for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                            __m256i *g = (__m256i *)&gauss16[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx * SNOVA_lr];
                            *g = _mm256_add_epi16(*g, _mm256_mullo_epi16(gfm[ti1], bm[ti2]));
                        }
                }
            }

        for (int ti = 0; ti < SNOVA_o * SNOVA_lr; ti++)
            for (int tj = 0; tj < SNOVA_o * SNOVA_lr; tj++)
                gauss[ti][tj] = (gf_t)((gauss[ti][tj] + gauss16[ti][tj]) % SNOVA_q);
        SNOVA_CLEAR_OBJ(gauss16);
    }
}
#endif

#if RCT_Q_SIMD && (SNOVA_r == SNOVA_l) && RCT_Q_HAVE_MAGIC && !defined(RCT_FVV_SCALAR)
static void rct_sign_fvv_oddq_sq(rct_sign_ctx *c) {
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm, *sum_t1 = c->sum_t1;
    gf_t *Fvv_in_GF16Matrix = c->Fvv;
    _Static_assert((uint32_t)SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "odd-q Fvv temp1 inner (s) uint16 accumulation may overflow");
    _Static_assert((uint32_t)SNOVA_alpha * SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "odd-q Fvv cross-alpha uint16 accumulation may overflow");
    {
        __m256i Fvvacc[SNOVA_o];
        for (int mi = 0; mi < SNOVA_o; ++mi) Fvvacc[mi] = _mm256_setzero_si256();
        for (int mi = 0; mi < SNOVA_o; ++mi) {
            for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
                int mi_prime = i_prime(mi, alpha);
                const gf_t *q1r = &q1[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                const gf_t *q2r = &q2[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                __m256i t1 = _mm256_setzero_si256();
                for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                    __m256i s = _mm256_setzero_si256();
                    for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                        __m256i blk = _mm256_cvtepu8_epi16(_mm_loadu_si128(
                            (const __m128i *)&sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2]));
                        s = _mm256_add_epi16(s, _mm256_mullo_epi16(_mm256_set1_epi16((short)q2r[b1]), blk));
                    }
#if SNOVA_l * SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) * (SNOVA_q - 1) >= 65536
                    s = rct_q_barrett16(s);
#endif
                    t1 = _mm256_add_epi16(t1, _mm256_mullo_epi16(_mm256_set1_epi16((short)q1r[a1]), s));
                }
                t1 = rct_q_barrett16(t1);
                __m256i bmv = _mm256_cvtepu8_epi16(_mm_loadu_si128(
                    (const __m128i *)&Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr]));
                __m256i t2 = rct_q_barrett16(rct_q_mm4_u16(t1, bmv));
                __m256i amv = _mm256_cvtepu8_epi16(_mm_loadu_si128(
                    (const __m128i *)&Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2]));
                Fvvacc[mi] = _mm256_add_epi16(Fvvacc[mi], rct_q_mm4_u16(amv, t2));
            }
            _Alignas(32) uint16_t fb[16];
            _mm256_store_si256((__m256i *)fb, rct_q_barrett16(Fvvacc[mi]));
            for (int k = 0; k < SNOVA_lr; ++k)
                Fvv_in_GF16Matrix[mi * SNOVA_lr + k] = (gf_t)fb[k];
        }
    }
}
#endif

#if RCT_Q_SIMD && (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC
static void rct_sign_fvv_oddq_rect(rct_sign_ctx *c) {
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm;
#if !RCT_OQDF
    const gf_t *sum_t1 = c->sum_t1;
#endif
    gf_t *Fvv_in_GF16Matrix = c->Fvv;
    {
        enum { RCT_R2P = ((SNOVA_r2 + 15) / 16) * 16 };
        _Alignas(32) static uint16_t Fvvacc[SNOVA_o * SNOVA_lr];
        memset(Fvvacc, 0, sizeof(Fvvacc));
        for (int mi = 0; mi < SNOVA_o; ++mi)
            for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
                int mi_prime = i_prime(mi, alpha);
                const gf_t *q1r = &q1[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                const gf_t *q2r = &q2[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                const gf_t *Bmr = &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr];
                const gf_t *Amr = &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2];
                _Alignas(32) uint16_t t1[RCT_R2P] = {0};
                for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                    __m256i s[RCT_R2P / 16];
                    for (int c = 0; c < RCT_R2P / 16; ++c) s[c] = _mm256_setzero_si256();
                    for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                        __m256i qv = _mm256_set1_epi16((short)q2r[b1]);
#if RCT_OQDF
                        const uint16_t *base = &rct_oq_sum_t1u[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2];
                        for (int c = 0; c < RCT_R2P / 16; ++c)
                            s[c] = _mm256_add_epi16(s[c], _mm256_mullo_epi16(qv,
                                _mm256_loadu_si256((const __m256i *)(base + c * 16))));
#else
                        const gf_t *base = &sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2];
                        for (int c = 0; c < RCT_R2P / 16; ++c)
                            s[c] = _mm256_add_epi16(s[c], _mm256_mullo_epi16(qv,
                                _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(base + c * 16)))));
#endif
                    }
                    __m256i qa = _mm256_set1_epi16((short)q1r[a1]);
                    for (int c = 0; c < RCT_R2P / 16; ++c) {
                        __m256i sc = rct_q_barrett16(s[c]);
                        __m256i tc = _mm256_load_si256((const __m256i *)&t1[c * 16]);
                        _mm256_store_si256((__m256i *)&t1[c * 16],
                                           _mm256_add_epi16(tc, _mm256_mullo_epi16(qa, sc)));
                    }
                }
                for (int c = 0; c < RCT_R2P / 16; ++c)
                    _mm256_store_si256((__m256i *)&t1[c * 16],
                                       rct_q_barrett16(_mm256_load_si256((const __m256i *)&t1[c * 16])));
                uint16_t t2[SNOVA_lr] = {0};
                for (int i1 = 0; i1 < SNOVA_r; ++i1)
                    for (int j1 = 0; j1 < SNOVA_l; ++j1) {
                        uint16_t acc = 0;
                        for (int k1 = 0; k1 < SNOVA_r; ++k1)
                            acc = (uint16_t)(acc + t1[i1 * SNOVA_r + k1] * Bmr[k1 * SNOVA_l + j1]);
                        t2[i1 * SNOVA_l + j1] = (uint16_t)(acc % SNOVA_q);
                    }
                for (int i1 = 0; i1 < SNOVA_r; ++i1)
                    for (int j1 = 0; j1 < SNOVA_l; ++j1) {
                        uint16_t acc = 0;
                        for (int k1 = 0; k1 < SNOVA_r; ++k1)
                            acc = (uint16_t)(acc + Amr[i1 * SNOVA_r + k1] * t2[k1 * SNOVA_l + j1]);
                        uint16_t *f = &Fvvacc[mi * SNOVA_lr + i1 * SNOVA_l + j1];
                        *f = (uint16_t)((*f + acc) % SNOVA_q);
                    }
            }
        for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i)
            Fvv_in_GF16Matrix[i] = (gf_t)Fvvacc[i];
        SNOVA_CLEAR_OBJ(Fvvacc);
    }
}
#endif

#if RCT_OQDF
static void rct_sign_whipbuild_oqdf(rct_sign_ctx *c) {
    const gf_t *signature_in_GF = c->signature_in_GF;
    _Static_assert((uint32_t)SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "RCT_OQ_DEVFLOW: u16 whip MAC overflow guard");
    memset(rct_oq_whip_w, 0, sizeof(rct_oq_whip_w));
    for (int ab = 0; ab < SNOVA_l; ++ab)
        for (int ni = 0; ni < SNOVA_v; ++ni)
            for (int i1 = 0; i1 < SNOVA_l; i1++)
                for (int k1 = 0; k1 < SNOVA_l; k1++) {
                    const uint16_t s = rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1];
                    const gf_t *sg = &signature_in_GF[ni * SNOVA_lr + k1 * SNOVA_r];
                    uint16_t *w = &rct_oq_whip_w[(i1 * SNOVA_v + ni) * RCT_Q_LRP + ab * SNOVA_r];
                    for (int j1 = 0; j1 < SNOVA_r; j1++)
                        w[j1] = (uint16_t)(w[j1] + (uint16_t)(s * (uint16_t)sg[j1]));
                }
    for (int i = 0; i < SNOVA_l * SNOVA_v * RCT_Q_LRP; i += 16)
        _mm256_store_si256((__m256i *)&rct_oq_whip_w[i],
                           rct_q_barrett16(_mm256_load_si256((const __m256i *)&rct_oq_whip_w[i])));
}
#endif

#if !RCT_OQDF && !(RCT_USE_SIMD && SNOVA_r <= 16) && RCT_OQWV
static void rct_sign_whipbuild_oqwv(rct_sign_ctx *c) {
    const gf_t *signature_in_GF = c->signature_in_GF;
    gf_t *whipped_sig = c->whipped_sig;
    _Static_assert((uint32_t)SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "RCT_OQ_WHIPVEC: u16 whip MAC overflow guard");
    {
        static uint16_t ws_acc[SNOVA_l * SNOVA_v * SNOVA_lr];
        memset(ws_acc, 0, sizeof(ws_acc));
        for (int ab = 0; ab < SNOVA_l; ++ab)
            for (int ni = 0; ni < SNOVA_v; ++ni)
                for (int i1 = 0; i1 < SNOVA_l; i1++)
                    for (int k1 = 0; k1 < SNOVA_l; k1++) {
                        const uint16_t s = rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1];
                        const gf_t *sg = &signature_in_GF[ni * SNOVA_lr + k1 * SNOVA_r];
                        uint16_t *w = &ws_acc[(ab * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r];
                        for (int j1 = 0; j1 < SNOVA_r; j1++)
                            w[j1] = (uint16_t)(w[j1] + (uint16_t)(s * (uint16_t)sg[j1]));
                    }
        const int rct_oqwv_n = SNOVA_l * SNOVA_v * SNOVA_lr;
#if RCT_Q_HAVE_MAGIC && defined(__AVX2__)
        int oi = 0;
        for (; oi + 16 <= rct_oqwv_n; oi += 16) {
            _Alignas(32) uint16_t rbuf[16];
            _mm256_store_si256((__m256i *)rbuf,
                               rct_q_barrett16(_mm256_loadu_si256((const __m256i *)&ws_acc[oi])));
            for (int t = 0; t < 16; ++t) whipped_sig[oi + t] = (gf_t)rbuf[t];
        }
        for (; oi < rct_oqwv_n; ++oi) whipped_sig[oi] = (gf_t)(ws_acc[oi] % SNOVA_q);
#else
        for (int oi = 0; oi < rct_oqwv_n; ++oi) whipped_sig[oi] = (gf_t)(ws_acc[oi] % SNOVA_q);
#endif
        SNOVA_CLEAR_OBJ(ws_acc);
    }
}
#endif

#if RCT_Q_SIMD && (SNOVA_r == SNOVA_l)
static void rct_sign_sumt_oddq_sq(rct_sign_ctx *c) {
    const gf_t *P11 = c->P11, *whipped_sig = c->whipped_sig;
    gf_t *sum_t1 = c->sum_t1;
    _Static_assert((uint32_t)SNOVA_v * SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "u16 accumulation overflow guard");
    {
        static uint16_t sum_t0u[SNOVA_m1 * SNOVA_l * SNOVA_v * SNOVA_lr];
        static uint16_t sum_t1u[SNOVA_m1 * SNOVA_l2 * SNOVA_r2];
        memset(sum_t0u, 0, sizeof(sum_t0u));
        memset(sum_t1u, 0, sizeof(sum_t1u));
        for (int mi = 0; mi < SNOVA_m1; ++mi)
            for (int nj = 0; nj < SNOVA_v; ++nj)
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        rct_q_matmul4_add(&sum_t0u[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr],
                                          &P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2],
                                          &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr]);
        static gf_t sum_t0b[SNOVA_m1 * SNOVA_l * SNOVA_v * SNOVA_lr + 16];
        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_l * SNOVA_lr; i1++)
            sum_t0b[i1] = (gf_t)(sum_t0u[i1] % SNOVA_q);
        for (int mi = 0; mi < SNOVA_m1; ++mi)
            for (int a1 = 0; a1 < SNOVA_l; ++a1)
                for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                    uint16_t *dst = &sum_t1u[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2];
                    for (int ni = 0; ni < SNOVA_v; ++ni)
                        rct_q_matmul4T_add(dst,
                                           &whipped_sig[(a1 * SNOVA_v + ni) * SNOVA_lr],
                                           &sum_t0b[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr]);
                }
        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l2 * SNOVA_r2; i1++)
            sum_t1[i1] = (gf_t)(sum_t1u[i1] % SNOVA_q);
        SNOVA_CLEAR_OBJ(sum_t0u);
        SNOVA_CLEAR_OBJ(sum_t1u);
        SNOVA_CLEAR_OBJ(sum_t0b);
    }
}
#endif

#if RCT_Q_SIMD && RCT_Q_HAVE_MAGIC && (SNOVA_r != SNOVA_l)
static void rct_sign_sumt_oddq_rect(rct_sign_ctx *c) {
#if !RCT_OQDF
    const gf_t *P11 = c->P11, *whipped_sig = c->whipped_sig;
    gf_t *sum_t1 = c->sum_t1;
#endif
    _Static_assert((uint32_t)SNOVA_v * SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "u16 accumulation overflow guard");
    {
        enum { RCT_LRP = RCT_Q_LRP, RCT_LR16 = RCT_Q_LR16 };
        _Alignas(32) static uint16_t sum_t0w[SNOVA_m1 * SNOVA_l * SNOVA_v * RCT_LRP];
        _Alignas(32) static uint16_t sum_t1w[SNOVA_m1 * SNOVA_l * SNOVA_r * RCT_LRP];
        memset(sum_t0w, 0, sizeof(sum_t0w));
        memset(sum_t1w, 0, sizeof(sum_t1w));
#if RCT_OQDF
        const uint16_t *whip_w = rct_oq_whip_w;
        const uint16_t *P11b = rct_oq_P11u;
#else
        const gf_t *P11b = P11;
        _Alignas(32) static uint16_t whip_w[SNOVA_l * SNOVA_v * RCT_LRP];
        memset(whip_w, 0, sizeof(whip_w));
        for (int ab = 0; ab < SNOVA_l; ++ab)
            for (int ni = 0; ni < SNOVA_v; ++ni)
                for (int i1 = 0; i1 < SNOVA_l; ++i1)
                    for (int j1 = 0; j1 < SNOVA_r; ++j1)
                        whip_w[(i1 * SNOVA_v + ni) * RCT_LRP + ab * SNOVA_r + j1] =
                            whipped_sig[(ab * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r + j1];
#endif
        for (int mi = 0; mi < SNOVA_m1; ++mi)
            for (int nj = 0; nj < SNOVA_v; ++nj)
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                        __m256i *s0 = (__m256i *)&sum_t0w[((mi * SNOVA_l + i1) * SNOVA_v + ni) * RCT_LRP];
                        for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                            __m256i qv = _mm256_set1_epi16((short)
                                P11b[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
                            const __m256i *wp = (const __m256i *)&whip_w[(k1 * SNOVA_v + nj) * RCT_LRP];
                            for (int b1 = 0; b1 < RCT_LR16; ++b1)
                                s0[b1] = _mm256_add_epi16(s0[b1], _mm256_mullo_epi16(qv, wp[b1]));
                        }
                    }
        for (int i = 0; i < SNOVA_m1 * SNOVA_l * SNOVA_v * RCT_LRP; i += 16)
            _mm256_store_si256((__m256i *)&sum_t0w[i],
                               rct_q_barrett16(_mm256_load_si256((const __m256i *)&sum_t0w[i])));
        for (int mi = 0; mi < SNOVA_m1; ++mi)
            for (int ni = 0; ni < SNOVA_v; ++ni)
                for (int a1 = 0; a1 < SNOVA_l; ++a1)
                    for (int i1 = 0; i1 < SNOVA_r; ++i1) {
                        __m256i *s1 = (__m256i *)&sum_t1w[((mi * SNOVA_l + a1) * SNOVA_r + i1) * RCT_LRP];
                        for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                            __m256i wv = _mm256_set1_epi16((short)
                                whip_w[(k1 * SNOVA_v + ni) * RCT_LRP + a1 * SNOVA_r + i1]);
                            const __m256i *s0 = (const __m256i *)&sum_t0w[((mi * SNOVA_l + k1) * SNOVA_v + ni) * RCT_LRP];
                            for (int b1 = 0; b1 < RCT_LR16; ++b1)
                                s1[b1] = _mm256_add_epi16(s1[b1], _mm256_mullo_epi16(wv, s0[b1]));
                        }
                    }
        for (int i = 0; i < SNOVA_m1 * SNOVA_l * SNOVA_r * RCT_LRP; i += 16)
            _mm256_store_si256((__m256i *)&sum_t1w[i],
                               rct_q_barrett16(_mm256_load_si256((const __m256i *)&sum_t1w[i])));
        for (int mi = 0; mi < SNOVA_m1; ++mi)
            for (int a1 = 0; a1 < SNOVA_l; ++a1)
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    for (int i1 = 0; i1 < SNOVA_r; ++i1)
                        for (int j1 = 0; j1 < SNOVA_r; ++j1)
#if RCT_OQDF
                            rct_oq_sum_t1u[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1] =
                                sum_t1w[((mi * SNOVA_l + a1) * SNOVA_r + i1) * RCT_LRP + b1 * SNOVA_r + j1];
#else
                            sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1] =
                                (gf_t)sum_t1w[((mi * SNOVA_l + a1) * SNOVA_r + i1) * RCT_LRP + b1 * SNOVA_r + j1];
#endif
#if !RCT_OQDF
        SNOVA_CLEAR_OBJ(whip_w);
#endif
        SNOVA_CLEAR_OBJ(sum_t0w);
        SNOVA_CLEAR_OBJ(sum_t1w);
    }
}
#endif

#if RCT_Q_SIMD && (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC
static void rct_sign_wF_gauss_oddq_rect(rct_sign_ctx *c) {
    const gf_t *F21 = c->F21, *F12 = c->F12, *whipped_sig = c->whipped_sig;
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm, *Q1 = c->Q1, *Q2 = c->Q2;
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
    (void)F21; (void)F12; (void)whipped_sig;
    (void)q1; (void)q2; (void)Am; (void)Bm; (void)Q1; (void)Q2;
        _Static_assert((uint32_t)SNOVA_v * SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                       "u16 accumulation overflow guard");
        _Static_assert((uint32_t)SNOVA_alpha * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                       "odd-q rect gausstmp cross-alpha uint16 accumulation may overflow");
        {
            enum { RCT_LRP = RCT_Q_LRP, RCT_LR16 = RCT_Q_LR16,
                   RCT_OLR16 = (SNOVA_o * SNOVA_lr / 16 + 1), RCT_OLR = RCT_OLR16 * 16 };
            _Alignas(32) static uint16_t wF21[SNOVA_m1 * SNOVA_l * SNOVA_o * RCT_LRP];
            _Alignas(32) static uint16_t wF12[SNOVA_m1 * SNOVA_l * SNOVA_o * RCT_LRP];
            memset(wF21, 0, sizeof(wF21));
            memset(wF12, 0, sizeof(wF12));
#if RCT_OQDF
            const uint16_t *whip_w = rct_oq_whip_w;
            const uint16_t *F21b = rct_oq_F21u;
            const uint16_t *F12b = rct_oq_F12u;
#else
            const gf_t *F21b = F21;
            const gf_t *F12b = F12;
            _Alignas(32) static uint16_t whip_w[SNOVA_l * SNOVA_v * RCT_LRP];
            memset(whip_w, 0, sizeof(whip_w));
            for (int ab = 0; ab < SNOVA_l; ++ab)
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int i1 = 0; i1 < SNOVA_l; ++i1)
                        for (int j1 = 0; j1 < SNOVA_r; ++j1)
                            whip_w[(i1 * SNOVA_v + ni) * RCT_LRP + ab * SNOVA_r + j1] =
                                whipped_sig[(ab * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r + j1];
#endif
            for (int mi = 0; mi < SNOVA_m1; ++mi)
                for (int idx = 0; idx < SNOVA_o; ++idx)
                    for (int nj = 0; nj < SNOVA_v; ++nj)
                        for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                            __m256i *s = (__m256i *)&wF21[((mi * SNOVA_l + i1) * SNOVA_o + idx) * RCT_LRP];
                            for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                __m256i qv = _mm256_set1_epi16((short)
                                    F21b[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
                                const __m256i *w = (const __m256i *)&whip_w[(k1 * SNOVA_v + nj) * RCT_LRP];
                                for (int b1 = 0; b1 < RCT_LR16; ++b1)
                                    s[b1] = _mm256_add_epi16(s[b1], _mm256_mullo_epi16(qv, w[b1]));
                            }
                        }
            for (int i = 0; i < SNOVA_m1 * SNOVA_l * SNOVA_o * RCT_LRP; i += 16)
                _mm256_store_si256((__m256i *)&wF21[i], rct_q_barrett16(_mm256_load_si256((const __m256i *)&wF21[i])));
            for (int mi = 0; mi < SNOVA_m1; ++mi)
                for (int nj = 0; nj < SNOVA_v; ++nj)
                    for (int idx = 0; idx < SNOVA_o; ++idx)
                        for (int i1 = 0; i1 < SNOVA_l; ++i1) {
                            __m256i *s = (__m256i *)&wF12[((mi * SNOVA_l + i1) * SNOVA_o + idx) * RCT_LRP];
                            for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                                __m256i qv = _mm256_set1_epi16((short)
                                    F12b[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1]);
                                const __m256i *w = (const __m256i *)&whip_w[(k1 * SNOVA_v + nj) * RCT_LRP];
                                for (int b1 = 0; b1 < RCT_LR16; ++b1)
                                    s[b1] = _mm256_add_epi16(s[b1], _mm256_mullo_epi16(qv, w[b1]));
                            }
                        }
            for (int i = 0; i < SNOVA_m1 * SNOVA_l * SNOVA_o * RCT_LRP; i += 16)
                _mm256_store_si256((__m256i *)&wF12[i], rct_q_barrett16(_mm256_load_si256((const __m256i *)&wF12[i])));

            _Alignas(32) static uint16_t gausstmp1[SNOVA_m1 * SNOVA_r * SNOVA_r * RCT_OLR];
            memset(gausstmp1, 0, sizeof(gausstmp1));
            for (int mi = 0; mi < SNOVA_o; ++mi)
                for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
                    int mi_prime = i_prime(mi, alpha);
                    const gf_t *q2r = &q2[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                    const gf_t *q1r = &q1[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                    const gf_t *Bmr = &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr];
                    const gf_t *Amr = &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2];
                    const gf_t *Q1r = &Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2];
                    (void)q1r;
                    _Alignas(32) uint16_t t0[SNOVA_o * SNOVA_lr] = {0};
                    _Alignas(32) uint16_t t1[SNOVA_o * SNOVA_l2] = {0};
                    _Alignas(32) uint16_t t2[RCT_OLR] = {0};
                    for (int idx = 0; idx < SNOVA_o; ++idx)
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            for (int i1 = 0; i1 < SNOVA_l; ++i1)
                                for (int j1 = 0; j1 < SNOVA_r; ++j1)
                                    t0[idx * SNOVA_lr + i1 * SNOVA_r + j1] +=
                                        wF21[((mi_prime * SNOVA_l + i1) * SNOVA_o + idx) * RCT_LRP + b1 * SNOVA_r + j1] * q2r[b1];
                    for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) t0[i] %= SNOVA_q;
                    for (int idx = 0; idx < SNOVA_o; ++idx)
                        for (int i1 = 0; i1 < SNOVA_l; ++i1)
                            for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                for (int k1 = 0; k1 < SNOVA_r; ++k1)
                                    t1[idx * SNOVA_l2 + i1 * SNOVA_l + j1] +=
                                        t0[idx * SNOVA_lr + i1 * SNOVA_r + k1] * Bmr[k1 * SNOVA_l + j1];
                    for (int i = 0; i < SNOVA_o * SNOVA_l2; ++i) t1[i] %= SNOVA_q;
                    for (int idx = 0; idx < SNOVA_o; ++idx)
                        for (int i1 = 0; i1 < SNOVA_l; ++i1)
                            for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                for (int k1 = 0; k1 < SNOVA_l; ++k1)
                                    t2[idx * SNOVA_l2 + i1 * SNOVA_l + j1] +=
                                        Q1r[i1 * SNOVA_l + k1] * t1[idx * SNOVA_l2 + k1 * SNOVA_l + j1];
                    for (int i = 0; i < SNOVA_o * SNOVA_l2; ++i) t2[i] %= SNOVA_q;
                    for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                        for (int tj2 = 0; tj2 < SNOVA_r; ++tj2) {
                            __m256i av = _mm256_set1_epi16((short)Amr[ti1 * SNOVA_r + tj2]);
                            __m256i *g = (__m256i *)&gausstmp1[((mi * SNOVA_r + ti1) * SNOVA_r + tj2) * RCT_OLR];
                            for (int c = 0; c < RCT_OLR / 16; ++c)
                                g[c] = _mm256_add_epi16(g[c],
                                    _mm256_mullo_epi16(av, _mm256_load_si256((const __m256i *)&t2[c * 16])));
                        }
                }
            for (int mi = 0; mi < SNOVA_o; ++mi)
                for (int idx = 0; idx < SNOVA_o; ++idx)
                    for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                        for (int tj1 = 0; tj1 < SNOVA_l; ++tj1)
                            for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                                for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                                    gauss[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx * SNOVA_lr + tj1 * SNOVA_r + tj2] =
                                        (gf_t)((gauss[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx * SNOVA_lr + tj1 * SNOVA_r + tj2] +
                                                gausstmp1[((mi * SNOVA_r + ti1) * SNOVA_r + tj2) * RCT_OLR +
                                                          idx * SNOVA_l2 + tj1 * SNOVA_l + ti2]) % SNOVA_q);

            _Alignas(32) static uint16_t gausstmp2[SNOVA_m1 * SNOVA_r * SNOVA_r * RCT_OLR];
            memset(gausstmp2, 0, sizeof(gausstmp2));
            for (int mi = 0; mi < SNOVA_o; ++mi)
                for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
                    int mi_prime = i_prime(mi, alpha);
                    const gf_t *q1r = &q1[(mi * SNOVA_alpha + alpha) * SNOVA_l];
                    const gf_t *Amr = &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2];
                    const gf_t *Bmr = &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr];
                    const gf_t *Q2r = &Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2];
                    _Alignas(32) uint16_t t0[SNOVA_o * SNOVA_lr] = {0};
                    _Alignas(32) uint16_t t1[SNOVA_o * SNOVA_lr] = {0};
                    _Alignas(32) uint16_t t2[RCT_OLR] = {0};
                    for (int idx = 0; idx < SNOVA_o; ++idx)
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            for (int i1 = 0; i1 < SNOVA_l; ++i1)
                                for (int j1 = 0; j1 < SNOVA_r; ++j1)
                                    t0[idx * SNOVA_lr + i1 * SNOVA_r + j1] +=
                                        wF12[((mi_prime * SNOVA_l + i1) * SNOVA_o + idx) * RCT_LRP + b1 * SNOVA_r + j1] * q1r[b1];
                    for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) t0[i] %= SNOVA_q;
                    for (int idx = 0; idx < SNOVA_o; ++idx)
                        for (int i1 = 0; i1 < SNOVA_r; ++i1)
                            for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                for (int k1 = 0; k1 < SNOVA_r; ++k1)
                                    t1[idx * SNOVA_lr + i1 * SNOVA_l + j1] +=
                                        Amr[i1 * SNOVA_r + k1] * t0[idx * SNOVA_lr + j1 * SNOVA_r + k1];
                    for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) t1[i] %= SNOVA_q;
                    for (int idx = 0; idx < SNOVA_o; ++idx)
                        for (int i1 = 0; i1 < SNOVA_r; ++i1)
                            for (int j1 = 0; j1 < SNOVA_l; ++j1)
                                for (int k1 = 0; k1 < SNOVA_l; ++k1)
                                    t2[idx * SNOVA_lr + i1 * SNOVA_l + j1] +=
                                        t1[idx * SNOVA_lr + i1 * SNOVA_l + k1] * Q2r[k1 * SNOVA_l + j1];
                    for (int i = 0; i < SNOVA_o * SNOVA_lr; ++i) t2[i] %= SNOVA_q;
                    for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                        for (int tj2 = 0; tj2 < SNOVA_r; ++tj2) {
                            __m256i bv = _mm256_set1_epi16((short)Bmr[tj2 * SNOVA_l + ti2]);
                            __m256i *g = (__m256i *)&gausstmp2[(mi * SNOVA_r * SNOVA_r + ti2 * SNOVA_r + tj2) * RCT_OLR];
                            for (int c = 0; c < RCT_OLR / 16; ++c)
                                g[c] = _mm256_add_epi16(g[c],
                                    _mm256_mullo_epi16(bv, _mm256_load_si256((const __m256i *)&t2[c * 16])));
                        }
                }
            for (int mi = 0; mi < SNOVA_o; ++mi)
                for (int idx = 0; idx < SNOVA_o; ++idx)
                    for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                        for (int tj1 = 0; tj1 < SNOVA_l; ++tj1)
                            for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                                for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                                    gauss[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx * SNOVA_lr + tj1 * SNOVA_r + tj2] =
                                        (gf_t)((gauss[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx * SNOVA_lr + tj1 * SNOVA_r + tj2] +
                                                gausstmp2[(mi * SNOVA_r * SNOVA_r + ti2 * SNOVA_r + tj2) * RCT_OLR +
                                                          idx * SNOVA_lr + ti1 * SNOVA_l + tj1]) % SNOVA_q);
#if !RCT_OQDF
            SNOVA_CLEAR_OBJ(whip_w);
#endif
            SNOVA_CLEAR_OBJ(wF21);
            SNOVA_CLEAR_OBJ(wF12);
            SNOVA_CLEAR_OBJ(gausstmp1);
            SNOVA_CLEAR_OBJ(gausstmp2);
        }
}
#endif

#if RCT_Q_SIMD
static void rct_skx_fold_F_oddq(gf_t *F21, gf_t *F12, const gf_t *T12, const gf_t *P11, const gf_t *P12, const gf_t *P21) {
    _Static_assert((uint32_t)SNOVA_v * 4u * (SNOVA_q - 1) * (SNOVA_q - 1)
                       + (SNOVA_q - 1) < 65536u,
                   "u16 accumulation overflow guard");
    {
        static uint16_t F21u[SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2];
        static uint16_t F12u[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2];

        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++) {
                __m256i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm256_setzero_si256();
                for (int j2 = 0; j2 < SNOVA_v; j2++) {
                    __m256i cperm[4];
                    rct_q_cperm(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2])), cperm);
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        acc[k1] = _mm256_add_epi16(acc[k1], rct_q_mm4_bc(
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&T12[(j2 * SNOVA_o + k1) * SNOVA_l2])), cperm));
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm256_storeu_si256((__m256i *)&F21u[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2], acc[k1]);
            }

        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++) {
                __m256i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm256_setzero_si256();
                for (int j2 = 0; j2 < SNOVA_v; j2++) {
                    __m256i bsh[4];
                    rct_q_bshuf(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&P11[((i1 * SNOVA_v + j1) * SNOVA_v + j2) * SNOVA_l2])), bsh);
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        acc[k1] = _mm256_add_epi16(acc[k1], rct_q_mm4_bs(bsh,
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&T12[(j2 * SNOVA_o + k1) * SNOVA_l2]))));
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm256_storeu_si256((__m256i *)&F12u[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2], acc[k1]);
            }

        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++)
            F12[i1] = (gf_t)((F12u[i1] + P12[i1]) % SNOVA_q);
        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2; i1++)
            F21[i1] = (gf_t)((F21u[i1] + P21[i1]) % SNOVA_q);
        SNOVA_CLEAR_OBJ(F21u);
        SNOVA_CLEAR_OBJ(F12u);
    }
}
#endif

#endif
