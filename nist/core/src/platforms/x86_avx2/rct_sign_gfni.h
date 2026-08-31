#ifndef RCT_SIGN_GFNI_H
#define RCT_SIGN_GFNI_H

#if SNOVA_Q == 16 && (RCT_USE_GFNI || RCT_HOT_QRP16) && !defined(RCT_GAUSS_SCALAR) && SNOVA_L == 4 \
    && !(RCT_HOT_QRP16 && !RCT_USE_GFNI && defined(RCT_T12_MULLO) && (RCT_T12_MULLO + 0))
static void rct_sign_apply_t12_gfni(rct_sign_ctx *c, const gf_t *solpad) {
    gf_t *signature_in_GF = c->signature_in_GF;
    const gf_t *T12 = c->T12;
    for (int index = 0; index < SNOVA_v; ++index) {
        gf_t *sigrow = &signature_in_GF[index * SNOVA_lr];
        for (int i1 = 0; i1 < SNOVA_l; ++i1) {
            __m128i acc = _mm_setzero_si128();
            const gf_t *tb = &T12[(index * SNOVA_o) * SNOVA_l2 + i1 * SNOVA_l];
            for (int mi = 0; mi < SNOVA_o; ++mi, tb += SNOVA_l2)
                for (int k1 = 0; k1 < SNOVA_l; ++k1)
                    acc = _mm_xor_si128(acc, RCT_GFMUL128(
                        _mm_set1_epi8((char)tb[k1]),
                        _mm_loadu_si128((const __m128i *)&solpad[mi * SNOVA_lr + k1 * SNOVA_r])));
            acc = rct_gfni_cleanup128(acc);
            uint8_t out16[16];
            _mm_storeu_si128((__m128i *)out16, acc);
            for (int j1 = 0; j1 < SNOVA_r; ++j1) sigrow[i1 * SNOVA_r + j1] ^= out16[j1];
        }
    }
}
#endif

#if SNOVA_Q == 16 && RCT_GFMUL_ANY && !defined(RCT_GAUSS_SCALAR)
static void rct_sign_backsub_gfni(rct_sign_ctx *c, gf_t *solpad, gf_t *solution) {
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
    memset(solpad, 0, (size_t)(RCT_SNB * 32 + 32) * sizeof(gf_t));
    for (int i = SNOVA_o * SNOVA_lr - 1; i >= 0; --i) {
        __m256i acc = _mm256_setzero_si256();
        for (int b = (i + 1) / 32; b < RCT_SNB; ++b)
            acc = _mm256_xor_si256(acc, RCT_GFMUL256(
                _mm256_loadu_si256((const __m256i *)&gauss[i][b * 32]),
                _mm256_load_si256((const __m256i *)&solpad[b * 32])));
#if RCT_USE_GFNI || RCT_HOT_QRP16
        acc = rct_gfni_cleanup256(acc);
#else
        acc = rct_sj_cleanup256(acc);
#endif
        __m128i x = _mm_xor_si128(_mm256_castsi256_si128(acc),
                                  _mm256_extracti128_si256(acc, 1));
        x = _mm_xor_si128(x, _mm_srli_si128(x, 8));
        x = _mm_xor_si128(x, _mm_srli_si128(x, 4));
        uint32_t w = (uint32_t)_mm_cvtsi128_si32(x);
        w ^= w >> 16; w ^= w >> 8;
        solpad[i] = (gf_t)(gauss[i][SNOVA_o * SNOVA_lr] ^ (w & 0x0fu));
    }
    for (int cc = 0; cc < SNOVA_o * SNOVA_lr; ++cc) solution[cc] = solpad[cc];
}
#endif

#if !(RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR))
static int rct_sign_gauss_q16(rct_sign_ctx *c) {
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
    int flag_redo;
    {
        enum { RCT_OLR = SNOVA_o * SNOVA_lr,
               RCT_W = SNOVA_o * SNOVA_lr + 1,
               RCT_WNB = (SNOVA_o * SNOVA_lr + 1 + 31) / 32 };
        uint32_t redo_acc = 0;
        for (int i = 0; i < RCT_OLR; ++i) {
#if (RCT_USE_SIMD || RCT_SIGN_JOG) && !defined(RCT_GAUSS_SCALAR)
            const int b0 = i / 32;
            __m256i rowi[RCT_WNB];
            for (int b = b0; b < RCT_WNB; ++b)
                rowi[b] = _mm256_loadu_si256((const __m256i *)&gauss[i][b * 32]);
            uint8_t ii = gauss[i][i];
            for (int j = i + 1; j < RCT_OLR; ++j) {
                uint32_t need = (1u - ct_gf_nz((uint32_t)ii)) & ct_gf_nz((uint32_t)gauss[j][i]);
                __m256i mv = _mm256_set1_epi8((char)(0u - need));
                for (int b = b0; b < RCT_WNB; ++b)
                    rowi[b] = _mm256_xor_si256(rowi[b],
                        _mm256_and_si256(_mm256_loadu_si256((const __m256i *)&gauss[j][b * 32]), mv));
                ii ^= (uint8_t)((0u - need) & (uint32_t)gauss[j][i]);
            }
            for (int b = b0; b < RCT_WNB; ++b)
                _mm256_storeu_si256((__m256i *)&gauss[i][b * 32], rowi[b]);
            redo_acc |= (1u - ct_gf_nz((uint32_t)ii));
            gf_t t_GF16 = gf_inv_sec(ii);
            rct_gauss_row_scale_sec(gauss[i], t_GF16, i, RCT_W);
            for (int j = i + 1; j < RCT_OLR; ++j)
                rct_gauss_row_axpy_sec(gauss[j], gauss[i], gauss[j][i], i, RCT_W);
#else
            for (int j = i + 1; j < RCT_OLR; ++j) {
                uint32_t need = (1u - ct_gf_nz((uint32_t)gauss[i][i])) & ct_gf_nz((uint32_t)gauss[j][i]);
                gf_t m = (gf_t)(0u - need);
                for (int k = i; k < RCT_W; ++k) gauss[i][k] = gf_add(gauss[i][k], (gf_t)(gauss[j][k] & m));
            }
            redo_acc |= (1u - ct_gf_nz((uint32_t)gauss[i][i]));
            gf_t t_GF16 = gf_inv_sec(gauss[i][i]);
            for (int k = i; k < RCT_W; ++k) gauss[i][k] = gf_mult_sec(gauss[i][k], t_GF16);
            for (int j = i + 1; j < RCT_OLR; ++j) {
                gf_t gji = gauss[j][i];
                for (int k = i; k < RCT_W; ++k)
                    gauss[j][k] = gf_sub(gauss[j][k], gf_mult_sec(gauss[i][k], gji));
            }
#endif
        }
        flag_redo = (int)redo_acc;
        SNOVA_CT_DECLASSIFY(&flag_redo, sizeof flag_redo);
    }
    return flag_redo;
}
#endif

#if !RCT_SIGN_JOG \
    && !(RCT_Q_SIMD && (SNOVA_r == SNOVA_l) && RCT_Q_HAVE_MAGIC && !defined(RCT_FVV_SCALAR)) \
    && !(RCT_Q_SIMD && (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC) \
    && !defined(RCT_CMS2_ONLY)
static void rct_sign_fvv_std(rct_sign_ctx *c) {
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm, *sum_t1 = c->sum_t1;
    gf_t *Fvv_in_GF16Matrix = c->Fvv;
#if RCT_USE_GFNI && (SNOVA_r == SNOVA_l) && (SNOVA_L == 4)
    const __m128i sD0 = _mm_setr_epi8(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12);
    const __m128i sD1 = _mm_setr_epi8(1,1,1,1,5,5,5,5,9,9,9,9,13,13,13,13);
    const __m128i sD2 = _mm_setr_epi8(2,2,2,2,6,6,6,6,10,10,10,10,14,14,14,14);
    const __m128i sD3 = _mm_setr_epi8(3,3,3,3,7,7,7,7,11,11,11,11,15,15,15,15);
    const __m128i sE0 = _mm_setr_epi8(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3);
    const __m128i sE1 = _mm_setr_epi8(4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7);
    const __m128i sE2 = _mm_setr_epi8(8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11);
    const __m128i sE3 = _mm_setr_epi8(12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15);
    for (int mi = 0; mi < SNOVA_o; ++mi) {
        __m128i Facc = _mm_setzero_si128();
        for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
            int mi_prime = i_prime(mi, alpha);
            __m128i t1 = _mm_setzero_si128();
            for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                __m128i t0 = _mm_setzero_si128();
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    t0 = _mm_xor_si128(t0, RCT_SV128(
                        RCT_BC128(q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]),
                        _mm_loadu_si128((const __m128i *)
                            &sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2])));
                t0 = rct_gfni_cleanup128(t0);
                t1 = _mm_xor_si128(t1, RCT_SV128(
                    RCT_BC128(q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]), t0));
            }
            t1 = rct_gfni_cleanup128(t1);
            __m128i Bc = _mm_loadu_si128((const __m128i *)&Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr]);
            __m128i acc = RCT_GFMUL128(_mm_shuffle_epi8(t1, sD0), _mm_shuffle_epi8(Bc, sE0));
            acc = _mm_xor_si128(acc, RCT_GFMUL128(_mm_shuffle_epi8(t1, sD1), _mm_shuffle_epi8(Bc, sE1)));
            acc = _mm_xor_si128(acc, RCT_GFMUL128(_mm_shuffle_epi8(t1, sD2), _mm_shuffle_epi8(Bc, sE2)));
            acc = _mm_xor_si128(acc, RCT_GFMUL128(_mm_shuffle_epi8(t1, sD3), _mm_shuffle_epi8(Bc, sE3)));
            __m128i t2 = rct_gfni_cleanup128(acc);
            __m128i Ac = _mm_loadu_si128((const __m128i *)&Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2]);
            Facc = _mm_xor_si128(Facc, RCT_GFMUL128(_mm_shuffle_epi8(Ac, sD0), _mm_shuffle_epi8(t2, sE0)));
            Facc = _mm_xor_si128(Facc, RCT_GFMUL128(_mm_shuffle_epi8(Ac, sD1), _mm_shuffle_epi8(t2, sE1)));
            Facc = _mm_xor_si128(Facc, RCT_GFMUL128(_mm_shuffle_epi8(Ac, sD2), _mm_shuffle_epi8(t2, sE2)));
            Facc = _mm_xor_si128(Facc, RCT_GFMUL128(_mm_shuffle_epi8(Ac, sD3), _mm_shuffle_epi8(t2, sE3)));
        }
        _mm_storeu_si128((__m128i *)&Fvv_in_GF16Matrix[mi * SNOVA_lr],
                         rct_gfni_cleanup128(Facc));
    }
#else
    for (int mi = 0; mi < SNOVA_o; ++mi)
        for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
            int mi_prime = i_prime(mi, alpha);
            gf_t gf16m_temp1[SNOVA_r2] = {0};
            gf_t gf16m_temp2[SNOVA_lr + 16] = {0};
#if RCT_HOT_SIMD && SNOVA_l == 4 && SNOVA_r2 <= 64
            {
                __m256i t1lo = _mm256_setzero_si256(), t1hi = _mm256_setzero_si256();
                for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                    __m256i slo = _mm256_setzero_si256(), shi = _mm256_setzero_si256();
                    for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                        const gf_t *base = &sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2];
                        __m256i qv = RCT_BC(q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
                        slo = _mm256_xor_si256(slo, RCT_SV(qv, _mm256_loadu_si256((const __m256i *)base)));
                        shi = _mm256_xor_si256(shi, RCT_SV(qv, _mm256_loadu_si256((const __m256i *)(base + 32))));
                    }
                    _Alignas(32) uint8_t sb[64];
                    _mm256_store_si256((__m256i *)sb, slo);
                    _mm256_store_si256((__m256i *)(sb + 32), shi);
                    for (int k = 0; k < SNOVA_r2; ++k) sb[k] = rct_gfni_cleanup(sb[k]);
                    __m256i qa = RCT_BC(q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]);
                    t1lo = _mm256_xor_si256(t1lo, RCT_SV(qa, _mm256_loadu_si256((const __m256i *)sb)));
                    t1hi = _mm256_xor_si256(t1hi, RCT_SV(qa, _mm256_loadu_si256((const __m256i *)(sb + 32))));
                }
                _Alignas(32) uint8_t t1b[64];
                _mm256_store_si256((__m256i *)t1b, t1lo);
                _mm256_store_si256((__m256i *)(t1b + 32), t1hi);
                for (int k = 0; k < SNOVA_r2; ++k) gf16m_temp1[k] = rct_gfni_cleanup(t1b[k]);
            }
#else
            for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                gf_t sumb[SNOVA_r2] = {0};
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    for (int i1 = 0; i1 < SNOVA_r; i1++)
                        for (int j1 = 0; j1 < SNOVA_r; j1++)
                            gf_set_add(&sumb[i1 * SNOVA_r + j1],
                                       gf_mult_sec(sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1],
                                               q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]));
                for (int i1 = 0; i1 < SNOVA_r; i1++)
                    for (int j1 = 0; j1 < SNOVA_r; j1++)
                        gf_set_add(&gf16m_temp1[i1 * SNOVA_r + j1],
                                   gf_mult_sec(sumb[i1 * SNOVA_r + j1], q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]));
            }
#endif
#if RCT_USE_SIMD && SNOVA_l == 4
            rct_matmul_l4rows_sec(gf16m_temp2, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_r, SNOVA_r);
            rct_matmul_l4rows_add(&Fvv_in_GF16Matrix[mi * SNOVA_lr], &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp2, SNOVA_r, SNOVA_r);
#else
            gf_mat_mul_add_lr_sec(gf16m_temp2, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_r, SNOVA_r, SNOVA_l);
            gf_mat_mul_add_lr_sec(&Fvv_in_GF16Matrix[mi * SNOVA_lr], &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp2,
                              SNOVA_r, SNOVA_r, SNOVA_l);
#endif
        }
#endif
}
#endif

#if !RCT_OQDF && RCT_USE_SIMD && SNOVA_r <= 16
static void rct_sign_whipbuild_q16(rct_sign_ctx *c) {
    const gf_t *signature_in_GF = c->signature_in_GF;
    gf_t *whipped_sig = c->whipped_sig;
    for (int ab = 0; ab < SNOVA_l; ++ab)
        for (int ni = 0; ni < SNOVA_v; ++ni)
            for (int i1 = 0; i1 < SNOVA_l; i1++) {
                __m128i acc = _mm_setzero_si128();
                for (int k1 = 0; k1 < SNOVA_l; k1++)
                    acc = _mm_xor_si128(acc, RCT_SV128(
                        RCT_BC128(rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1]),
                        _mm_loadu_si128((const __m128i *)&signature_in_GF[ni * SNOVA_lr + k1 * SNOVA_r])));
                acc = rct_gfni_cleanup128(acc);
                _Alignas(16) uint8_t tmpws[16];
                _mm_store_si128((__m128i *)tmpws, acc);
                memcpy(&whipped_sig[(ab * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r], tmpws, SNOVA_r);
            }
}
#endif

#if RCT_USE_SIMD
static void rct_sign_s1whip_simd(rct_sign_ctx *c) {
    const gf_t *whipped_sig = c->whipped_sig;
    const gf_t *P11 = c->P11;
    gf_t *sum_t1 = c->sum_t1;
#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
    uint16_t *rct_cm_whip = c->cm->cm_whip;
#if RCT_USE_GFNI
    uint8_t *rct_cm_whipb = c->cm->cm_whipb;
#endif
#endif
        {
            RCT_SCRATCH _Alignas(32) uint8_t whipped_sig2[SNOVA_l * SNOVA_v * SNOVA_lr32];
            memset(whipped_sig2, 0, sizeof(whipped_sig2));
#if RCT_USE_SIMD && SNOVA_r <= 16 && !RCT_OQDF
            for (int idx = 0; idx < SNOVA_v; ++idx)
                for (int i1 = 0; i1 < SNOVA_l; i1++)
                    for (int ab = 0; ab < SNOVA_l; ++ab)
                        memcpy(&whipped_sig2[idx * SNOVA_l * SNOVA_lr32 + i1 * SNOVA_lr32 + ab * SNOVA_r],
                               &whipped_sig[(ab * SNOVA_v + idx) * SNOVA_lr + i1 * SNOVA_r], SNOVA_r);
#else
            const gf_t *signature_in_GF = c->signature_in_GF;
            for (int idx = 0; idx < SNOVA_v; ++idx)
                for (int i1 = 0; i1 < SNOVA_l; i1++)
                    for (int ab = 0; ab < SNOVA_l; ++ab) {
                        __m256i acc = _mm256_setzero_si256();
                        for (int k1 = 0; k1 < SNOVA_l; k1++)
                            acc = _mm256_xor_si256(acc, RCT_SV(
                                RCT_BC(rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1]),
                                _mm256_loadu_si256((const __m256i *)&signature_in_GF[idx * SNOVA_lr + k1 * SNOVA_r])));
                        _Alignas(32) uint8_t tmpws[32];
                        _mm256_store_si256((__m256i *)tmpws, rct_gfni_cleanup256(acc));
                        for (int j1 = 0; j1 < SNOVA_r; j1++)
                            whipped_sig2[idx * SNOVA_l * SNOVA_lr32 + i1 * SNOVA_lr32 + ab * SNOVA_r + j1] = tmpws[j1];
                    }
#endif
#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
            for (int idx = 0; idx < SNOVA_v; ++idx)
                for (int i1 = 0; i1 < SNOVA_l; i1++)
                    for (int c = 0; c < RCT_CMW; c += 16)
                        _mm256_store_si256(
                            (__m256i *)&rct_cm_whip[(i1 * SNOVA_v + idx) * RCT_CMW + c],
                            gf16_expand_u16x16(_mm256_cvtepu8_epi16(_mm_load_si128((const __m128i *)
                                &whipped_sig2[idx * SNOVA_l * SNOVA_lr32 + i1 * SNOVA_lr32 + c]))));
#if RCT_USE_GFNI
            for (int idx = 0; idx < SNOVA_v; ++idx)
                for (int i1 = 0; i1 < SNOVA_l; i1++)
                    memcpy(&rct_cm_whipb[(i1 * SNOVA_v + idx) * RCT_CMW],
                           &whipped_sig2[idx * SNOVA_l * SNOVA_lr32 + i1 * SNOVA_lr32], RCT_CMW);
#endif
#endif

#ifdef RCT_RD_ACTIVE
            enum { RCT_RD_VG = SNOVA_v / 2 };
            RCT_SCRATCH _Alignas(32) uint8_t whipped_sig2d[RCT_RD_VG * SNOVA_l * 32];
            for (int g = 0; g < RCT_RD_VG; ++g)
                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                    __m128i a2a = _mm_load_si128((const __m128i *)
                        &whipped_sig2[(2 * g) * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32]);
                    __m128i a2b = _mm_load_si128((const __m128i *)
                        &whipped_sig2[(2 * g + 1) * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32]);
                    _mm256_store_si256((__m256i *)&whipped_sig2d[(g * SNOVA_l + k1) * 32],
                        _mm256_inserti128_si256(_mm256_castsi128_si256(a2a), a2b, 1));
                }
#endif

#if RCT_USE_GFNI && (SNOVA_r == SNOVA_l) && (SNOVA_L == 4)
            {
                const __m256i s1_Pcol_y[4] = {
                    _mm256_broadcastsi128_si256(_mm_setr_epi8(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12)),
                    _mm256_broadcastsi128_si256(_mm_setr_epi8(1,1,1,1,5,5,5,5,9,9,9,9,13,13,13,13)),
                    _mm256_broadcastsi128_si256(_mm_setr_epi8(2,2,2,2,6,6,6,6,10,10,10,10,14,14,14,14)),
                    _mm256_broadcastsi128_si256(_mm_setr_epi8(3,3,3,3,7,7,7,7,11,11,11,11,15,15,15,15)),
                };
                const __m256i s1_Prow_y[4] = {
                    _mm256_broadcastsi128_si256(_mm_setr_epi8(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3)),
                    _mm256_broadcastsi128_si256(_mm_setr_epi8(4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7)),
                    _mm256_broadcastsi128_si256(_mm_setr_epi8(8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11)),
                    _mm256_broadcastsi128_si256(_mm_setr_epi8(12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15)),
                };
                const __m128i s1_Prow[4] = {
                    _mm_setr_epi8(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3),
                    _mm_setr_epi8(4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7),
                    _mm_setr_epi8(8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11),
                    _mm_setr_epi8(12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15),
                };
                const __m128i s1_Pcolt[4] = {
                    _mm_setr_epi8(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3),
                    _mm_setr_epi8(4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7),
                    _mm_setr_epi8(8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11),
                    _mm_setr_epi8(12,12,12,12,13,13,13,13,14,14,14,14,15,15,15,15),
                };
                RCT_SCRATCH __m128i s1_wsig[SNOVA_l][SNOVA_v];
                RCT_SCRATCH __m256i s1_Brow2[SNOVA_l / 2][SNOVA_v][SNOVA_l];
                for (int nj = 0; nj < SNOVA_v; ++nj)
                    for (int g = 0; g < SNOVA_l / 2; ++g) {
                        __m128i c0 = _mm_loadu_si128((const __m128i *)&whipped_sig[((2 * g) * SNOVA_v + nj) * SNOVA_lr]);
                        __m128i c1 = _mm_loadu_si128((const __m128i *)&whipped_sig[((2 * g + 1) * SNOVA_v + nj) * SNOVA_lr]);
                        s1_wsig[2 * g][nj]     = c0;
                        s1_wsig[2 * g + 1][nj] = c1;
                        __m256i cc = _mm256_inserti128_si256(_mm256_castsi128_si256(c0), c1, 1);
                        for (int k1 = 0; k1 < SNOVA_l; ++k1)
                            s1_Brow2[g][nj][k1] = _mm256_shuffle_epi8(cc, s1_Prow_y[k1]);
                    }
                RCT_SCRATCH __m128i s1_sum_t0[SNOVA_l][SNOVA_v];
                for (int mi = 0; mi < SNOVA_m1; ++mi) {
                    for (int ni = 0; ni < SNOVA_v; ++ni) {
                        __m256i acc[SNOVA_l / 2];
                        for (int g = 0; g < SNOVA_l / 2; ++g) acc[g] = _mm256_setzero_si256();
                        for (int nj = 0; nj < SNOVA_v; ++nj) {
                            __m256i fy = _mm256_broadcastsi128_si256(_mm_loadu_si128(
                                (const __m128i *)&P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2]));
                            __m256i Acol[SNOVA_l];
                            for (int k1 = 0; k1 < SNOVA_l; ++k1)
                                Acol[k1] = _mm256_shuffle_epi8(fy, s1_Pcol_y[k1]);
                            for (int g = 0; g < SNOVA_l / 2; ++g)
                                for (int k1 = 0; k1 < SNOVA_l; ++k1)
                                    acc[g] = _mm256_xor_si256(acc[g],
                                        RCT_GFMUL256(Acol[k1], s1_Brow2[g][nj][k1]));
                        }
                        for (int g = 0; g < SNOVA_l / 2; ++g) {
                            __m256i cy = rct_gfni_cleanup256(acc[g]);
                            s1_sum_t0[2 * g][ni]     = _mm256_castsi256_si128(cy);
                            s1_sum_t0[2 * g + 1][ni] = _mm256_extracti128_si256(cy, 1);
                        }
                    }
                    for (int a1 = 0; a1 < SNOVA_l; ++a1)
                        for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                            __m128i acc = _mm_setzero_si128();
                            for (int ni = 0; ni < SNOVA_v; ++ni) {
                                __m128i w = s1_wsig[a1][ni], s = s1_sum_t0[b1][ni];
                                for (int k1 = 0; k1 < SNOVA_l; ++k1)
                                    acc = _mm_xor_si128(acc, RCT_GFMUL128(
                                        _mm_shuffle_epi8(w, s1_Pcolt[k1]),
                                        _mm_shuffle_epi8(s, s1_Prow[k1])));
                            }
                            _mm_storeu_si128((__m128i *)&sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2],
                                             rct_gfni_cleanup128(acc));
                        }
                }
                SNOVA_CLEAR_OBJ(s1_wsig);
                SNOVA_CLEAR_OBJ(s1_Brow2);
                SNOVA_CLEAR_OBJ(s1_sum_t0);
            }
#else
            RCT_SCRATCH _Alignas(32) uint8_t sum_t1p[SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32];
            memset(sum_t1p, 0, sizeof(sum_t1p));
            for (int mi = 0; mi < SNOVA_m1; ++mi) {
                RCT_SCRATCH _Alignas(32) uint8_t sum_t0[SNOVA_v * SNOVA_l * SNOVA_lr32];
#ifndef RCT_RD_ACTIVE
                memset(sum_t0, 0, sizeof(sum_t0));
#endif
#ifdef RCT_RD_ACTIVE
                for (int ni = 0; ni < SNOVA_v; ++ni) {
                    const gf_t *pbase = &P11[((size_t)(mi * SNOVA_v + ni) * SNOVA_v) * SNOVA_l2];
                    for (int i1 = 0; i1 < SNOVA_l; i1++) {
                        __m256i acc = _mm256_setzero_si256();
                        for (int g = 0; g < RCT_RD_VG; ++g) {
                            const gf_t *p0 = pbase + (size_t)(2 * g) * SNOVA_l2 + i1 * SNOVA_l;
                            const gf_t *p1 = p0 + SNOVA_l2;
                            const __m256i *wd = (const __m256i *)&whipped_sig2d[(g * SNOVA_l) * 32];
                            for (int k1 = 0; k1 < SNOVA_l; k1++) {
                                __m256i tp = _mm256_inserti128_si256(
                                    _mm256_castsi128_si256(rct_vtl128(p0[k1])), rct_vtl128(p1[k1]), 1);
                                acc = _mm256_xor_si256(acc, _mm256_shuffle_epi8(tp, wd[k1]));
                            }
                        }
#if (SNOVA_v & 1)
                        {
                            const gf_t *pt = pbase + (size_t)(SNOVA_v - 1) * SNOVA_l2 + i1 * SNOVA_l;
                            const __m256i *wp = (const __m256i *)
                                &whipped_sig2[(SNOVA_v - 1) * SNOVA_l * SNOVA_lr32];
                            for (int k1 = 0; k1 < SNOVA_l; k1++)
                                acc = _mm256_xor_si256(acc, RCT_SV(RCT_BC(pt[k1]), wp[k1]));
                        }
#endif
                        _mm256_store_si256((__m256i *)&sum_t0[(ni * SNOVA_l + i1) * SNOVA_lr32],
                            _mm256_zextsi128_si256(_mm_xor_si128(_mm256_castsi256_si128(acc),
                                          _mm256_extracti128_si256(acc, 1))));
                    }
                }
#else
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int nj = 0; nj < SNOVA_v; ++nj)
                        for (int i1 = 0; i1 < SNOVA_l; i1++)
                            for (int k1 = 0; k1 < SNOVA_l; k1++)
                                for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
                                    __m256i *s0 = (__m256i *)&sum_t0[(ni * SNOVA_l + i1) * SNOVA_lr32];
                                    __m256i q1v = RCT_BC(
                                        P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
                                    __m256i *wp = (__m256i *)&whipped_sig2[nj * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32];
                                    s0[b1] = _mm256_xor_si256(s0[b1], RCT_SV(q1v, wp[b1]));
                                }
#endif
                for (int i = 0; i < SNOVA_v * SNOVA_l * SNOVA_lr32; ++i) sum_t0[i] = rct_gfni_cleanup(sum_t0[i]);
#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
                RCT_SCRATCH _Alignas(32) uint16_t sum_t0u[SNOVA_v * SNOVA_l * RCT_CMW];
#if !((SNOVA_r == SNOVA_l) && defined(RCT_SQ_CM_TIGHT) && (RCT_SQ_CM_TIGHT + 0))
                for (int i = 0; i < SNOVA_v * SNOVA_l * SNOVA_lr32; i += 32) {
                    __m256i b = _mm256_load_si256((const __m256i *)&sum_t0[i]);
                    _mm256_store_si256((__m256i *)&sum_t0u[i],
                                       _mm256_cvtepu8_epi16(_mm256_castsi256_si128(b)));
                    _mm256_store_si256((__m256i *)&sum_t0u[i + 16],
                                       _mm256_cvtepu8_epi16(_mm256_extracti128_si256(b, 1)));
                }
#else
                for (int row = 0; row < SNOVA_v * SNOVA_l; ++row)
                    _mm256_store_si256((__m256i *)&sum_t0u[row * RCT_CMW],
                        _mm256_cvtepu8_epi16(_mm_load_si128((const __m128i *)&sum_t0[row * SNOVA_lr32])));
#endif
                for (int a1 = 0; a1 < SNOVA_l; ++a1)
                    for (int i1 = 0; i1 < SNOVA_r; i1++) {
                        __m256i acc[RCT_CMW16];
                        for (int c = 0; c < RCT_CMW16; ++c) acc[c] = _mm256_setzero_si256();
                        for (int ni = 0; ni < SNOVA_v; ++ni)
                            for (int k1 = 0; k1 < SNOVA_l; k1++) {
                                __m256i wv = _mm256_set1_epi16((short)
                                    rct_cm_whip[(k1 * SNOVA_v + ni) * RCT_CMW + a1 * SNOVA_r + i1]);
                                const __m256i *s0 = (const __m256i *)&sum_t0u[(ni * SNOVA_l + k1) * RCT_CMW];
                                for (int c = 0; c < RCT_CMW16; ++c)
                                    acc[c] = _mm256_xor_si256(acc[c], _mm256_mullo_epi16(wv, s0[c]));
                            }
                        uint8_t *dst = &sum_t1p[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32];
                        for (int c = 0; c < RCT_CMW16; ++c)
                            _mm_store_si128((__m128i *)(dst + c * 16),
                                gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc[c])));
                    }
#else
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int a1 = 0; a1 < SNOVA_l; ++a1)
                        for (int k1 = 0; k1 < SNOVA_l; k1++)
                            for (int i1 = 0; i1 < SNOVA_r; i1++) {
                                __m256i wp = RCT_BC_SEC(
                                    whipped_sig2[ni * SNOVA_l * SNOVA_lr32 + k1 * SNOVA_lr32 + a1 * SNOVA_r + i1]);
                                for (int b1 = 0; b1 < SNOVA_lr16; ++b1) {
                                    __m256i *s1 = (__m256i *)&sum_t1p[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32];
                                    __m256i *s0 = (__m256i *)&sum_t0[(ni * SNOVA_l + k1) * SNOVA_lr32];
                                    s1[b1] = _mm256_xor_si256(s1[b1], RCT_SV(wp, s0[b1]));
                                }
                            }
#endif
                SNOVA_CLEAR_OBJ(sum_t0);
#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
                SNOVA_CLEAR_OBJ(sum_t0u);
#endif
            }
#if !defined(RCT_CM_ACTIVE) && !defined(RCT_CML_ONLY)
            for (int i = 0; i < SNOVA_m1 * SNOVA_l * SNOVA_r * SNOVA_lr32; ++i) sum_t1p[i] = rct_gfni_cleanup(sum_t1p[i]);
#endif
            for (int mi = 0; mi < SNOVA_m1; ++mi)
                for (int a1 = 0; a1 < SNOVA_l; ++a1)
                    for (int i1 = 0; i1 < SNOVA_r; i1++)
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            for (int j1 = 0; j1 < SNOVA_r; j1++)
                                sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1] =
                                    sum_t1p[(mi * SNOVA_l + a1) * SNOVA_r * SNOVA_lr32 + i1 * SNOVA_lr32 + b1 * SNOVA_r + j1];
            SNOVA_CLEAR_OBJ(sum_t1p);
#endif
            SNOVA_CLEAR_OBJ(whipped_sig2);
#ifdef RCT_RD_ACTIVE
            SNOVA_CLEAR_OBJ(whipped_sig2d);
#endif
        }
}
#endif

#if !RCT_Q_SIMD && RCT_USE_GFNI && (SNOVA_l == 4)
static void rct_skx_fold_F_gfni(gf_t *F21, gf_t *F12, const gf_t *T12, const gf_t *P11) {
    {
        const __m128i sA0 = _mm_setr_epi8(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12);
        const __m128i sA1 = _mm_setr_epi8(1,1,1,1,5,5,5,5,9,9,9,9,13,13,13,13);
        const __m128i sA2 = _mm_setr_epi8(2,2,2,2,6,6,6,6,10,10,10,10,14,14,14,14);
        const __m128i sA3 = _mm_setr_epi8(3,3,3,3,7,7,7,7,11,11,11,11,15,15,15,15);
        const __m128i sB0 = _mm_setr_epi8(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3);
        const __m128i sB1 = _mm_setr_epi8(4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7);
        const __m128i sB2 = _mm_setr_epi8(8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11);
        const __m128i sB3 = _mm_setr_epi8(12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15);
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++) {
                __m128i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm_setzero_si128();
                for (int j2 = 0; j2 < SNOVA_v; j2++) {
                    __m128i bv = _mm_loadu_si128((const __m128i *)&P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2]);
                    __m128i bk0 = _mm_shuffle_epi8(bv, sB0), bk1 = _mm_shuffle_epi8(bv, sB1),
                            bk2 = _mm_shuffle_epi8(bv, sB2), bk3 = _mm_shuffle_epi8(bv, sB3);
                    for (int k1 = 0; k1 < SNOVA_o; k1++) {
                        __m128i av = _mm_loadu_si128((const __m128i *)&T12[(j2 * SNOVA_o + k1) * SNOVA_l2]);
                        __m128i p = RCT_GFMUL128(_mm_shuffle_epi8(av, sA0), bk0);
                        p = _mm_xor_si128(p, RCT_GFMUL128(_mm_shuffle_epi8(av, sA1), bk1));
                        p = _mm_xor_si128(p, RCT_GFMUL128(_mm_shuffle_epi8(av, sA2), bk2));
                        p = _mm_xor_si128(p, RCT_GFMUL128(_mm_shuffle_epi8(av, sA3), bk3));
                        acc[k1] = _mm_xor_si128(acc[k1], p);
                    }
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2],
                                     rct_gfni_cleanup128(acc[k1]));
            }
        rct_l4g_fold_bsec(F12, P11, T12, SNOVA_v, 0);
    }
}
#endif

#endif
