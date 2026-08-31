#ifndef RCT_SIGN_MULLO_H
#define RCT_SIGN_MULLO_H

#if SNOVA_Q == 16 && RCT_HOT_QRP16 && !RCT_USE_GFNI && defined(RCT_T12_MULLO) \
    && (RCT_T12_MULLO + 0) && !defined(RCT_GAUSS_SCALAR) && SNOVA_L == 4
static void rct_sign_apply_t12_mullo(rct_sign_ctx *c, const gf_t *solpad) {
    gf_t *signature_in_GF = c->signature_in_GF;
    const gf_t *T12 = c->T12;
    for (int index = 0; index < SNOVA_v; ++index) {
        gf_t *sigrow = &signature_in_GF[index * SNOVA_lr];
        for (int i1 = 0; i1 < SNOVA_l; ++i1) {
            __m256i acc = _mm256_setzero_si256();
            const gf_t *tb = &T12[(index * SNOVA_o) * SNOVA_l2 + i1 * SNOVA_l];
            for (int mi = 0; mi < SNOVA_o; ++mi, tb += SNOVA_l2)
                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                    uint16_t te = (uint16_t)((tb[k1] | ((uint16_t)tb[k1] << 3)
                        | ((uint16_t)tb[k1] << 6) | ((uint16_t)tb[k1] << 9)) & 0x1111);
                    __m256i cv = _mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&solpad[mi * SNOVA_lr + k1 * SNOVA_r]));
                    acc = _mm256_xor_si256(acc,
                        _mm256_mullo_epi16(_mm256_set1_epi16((short)te), cv));
                }
            __m128i out = gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc));
            uint8_t out16[16];
            _mm_storeu_si128((__m128i *)out16, out);
            for (int j1 = 0; j1 < SNOVA_r; ++j1) sigrow[i1 * SNOVA_r + j1] ^= out16[j1];
        }
    }
}
#endif

#if defined(RCT_CMS2_ONLY)
static void rct_sign_fvv_cms2(rct_sign_ctx *c) {
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm, *sum_t1 = c->sum_t1;
    gf_t *Fvv_in_GF16Matrix = c->Fvv;
    for (int mi = 0; mi < SNOVA_o; ++mi) {
        __m256i Facc = _mm256_setzero_si256();
        for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
            int mi_prime = i_prime(mi, alpha);
            __m128i t1 = _mm_setzero_si128();
            for (int a1 = 0; a1 < SNOVA_l; ++a1) {
                __m128i t0 = _mm_setzero_si128();
                for (int b1 = 0; b1 < SNOVA_l; ++b1)
                    t0 = _mm_xor_si128(t0, _mm_shuffle_epi8(
                        rct_vtl128(q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]),
                        _mm_loadu_si128((const __m128i *)
                            &sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2])));
                t1 = _mm_xor_si128(t1, _mm_shuffle_epi8(
                    rct_vtl128(q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1]), t0));
            }
            __m128i Bc = _mm_loadu_si128((const __m128i *)&Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr]);
            __m256i acc = _mm256_setzero_si256();
            for (int k1 = 0; k1 < SNOVA_l; ++k1)
                acc = _mm256_xor_si256(acc, _mm256_mullo_epi16(
                    _mm256_cvtepu8_epi16(_mm_shuffle_epi8(t1, rct_s2_pcol(k1))),
                    gf16_expand_u16x16(_mm256_cvtepu8_epi16(_mm_shuffle_epi8(Bc, rct_s2_prow(k1))))));
            __m128i t2 = gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc));
            __m128i Ac = _mm_loadu_si128((const __m128i *)&Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2]);
            for (int k1 = 0; k1 < SNOVA_l; ++k1)
                Facc = _mm256_xor_si256(Facc, _mm256_mullo_epi16(
                    _mm256_cvtepu8_epi16(_mm_shuffle_epi8(Ac, rct_s2_pcol(k1))),
                    gf16_expand_u16x16(_mm256_cvtepu8_epi16(_mm_shuffle_epi8(t2, rct_s2_prow(k1))))));
        }
        _mm_storeu_si128((__m128i *)&Fvv_in_GF16Matrix[mi * SNOVA_lr],
            gf16_pack_u16_to_bytes(gf16_compress_u16x16(Facc)));
    }
}
#endif

#if !RCT_Q_SIMD && !(RCT_USE_GFNI && (SNOVA_l == 4)) && RCT_USE_PSHUFB && (SNOVA_l == 4) && defined(RCT_MULLO) && (RCT_MULLO + 0)
static void rct_skx_fold_F_mullo(gf_t *F21, gf_t *F12, const gf_t *T12, const gf_t *P11) {
    {
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++) {
                __m256i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm256_setzero_si256();
                for (int j2 = 0; j2 < SNOVA_v; j2++) {
                    __m256i cperm[4];
                    rct_gf16_cperm_exp(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2])), cperm);
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        acc[k1] = _mm256_xor_si256(acc[k1], rct_gf16_mm4_bc(
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&T12[(j2 * SNOVA_o + k1) * SNOVA_l2])), cperm));
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2],
                                     gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc[k1])));
            }
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++) {
                __m256i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm256_setzero_si256();
                for (int j2 = 0; j2 < SNOVA_v; j2++) {
                    __m256i bsh[4];
                    rct_gf16_bshuf_exp(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&P11[((i1 * SNOVA_v + j1) * SNOVA_v + j2) * SNOVA_l2])), bsh);
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        acc[k1] = _mm256_xor_si256(acc[k1], rct_gf16_mm4_bs(bsh,
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&T12[(j2 * SNOVA_o + k1) * SNOVA_l2]))));
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
                                     gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc[k1])));
            }
    }
}
#endif

#if !RCT_Q_SIMD && !(RCT_USE_GFNI && (SNOVA_l == 4)) && RCT_HOT_QRP16 && (SNOVA_l == 4) && !(defined(RCT_MULLO) && (RCT_MULLO + 0))
static void rct_skx_fold_F_qrp16(gf_t *F21, gf_t *F12, const gf_t *T12, const gf_t *P11) {
    {
        const __m128i sA0 = _mm_setr_epi8(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12);
        const __m128i sA1 = _mm_setr_epi8(1,1,1,1,5,5,5,5,9,9,9,9,13,13,13,13);
        const __m128i sA2 = _mm_setr_epi8(2,2,2,2,6,6,6,6,10,10,10,10,14,14,14,14);
        const __m128i sA3 = _mm_setr_epi8(3,3,3,3,7,7,7,7,11,11,11,11,15,15,15,15);
        const __m128i sB0 = _mm_setr_epi8(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3);
        const __m128i sB1 = _mm_setr_epi8(4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7);
        const __m128i sB2 = _mm_setr_epi8(8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11);
        const __m128i sB3 = _mm_setr_epi8(12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15);
        const __m256i SA01 = _mm256_set_m128i(sA1, sA0), SA23 = _mm256_set_m128i(sA3, sA2);
        const __m256i SB01 = _mm256_set_m128i(sB1, sB0), SB23 = _mm256_set_m128i(sB3, sB2);
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++) {
                __m128i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm_setzero_si128();
                for (int j2 = 0; j2 < SNOVA_v; j2++) {
                    __m256i bv = _mm256_broadcastsi128_si256(_mm_loadu_si128(
                        (const __m128i *)&P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2]));
                    __m256i B01 = _mm256_shuffle_epi8(bv, SB01), B23 = _mm256_shuffle_epi8(bv, SB23);
                    for (int k1 = 0; k1 < SNOVA_o; k1++) {
                        __m256i av = _mm256_broadcastsi128_si256(_mm_loadu_si128(
                            (const __m128i *)&T12[(j2 * SNOVA_o + k1) * SNOVA_l2]));
                        __m256i P = _mm256_xor_si256(
                            gf16_qrp16_256_byte_mul(_mm256_shuffle_epi8(av, SA01), B01),
                            gf16_qrp16_256_byte_mul(_mm256_shuffle_epi8(av, SA23), B23));
                        acc[k1] = _mm_xor_si128(acc[k1], _mm_xor_si128(
                            _mm256_castsi256_si128(P), _mm256_extracti128_si256(P, 1)));
                    }
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2],
                                     rct_gfni_cleanup128(acc[k1]));
            }
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++) {
                __m128i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm_setzero_si128();
                for (int j2 = 0; j2 < SNOVA_v; j2++) {
                    __m256i av = _mm256_broadcastsi128_si256(_mm_loadu_si128(
                        (const __m128i *)&P11[((i1 * SNOVA_v + j1) * SNOVA_v + j2) * SNOVA_l2]));
                    __m256i A01 = _mm256_shuffle_epi8(av, SA01), A23 = _mm256_shuffle_epi8(av, SA23);
                    for (int k1 = 0; k1 < SNOVA_o; k1++) {
                        __m256i bv = _mm256_broadcastsi128_si256(_mm_loadu_si128(
                            (const __m128i *)&T12[(j2 * SNOVA_o + k1) * SNOVA_l2]));
                        __m256i P = _mm256_xor_si256(
                            gf16_qrp16_256_byte_mul(A01, _mm256_shuffle_epi8(bv, SB01)),
                            gf16_qrp16_256_byte_mul(A23, _mm256_shuffle_epi8(bv, SB23)));
                        acc[k1] = _mm_xor_si128(acc[k1], _mm_xor_si128(
                            _mm256_castsi256_si128(P), _mm256_extracti128_si256(P, 1)));
                    }
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
                                     rct_gfni_cleanup128(acc[k1]));
            }
    }
}
#endif

#endif
