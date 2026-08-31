#ifndef RCT_SIGN_SCALAR_H
#define RCT_SIGN_SCALAR_H

#if !(RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC) \
    && !(SNOVA_Q == 16 && (RCT_USE_GFNI || RCT_HOT_QRP16) && !defined(RCT_GAUSS_SCALAR) && SNOVA_L == 4) \
    && !(RCT_SIGN_JOG && SNOVA_Q == 16 && RCT_GFMUL_ANY && !defined(RCT_GAUSS_SCALAR))
static void rct_sign_apply_t12_scalar(rct_sign_ctx *c, const gf_t *solution) {
    gf_t *signature_in_GF = c->signature_in_GF;
    const gf_t *T12 = c->T12;
    for (int index = 0; index < SNOVA_v; ++index)
        for (int mi = 0; mi < SNOVA_o; ++mi)
            gf_mat_mul_add_lr_sec(&signature_in_GF[index * SNOVA_lr], &T12[(index * SNOVA_o + mi) * SNOVA_l2],
                              &solution[mi * SNOVA_lr], SNOVA_l, SNOVA_l, SNOVA_r);
}
#endif

#if !(RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC) \
    && !(SNOVA_Q == 16 && RCT_GFMUL_ANY && !defined(RCT_GAUSS_SCALAR))
static void rct_sign_backsub_scalar(rct_sign_ctx *c, gf_t *solution) {
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
    memset(solution, 0, (size_t)(SNOVA_o * SNOVA_lr) * sizeof(gf_t));
    for (int i = SNOVA_o * SNOVA_lr - 1; i >= 0; --i) {
        gf_t sum = 0;
        for (int k = i + 1; k < SNOVA_o * SNOVA_lr; ++k)
            gf_set_add(&sum, gf_mult_sec(gauss[i][k], solution[k]));
        solution[i] = gf_sub(gauss[i][SNOVA_o * SNOVA_lr], sum);
    }
}
#endif

#if !RCT_SIGN_JOG && !defined(RCT_CM_ACTIVE) && !defined(RCT_CMS3_ONLY) \
    && !(RCT_Q_SIMD && (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC) \
    && !(RCT_Q_SIMD && (SNOVA_r == SNOVA_l))
static void rct_sign_wf_std(rct_sign_ctx *c) {
    const gf_t *F21 = c->F21, *F12 = c->F12, *whipped_sig = c->whipped_sig;
    gf_t *whipped_F21 = c->whipped_F21, *whipped_F12 = c->whipped_F12;
    for (int mi = 0; mi < SNOVA_m1; mi++) {
#if RCT_HOT_SIMD && SNOVA_l == 4
        for (int idx = 0; idx < SNOVA_o; idx++)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int i1 = 0; i1 < SNOVA_l; i1++) {
                    __m128i acc = _mm_setzero_si128();
                    for (int nj = 0; nj < SNOVA_v; ++nj)
                        for (int k1 = 0; k1 < SNOVA_l; k1++) {
                            __m128i s = RCT_BC128_SEC(
                                F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1]);
                            __m128i w = _mm_loadu_si128(
                                (const __m128i *)&whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr + k1 * SNOVA_r]);
                            acc = _mm_xor_si128(acc, RCT_SV128(s, w));
                        }
                    acc = rct_gfni_cleanup128(acc);
                    _Alignas(16) uint8_t tmp[16];
                    _mm_store_si128((__m128i *)tmp, acc);
                    memcpy(&whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r],
                           tmp, SNOVA_r);
                }
#else
        for (int idx = 0; idx < SNOVA_o; idx++)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int nj = 0; nj < SNOVA_v; ++nj)
                    gf_mat_mul_add_lr_sec(&whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr],
                                      &F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2],
                                      &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr], SNOVA_l, SNOVA_l, SNOVA_r);
#endif

#if RCT_HOT_SIMD && SNOVA_l == 4
        for (int idx = 0; idx < SNOVA_o; idx++)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int i1 = 0; i1 < SNOVA_l; i1++) {
                    __m128i acc = _mm_setzero_si128();
                    for (int nj = 0; nj < SNOVA_v; ++nj)
                        for (int k1 = 0; k1 < SNOVA_l; k1++) {
                            __m128i s = RCT_BC128_SEC(
                                F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1]);
                            __m128i w = _mm_loadu_si128(
                                (const __m128i *)&whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr + k1 * SNOVA_r]);
                            acc = _mm_xor_si128(acc, RCT_SV128(s, w));
                        }
                    acc = rct_gfni_cleanup128(acc);
                    _Alignas(16) uint8_t tmp[16];
                    _mm_store_si128((__m128i *)tmp, acc);
                    memcpy(&whipped_F12[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r],
                           tmp, SNOVA_r);
                }
#else
        for (int idx = 0; idx < SNOVA_o; idx++)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int nj = 0; nj < SNOVA_v; ++nj)
                    for (int i1 = 0; i1 < SNOVA_l; i1++)
                        for (int j1 = 0; j1 < SNOVA_r; j1++)
                            for (int k1 = 0; k1 < SNOVA_l; k1++)
                                gf_set_add(&whipped_F12[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r + j1],
                                           gf_mult_sec(F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1],
                                                   whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr + k1 * SNOVA_r + j1]));
#endif
    }
}
#endif

#if !defined(RCT_CM_ACTIVE) && !defined(RCT_CMS3_ONLY) \
    && !(RCT_Q_SIMD && (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC) \
    && !(RCT_Q_SIMD && (SNOVA_r == SNOVA_l))
#ifndef RCT_SIGN_P2_RSWAP
#define RCT_SIGN_P2_RSWAP 1
#endif
#if RCT_SIGN_P2_RSWAP && (RCT_SJ_A4 || RCT_SJ_M4)
#define RCT_SJ_RSWAP 1
#else
#define RCT_SJ_RSWAP 0
#endif

#if RCT_SJ_RSWAP
#if RCT_SJ_M4
#define RCT_RS_BC_T rct_m4_bc_t
#if SNOVA_r == 8
#define RCT_RS_OVW 48
#else
#define RCT_RS_OVW 32
#endif
#define RCT_RS_G16W ((SNOVA_o * SNOVA_lr + RCT_RS_OVW + 15) / 16 * 16)
static inline rct_m4_bc_t rct_rs_m4_bc16(const uint16_t *ev5) {
    __m256i evb = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)ev5));
    rct_m4_bc_t b;
    b.b0 = _mm256_shuffle_epi8(evb, RCT_M4_PAT0);
    b.b1 = _mm256_shuffle_epi8(evb, RCT_M4_PAT1);
#if SNOVA_r == 8
    b.b2 = _mm256_set1_epi16((short)ev5[4]);
#endif
    return b;
}
static inline void rct_rs_axpy16(uint16_t *g16, rct_m4_bc_t b, const gf_t *tile) {
    rct_m4_acc_t A = rct_m4_zero();
    rct_m4_mac(&A, b, tile);
    _mm256_storeu_si256((__m256i *)g16,
        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)g16), A.a0));
    _mm256_storeu_si256((__m256i *)&g16[16],
        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)&g16[16]), A.a1));
#if SNOVA_r == 8
    _mm256_storeu_si256((__m256i *)&g16[32],
        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)&g16[32]), A.a2));
#endif
}
#else
#if SNOVA_r == 8
#define RCT_RS_BC_T rct_a4_bc_t
#define RCT_RS_BC_FN(sc) rct_a4_bc(sc)
#define RCT_RS_OVW 48
static inline void rct_rs_axpy(gf_t *g, rct_a4_bc_t b, const gf_t *tile) {
    rct_a4_acc_t A = rct_a4_zero();
    rct_a4_mac(&A, b, tile);
    _mm256_storeu_si256((__m256i *)g,
        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)g), A.a));
    _mm_storeu_si128((__m128i *)&g[32],
        _mm_xor_si128(_mm_loadu_si128((const __m128i *)&g[32]), A.t));
}
#else
#define RCT_RS_BC_T __m256i
#define RCT_RS_BC_FN(sc) rct_rs_bc(sc)
#define RCT_RS_OVW 32
#define RCT_RS_BL(p) (char)((p) < SNOVA_lr ? (p) / SNOVA_r : -1)
static inline __m256i rct_rs_bc(const gf_t *sc) {
    const __m256i P = _mm256_setr_epi8(
        RCT_RS_BL(0), RCT_RS_BL(1), RCT_RS_BL(2), RCT_RS_BL(3),
        RCT_RS_BL(4), RCT_RS_BL(5), RCT_RS_BL(6), RCT_RS_BL(7),
        RCT_RS_BL(8), RCT_RS_BL(9), RCT_RS_BL(10), RCT_RS_BL(11),
        RCT_RS_BL(12), RCT_RS_BL(13), RCT_RS_BL(14), RCT_RS_BL(15),
        RCT_RS_BL(16), RCT_RS_BL(17), RCT_RS_BL(18), RCT_RS_BL(19),
        RCT_RS_BL(20), RCT_RS_BL(21), RCT_RS_BL(22), RCT_RS_BL(23),
        RCT_RS_BL(24), RCT_RS_BL(25), RCT_RS_BL(26), RCT_RS_BL(27),
        RCT_RS_BL(28), RCT_RS_BL(29), RCT_RS_BL(30), RCT_RS_BL(31));
    return _mm256_shuffle_epi8(_mm256_broadcastsi128_si256(
               _mm_loadu_si128((const __m128i *)sc)), P);
}
static inline void rct_rs_axpy(gf_t *g, __m256i bc, const gf_t *tile) {
    __m256i prod = _mm256_gf2p8mul_epi8(bc, _mm256_loadu_si256((const __m256i *)tile));
    _mm256_storeu_si256((__m256i *)g,
        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)g), prod));
}
#endif
#endif
_Static_assert((SNOVA_o - 1) * SNOVA_lr + RCT_RS_OVW <= SNOVA_o * SNOVA_lr + 1 + 64,
    "rswap: overlapped storeu must stay within gauss row (OLR+1+guard)");
#define RCT_RS_PL(C, p) (char)((16 * (C) + (p)) < SNOVA_lr ? ((16 * (C) + (p)) % SNOVA_r) : -1)
#define RCT_RS_PAT(C) _mm_setr_epi8( \
    RCT_RS_PL(C, 0), RCT_RS_PL(C, 1), RCT_RS_PL(C, 2), RCT_RS_PL(C, 3), \
    RCT_RS_PL(C, 4), RCT_RS_PL(C, 5), RCT_RS_PL(C, 6), RCT_RS_PL(C, 7), \
    RCT_RS_PL(C, 8), RCT_RS_PL(C, 9), RCT_RS_PL(C, 10), RCT_RS_PL(C, 11), \
    RCT_RS_PL(C, 12), RCT_RS_PL(C, 13), RCT_RS_PL(C, 14), RCT_RS_PL(C, 15))
#define RCT_RS_NCH ((SNOVA_lr + 16 + 15) / 16)
static inline void rct_rs_tile_build(uint8_t *tile, __m128i src) {
    _mm_store_si128((__m128i *)tile, _mm_shuffle_epi8(src, RCT_RS_PAT(0)));
    _mm_store_si128((__m128i *)&tile[16], _mm_shuffle_epi8(src, RCT_RS_PAT(1)));
    _mm_store_si128((__m128i *)&tile[32], _mm_shuffle_epi8(src, RCT_RS_PAT(2)));
#if RCT_RS_NCH > 3
    _mm_store_si128((__m128i *)&tile[48], _mm_shuffle_epi8(src, RCT_RS_PAT(3)));
#endif
}
#endif
static void rct_sign_gauss_scatter_std(rct_sign_ctx *c) {
    const gf_t *whipped_F21 = c->whipped_F21, *whipped_F12 = c->whipped_F12;
    const gf_t *q1 = c->q1, *q2 = c->q2, *Am = c->Am, *Bm = c->Bm, *Q1 = c->Q1, *Q2 = c->Q2;
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64] = c->gauss;
#if RCT_SJ_RSWAP
    _Alignas(16) gf_t rs_xT[SNOVA_l2 + 16];
    memset(rs_xT, 0, sizeof(rs_xT));
#if RCT_SJ_M4
    RCT_SCRATCH _Alignas(32) uint16_t rs_g16[SNOVA_o * SNOVA_lr][RCT_RS_G16W];
    memset(rs_g16, 0, sizeof(rs_g16));
    _Alignas(32) uint16_t rs_xe16[(SNOVA_lr + 16 + 15) / 16 * 16];
    memset(rs_xe16, 0, sizeof(rs_xe16));
#define RCT_RS_AXPY_ROW(row, ix, bcv, tile) rct_rs_axpy16(&rs_g16[row][(ix) * SNOVA_lr], bcv, tile)
#else
#define RCT_RS_AXPY_ROW(row, ix, bcv, tile) rct_rs_axpy(&gauss[row][(ix) * SNOVA_lr], bcv, tile)
#endif

    for (int mi = 0; mi < SNOVA_o; mi++)
        for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
            int mi_prime = i_prime(mi, alpha);
            _Alignas(16) uint8_t amtile[SNOVA_r][RCT_RS_NCH * 16];
            for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                rct_rs_tile_build(amtile[ti1], _mm_loadu_si128(
                    (const __m128i *)&Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r]));
            for (int idx = 0; idx < SNOVA_o; idx++) {
                gf_t gf16m_temp0[SNOVA_l2 + 16] = {0};
                gf_t gf16m_temp1[SNOVA_lr + 16] = {0};
                gf_t X_tmp[SNOVA_l2 + 16] = {0};
#ifdef RCT_PROFILE
                uint64_t _t0 = __rdtsc(), _t1, _t2;
#endif
                for (int cc = 0; cc < SNOVA_lr; cc += 16) {
                    __m128i acc = _mm_setzero_si128();
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        acc = _mm_xor_si128(acc, rct_sj_sv128(
                            rct_sj_bc128_pub(q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]),
                            _mm_loadu_si128((const __m128i *)&whipped_F21[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + cc])));
                    _Alignas(16) uint8_t pb[16]; _mm_store_si128((__m128i *)pb, RCT_SJ_CLEAN(acc));
                    int nn = SNOVA_lr - cc; if (nn > 16) nn = 16;
                    for (int j = 0; j < nn; ++j) gf16m_temp1[cc + j] = pb[j];
                }
#if RCT_SJ_M4
                rct_m4_mm_add(gf16m_temp0, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_l, SNOVA_r, SNOVA_l);
#else
                rct_sj_mm_add(gf16m_temp0, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_l, SNOVA_r, SNOVA_l);
#endif
                rct_sj_mm_add_pub(X_tmp, &Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2], gf16m_temp0, SNOVA_l, SNOVA_l, SNOVA_l);
#ifdef RCT_PROFILE
                _t1 = __rdtsc(); RCT_PACC6B(0, _t0, _t1);
#endif
                for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                    for (int tj1 = 0; tj1 < SNOVA_l; ++tj1)
                        rs_xT[ti2 * SNOVA_l + tj1] = X_tmp[tj1 * SNOVA_l + ti2];
#if RCT_SJ_M4
                rct_m4_expand_buf(rs_xe16, rs_xT, SNOVA_l2);
#endif
                for (int ti2 = 0; ti2 < SNOVA_l; ++ti2) {
#if RCT_SJ_M4
                    RCT_RS_BC_T bcv = rct_rs_m4_bc16(&rs_xe16[ti2 * SNOVA_l]);
#else
                    RCT_RS_BC_T bcv = RCT_RS_BC_FN(&rs_xT[ti2 * SNOVA_l]);
#endif
                    for (int ti1 = 0; ti1 < SNOVA_r; ++ti1)
                        RCT_RS_AXPY_ROW(mi * SNOVA_lr + ti1 * SNOVA_l + ti2, idx, bcv, amtile[ti1]);
                }
#ifdef RCT_PROFILE
                _t2 = __rdtsc(); RCT_PACC6B(1, _t1, _t2);
#endif
            }
        }

    for (int mi = 0; mi < SNOVA_o; mi++)
        for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
            int mi_prime = i_prime(mi, alpha);
            _Alignas(16) uint8_t bmtile[SNOVA_l][RCT_RS_NCH * 16];
            {
                _Alignas(16) uint8_t bmT[SNOVA_l * SNOVA_r + 16];
                memset(&bmT[SNOVA_l * SNOVA_r], 0, 16);
                for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                    for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                        bmT[ti2 * SNOVA_r + tj2] = Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2];
                for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                    rct_rs_tile_build(bmtile[ti2], _mm_loadu_si128((const __m128i *)&bmT[ti2 * SNOVA_r]));
            }
            for (int idx = 0; idx < SNOVA_o; idx++) {
                gf_t gf16m_temp0[SNOVA_lr + 16] = {0};
                gf_t gf16m_temp1[SNOVA_lr + 16] = {0};
                gf_t X_tmp[SNOVA_lr + 16] = {0};
#ifdef RCT_PROFILE
                uint64_t _t0 = __rdtsc(), _t1, _t2;
#endif
                {
                    _Alignas(16) uint8_t fbuf[SNOVA_lr + 16];
                    for (int cc = 0; cc < SNOVA_lr; cc += 16) {
                        __m128i acc = _mm_setzero_si128();
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            acc = _mm_xor_si128(acc, rct_sj_sv128(
                                rct_sj_bc128_pub(q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]),
                                _mm_loadu_si128((const __m128i *)&whipped_F12[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + cc])));
                        _mm_store_si128((__m128i *)&fbuf[cc], RCT_SJ_CLEAN(acc));
                    }
                    for (int i1 = 0; i1 < SNOVA_l; ++i1)
                        for (int j1 = 0; j1 < SNOVA_r; ++j1)
                            gf16m_temp1[j1 * SNOVA_l + i1] = fbuf[i1 * SNOVA_r + j1];
                }
                rct_sj_mm_add_pub(gf16m_temp0, &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp1, SNOVA_r, SNOVA_r, SNOVA_l);
#if RCT_SJ_M4
                rct_m4_mm_add(X_tmp, gf16m_temp0, &Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2], SNOVA_r, SNOVA_l, SNOVA_l);
#else
                rct_sj_mm_add(X_tmp, gf16m_temp0, &Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2], SNOVA_r, SNOVA_l, SNOVA_l);
#endif
#ifdef RCT_PROFILE
                _t1 = __rdtsc(); RCT_PACC6B(0, _t0, _t1);
#endif
#if RCT_SJ_M4
                rct_m4_expand_buf(rs_xe16, X_tmp, SNOVA_lr);
#endif
                for (int ti1 = 0; ti1 < SNOVA_r; ++ti1) {
#if RCT_SJ_M4
                    RCT_RS_BC_T bcv = rct_rs_m4_bc16(&rs_xe16[ti1 * SNOVA_l]);
#else
                    RCT_RS_BC_T bcv = RCT_RS_BC_FN(&X_tmp[ti1 * SNOVA_l]);
#endif
                    for (int ti2 = 0; ti2 < SNOVA_l; ++ti2)
                        RCT_RS_AXPY_ROW(mi * SNOVA_lr + ti1 * SNOVA_l + ti2, idx, bcv, bmtile[ti2]);
                }
#ifdef RCT_PROFILE
                _t2 = __rdtsc(); RCT_PACC6B(1, _t1, _t2);
#endif
            }
        }

#ifdef RCT_PROFILE
    uint64_t _tc0 = __rdtsc(), _tc1;
#endif
#if !RCT_SJ_M4
    for (int rr = 0; rr < SNOVA_o * SNOVA_lr; ++rr)
        for (int cc = 0; cc < SNOVA_o * SNOVA_lr; cc += 32)
            _mm256_storeu_si256((__m256i *)&gauss[rr][cc],
                rct_sj_cleanup256(_mm256_loadu_si256((const __m256i *)&gauss[rr][cc])));
#else
    for (int rr = 0; rr < SNOVA_o * SNOVA_lr; ++rr)
        for (int cc = 0; cc < SNOVA_o * SNOVA_lr; cc += 16) {
            __m128i pb = cl_gf16_pack_u16_to_bytes(cl_gf16_compress_u16x16(
                _mm256_loadu_si256((const __m256i *)&rs_g16[rr][cc])));
            _mm_storeu_si128((__m128i *)&gauss[rr][cc],
                _mm_xor_si128(_mm_loadu_si128((const __m128i *)&gauss[rr][cc]), pb));
        }
#endif
#ifdef RCT_PROFILE
    _tc1 = __rdtsc(); RCT_PACC6B(2, _tc0, _tc1);
#endif
    SNOVA_CLEAR_OBJ(rs_xT);
#if RCT_SJ_M4
    SNOVA_CLEAR_OBJ(rs_g16);
    SNOVA_CLEAR_OBJ(rs_xe16);
#endif
#undef RCT_RS_AXPY_ROW
#else
#if RCT_HOT_SIMD && SNOVA_l == 4 && SNOVA_r <= 7
        #define RCT_GPW 32
        RCT_SCRATCH _Alignas(32) uint8_t gaussp[SNOVA_o * SNOVA_lr][SNOVA_o][RCT_GPW];
        memset(gaussp, 0, sizeof(gaussp));
#elif RCT_SJ_M4
        #define RCT_GPW (SNOVA_l * 16)
        RCT_SCRATCH _Alignas(32) uint16_t gaussp16[SNOVA_o * SNOVA_lr][SNOVA_o][SNOVA_l * 8];
        memset(gaussp16, 0, sizeof(gaussp16));
#elif RCT_SIGN_JOG
        #define RCT_GPW (SNOVA_l * 16)
        RCT_SCRATCH _Alignas(32) uint8_t gaussp[SNOVA_o * SNOVA_lr][SNOVA_o][RCT_GPW];
        memset(gaussp, 0, sizeof(gaussp));
#endif

        for (int mi = 0; mi < SNOVA_o; mi++)
            for (int idx = 0; idx < SNOVA_o; idx++)
                for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
                    gf_t gf16m_temp0[SNOVA_l2 + 16] = {0};
                    gf_t gf16m_temp1[SNOVA_lr + 16] = {0};
                    gf_t X_tmp[SNOVA_l2 + 16] = {0};
                    int mi_prime = i_prime(mi, alpha);
#ifdef RCT_PROFILE
                    uint64_t _t0 = __rdtsc(), _t1, _t2;
#endif
#if RCT_HOT_SIMD && SNOVA_l == 4 && SNOVA_r <= 7
                    {
                        __m256i acc = _mm256_setzero_si256();
                        for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                            __m256i qv = RCT_BC(q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
                            __m256i wv = _mm256_loadu_si256(
                                (const __m256i *)&whipped_F21[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr]);
                            acc = _mm256_xor_si256(acc, RCT_SV(qv, wv));
                        }
                        _Alignas(32) uint8_t accbuf[32];
                        _mm256_store_si256((__m256i *)accbuf, acc);
                        for (int i1 = 0; i1 < SNOVA_lr; i1++) gf16m_temp1[i1] = rct_gfni_cleanup(accbuf[i1]);
                    }
#elif RCT_SIGN_JOG
                    for (int cc = 0; cc < SNOVA_lr; cc += 16) {
                        __m128i acc = _mm_setzero_si128();
                        for (int b1 = 0; b1 < SNOVA_l; ++b1)
                            acc = _mm_xor_si128(acc, rct_sj_sv128(
                                rct_sj_bc128_pub(q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]),
                                _mm_loadu_si128((const __m128i *)&whipped_F21[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + cc])));
                        _Alignas(16) uint8_t pb[16]; _mm_store_si128((__m128i *)pb, RCT_SJ_CLEAN(acc));
                        int nn = SNOVA_lr - cc; if (nn > 16) nn = 16;
                        for (int j = 0; j < nn; ++j) gf16m_temp1[cc + j] = pb[j];
                    }
#else
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        for (int i1 = 0; i1 < SNOVA_lr; i1++)
                            gf_set_add(&gf16m_temp1[i1],
                                       gf_mult_sec(whipped_F21[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1],
                                               q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]));
#endif
#if RCT_USE_SIMD && SNOVA_l == 4
                    rct_matmul_l4rows_sec(gf16m_temp0, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_l, SNOVA_r);
#elif RCT_SJ_M4
                    rct_m4_mm_add(gf16m_temp0, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_l, SNOVA_r, SNOVA_l);
#elif RCT_SIGN_JOG
                    rct_sj_mm_add(gf16m_temp0, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_l, SNOVA_r, SNOVA_l);
#else
                    gf_mat_mul_add_lr_sec(gf16m_temp0, gf16m_temp1, &Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr], SNOVA_l, SNOVA_r, SNOVA_l);
#endif
#if RCT_SIGN_JOG
                    rct_sj_mm_add_pub(X_tmp, &Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2], gf16m_temp0, SNOVA_l, SNOVA_l, SNOVA_l);
#else
                    RCT_MATMUL_ADD_BSEC(X_tmp, &Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2], gf16m_temp0);
#endif
#ifdef RCT_PROFILE
                    _t1 = __rdtsc(); RCT_PACC6B(0, _t0, _t1);
#endif
#if RCT_HOT_SIMD && SNOVA_l == 4 && SNOVA_r <= 7
                    const __m128i _m7 = _mm_set_epi64x(0LL, 0x00FFFFFFFFFFFFFFLL);
                    for (int ti1 = 0; ti1 < SNOVA_r; ti1++) {
                        __m128i amv = _mm_and_si128(_mm_loadu_si128(
                            (const __m128i *)&Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r]), _m7);
                        for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                            __m128i p0 = RCT_SV128(RCT_BC128_SEC(X_tmp[0 * SNOVA_l + ti2]), amv);
                            __m128i p1 = RCT_SV128(RCT_BC128_SEC(X_tmp[1 * SNOVA_l + ti2]), amv);
                            __m128i p2 = RCT_SV128(RCT_BC128_SEC(X_tmp[2 * SNOVA_l + ti2]), amv);
                            __m128i p3 = RCT_SV128(RCT_BC128_SEC(X_tmp[3 * SNOVA_l + ti2]), amv);
                            __m256i comb = _mm256_set_m128i(_mm_unpacklo_epi64(p2, p3), _mm_unpacklo_epi64(p0, p1));
                            __m256i *g = (__m256i *)&gaussp[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx][0];
                            *g = _mm256_xor_si256(*g, comb);
                        }
                    }
#elif RCT_SJ_M4
                    {
                        _Alignas(32) uint16_t xe[SNOVA_l2 + 16];
                        rct_m4_expand_buf(xe, X_tmp, SNOVA_l2);
                        for (int ti1 = 0; ti1 < SNOVA_r; ti1++) {
                            __m128i amv = _mm256_castsi256_si128(_mm256_cvtepu8_epi16(
                                _mm_loadu_si128((const __m128i *)
                                    &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r])));
                            for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                                __m128i *gs = (__m128i *)gaussp16[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx];
                                for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                                    gs[tj1] = _mm_xor_si128(gs[tj1], _mm_mullo_epi16(
                                        _mm_set1_epi16((short)xe[tj1 * SNOVA_l + ti2]), amv));
                            }
                        }
                    }
#elif RCT_SIGN_JOG
                    {
                        __m128i bcx[SNOVA_l * SNOVA_l];
                        for (int k = 0; k < SNOVA_l * SNOVA_l; ++k) bcx[k] = rct_sj_bc128(X_tmp[k]);
                        for (int ti1 = 0; ti1 < SNOVA_r; ti1++) {
                            _Alignas(16) uint8_t amr[16] = {0};
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                                amr[tj2] = Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r + tj2];
                            __m128i amv = _mm_load_si128((const __m128i *)amr);
                            for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                                __m128i *gs = (__m128i *)gaussp[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx];
                                for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                                    gs[tj1] = _mm_xor_si128(gs[tj1], rct_sj_sv128(bcx[tj1 * SNOVA_l + ti2], amv));
                            }
                        }
                    }
#else
                    for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
                        for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
                            for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                                for (int tj2 = 0; tj2 < SNOVA_r; tj2++) {
                                    int ti = ti1 * SNOVA_l + ti2;
                                    int tj = tj1 * SNOVA_r + tj2;
                                    gf_set_add(&gauss[mi * SNOVA_lr + ti][idx * SNOVA_lr + tj],
                                               gf_mult_sec(X_tmp[tj1 * SNOVA_l + ti2],
                                                       Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + ti1 * SNOVA_r + tj2]));
                                }
#endif
#ifdef RCT_PROFILE
                    _t2 = __rdtsc(); RCT_PACC6B(1, _t1, _t2);
#endif
                }

        for (int mi = 0; mi < SNOVA_o; mi++)
            for (int idx = 0; idx < SNOVA_o; idx++)
                for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
                    gf_t gf16m_temp0[SNOVA_lr + 16] = {0};
                    gf_t gf16m_temp1[SNOVA_lr + 16] = {0};
                    gf_t X_tmp[SNOVA_lr + 16] = {0};
                    int mi_prime = i_prime(mi, alpha);
#ifdef RCT_PROFILE
                    uint64_t _t0 = __rdtsc(), _t1, _t2;
#endif
#if RCT_HOT_SIMD && SNOVA_l == 4 && SNOVA_r <= 7
                    {
                        __m256i acc = _mm256_setzero_si256();
                        for (int b1 = 0; b1 < SNOVA_l; ++b1) {
                            __m256i qv = RCT_BC(q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
                            __m256i wv = _mm256_loadu_si256(
                                (const __m256i *)&whipped_F12[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr]);
                            acc = _mm256_xor_si256(acc, RCT_SV(qv, wv));
                        }
                        _Alignas(32) uint8_t accbuf[32];
                        _mm256_store_si256((__m256i *)accbuf, acc);
                        for (int i1 = 0; i1 < SNOVA_l; i1++)
                            for (int j1 = 0; j1 < SNOVA_r; j1++)
                                gf16m_temp1[j1 * SNOVA_l + i1] = rct_gfni_cleanup(accbuf[i1 * SNOVA_r + j1]);
                    }
#elif RCT_SIGN_JOG
                    {
                        _Alignas(16) uint8_t fbuf[SNOVA_lr + 16];
                        for (int cc = 0; cc < SNOVA_lr; cc += 16) {
                            __m128i acc = _mm_setzero_si128();
                            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                                acc = _mm_xor_si128(acc, rct_sj_sv128(
                                    rct_sj_bc128_pub(q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]),
                                    _mm_loadu_si128((const __m128i *)&whipped_F12[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + cc])));
                            _mm_store_si128((__m128i *)&fbuf[cc], RCT_SJ_CLEAN(acc));
                        }
                        for (int i1 = 0; i1 < SNOVA_l; ++i1)
                            for (int j1 = 0; j1 < SNOVA_r; ++j1)
                                gf16m_temp1[j1 * SNOVA_l + i1] = fbuf[i1 * SNOVA_r + j1];
                    }
#else
                    for (int b1 = 0; b1 < SNOVA_l; ++b1)
                        for (int i1 = 0; i1 < SNOVA_l; i1++)
                            for (int j1 = 0; j1 < SNOVA_r; j1++)
                                gf_set_add(&gf16m_temp1[j1 * SNOVA_l + i1],
                                           gf_mult_sec(whipped_F12[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_lr + i1 * SNOVA_r + j1],
                                                   q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]));
#endif
#if RCT_USE_SIMD && SNOVA_l == 4
                    rct_matmul_l4rows(gf16m_temp0, &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp1, SNOVA_r, SNOVA_r);
                    rct_matmul_l4rows_sec(X_tmp, gf16m_temp0, &Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2], SNOVA_r, SNOVA_l);
#elif RCT_SIGN_JOG
                    rct_sj_mm_add_pub(gf16m_temp0, &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp1, SNOVA_r, SNOVA_r, SNOVA_l);
#if RCT_SJ_M4
                    rct_m4_mm_add(X_tmp, gf16m_temp0, &Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2], SNOVA_r, SNOVA_l, SNOVA_l);
#else
                    rct_sj_mm_add(X_tmp, gf16m_temp0, &Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2], SNOVA_r, SNOVA_l, SNOVA_l);
#endif
#else
                    gf_mat_mul_add_lr_sec(gf16m_temp0, &Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2], gf16m_temp1, SNOVA_r, SNOVA_r, SNOVA_l);
                    gf_mat_mul_add_lr_sec(X_tmp, gf16m_temp0, &Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2], SNOVA_r, SNOVA_l, SNOVA_l);
#endif
#ifdef RCT_PROFILE
                    _t1 = __rdtsc(); RCT_PACC6B(0, _t0, _t1);
#endif
#if RCT_HOT_SIMD && SNOVA_l == 4 && SNOVA_r <= 7
                    __m128i bmcol[SNOVA_l];
                    for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                        _Alignas(16) uint8_t cc[16] = {0};
                        for (int tj2 = 0; tj2 < SNOVA_r; tj2++)
                            cc[tj2] = Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2];
                        bmcol[ti2] = _mm_load_si128((const __m128i *)cc);
                    }
                    for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
                        for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                            __m128i b = bmcol[ti2];
                            __m128i p0 = RCT_SV128(RCT_BC128_SEC(X_tmp[ti1 * SNOVA_l + 0]), b);
                            __m128i p1 = RCT_SV128(RCT_BC128_SEC(X_tmp[ti1 * SNOVA_l + 1]), b);
                            __m128i p2 = RCT_SV128(RCT_BC128_SEC(X_tmp[ti1 * SNOVA_l + 2]), b);
                            __m128i p3 = RCT_SV128(RCT_BC128_SEC(X_tmp[ti1 * SNOVA_l + 3]), b);
                            __m256i comb = _mm256_set_m128i(_mm_unpacklo_epi64(p2, p3), _mm_unpacklo_epi64(p0, p1));
                            __m256i *g = (__m256i *)&gaussp[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx][0];
                            *g = _mm256_xor_si256(*g, comb);
                        }
#elif RCT_SJ_M4
                    {
                        _Alignas(32) uint16_t xe2[SNOVA_lr + 16];
                        rct_m4_expand_buf(xe2, X_tmp, SNOVA_lr);
                        for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                            _Alignas(16) uint8_t bmc[16] = {0};
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                                bmc[tj2] = Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2];
                            __m128i bmv = _mm256_castsi256_si128(_mm256_cvtepu8_epi16(
                                _mm_load_si128((const __m128i *)bmc)));
                            for (int ti1 = 0; ti1 < SNOVA_r; ti1++) {
                                __m128i *gs = (__m128i *)gaussp16[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx];
                                for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                                    gs[tj1] = _mm_xor_si128(gs[tj1], _mm_mullo_epi16(
                                        _mm_set1_epi16((short)xe2[ti1 * SNOVA_l + tj1]), bmv));
                            }
                        }
                    }
#elif RCT_SIGN_JOG
                    {
                        __m128i bcx[SNOVA_lr];
                        for (int k = 0; k < SNOVA_lr; ++k) bcx[k] = rct_sj_bc128(X_tmp[k]);
                        for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
                            _Alignas(16) uint8_t bmc[16] = {0};
                            for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                                bmc[tj2] = Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2];
                            __m128i bmv = _mm_load_si128((const __m128i *)bmc);
                            for (int ti1 = 0; ti1 < SNOVA_r; ti1++) {
                                __m128i *gs = (__m128i *)gaussp[mi * SNOVA_lr + ti1 * SNOVA_l + ti2][idx];
                                for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                                    gs[tj1] = _mm_xor_si128(gs[tj1], rct_sj_sv128(bcx[ti1 * SNOVA_l + tj1], bmv));
                            }
                        }
                    }
#else
                    for (int ti1 = 0; ti1 < SNOVA_r; ti1++)
                        for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
                            for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
                                for (int tj2 = 0; tj2 < SNOVA_r; tj2++) {
                                    int ti = ti1 * SNOVA_l + ti2;
                                    int tj = tj1 * SNOVA_r + tj2;
                                    gf_set_add(&gauss[mi * SNOVA_lr + ti][idx * SNOVA_lr + tj],
                                               gf_mult_sec(X_tmp[ti1 * SNOVA_l + tj1],
                                                       Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + tj2 * SNOVA_l + ti2]));
                                }
#endif
#ifdef RCT_PROFILE
                    _t2 = __rdtsc(); RCT_PACC6B(1, _t1, _t2);
#endif
                }

#if RCT_HOT_SIMD && SNOVA_l == 4 && SNOVA_r <= 7
#ifdef RCT_PROFILE
        uint64_t _tc0 = __rdtsc(), _tc1;
#endif
        for (int rr = 0; rr < SNOVA_o * SNOVA_lr; ++rr)
            for (int ix = 0; ix < SNOVA_o; ++ix)
                for (int tj1 = 0; tj1 < SNOVA_l; ++tj1)
                    for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                        gauss[rr][ix * SNOVA_lr + tj1 * SNOVA_r + tj2] =
                            rct_gfni_cleanup(gaussp[rr][ix][tj1 * 8 + tj2]);
#ifdef RCT_PROFILE
        _tc1 = __rdtsc(); RCT_PACC6B(2, _tc0, _tc1);
#endif
        #undef RCT_GPW
#elif RCT_SJ_M4
#ifdef RCT_PROFILE
        uint64_t _tc0m = __rdtsc(), _tc1m;
#endif
        for (int rr = 0; rr < SNOVA_o * SNOVA_lr; ++rr)
            for (int ix = 0; ix < SNOVA_o; ++ix)
                for (int tj1 = 0; tj1 < SNOVA_l; ++tj1) {
                    __m128i c16 = rct_m4_compress128(
                        _mm_load_si128((const __m128i *)gaussp16[rr][ix] + tj1));
                    _Alignas(16) uint8_t cb[16];
                    _mm_store_si128((__m128i *)cb, _mm_packus_epi16(c16, _mm_setzero_si128()));
                    for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                        gauss[rr][ix * SNOVA_lr + tj1 * SNOVA_r + tj2] = cb[tj2];
                }
#ifdef RCT_PROFILE
        _tc1m = __rdtsc(); RCT_PACC6B(2, _tc0m, _tc1m);
#endif
        #undef RCT_GPW
#elif RCT_SIGN_JOG
#ifdef RCT_PROFILE
        uint64_t _tc0 = __rdtsc(), _tc1;
#endif
        for (int rr = 0; rr < SNOVA_o * SNOVA_lr; ++rr)
            for (int ix = 0; ix < SNOVA_o; ++ix)
                for (int tj1 = 0; tj1 < SNOVA_l; ++tj1) {
                    _Alignas(16) uint8_t cb[16];
                    _mm_store_si128((__m128i *)cb, RCT_SJ_CLEAN(
                        _mm_load_si128((const __m128i *)&gaussp[rr][ix][tj1 * 16])));
                    for (int tj2 = 0; tj2 < SNOVA_r; ++tj2)
                        gauss[rr][ix * SNOVA_lr + tj1 * SNOVA_r + tj2] = cb[tj2];
                }
#ifdef RCT_PROFILE
        _tc1 = __rdtsc(); RCT_PACC6B(2, _tc0, _tc1);
#endif
        #undef RCT_GPW
#endif
#if RCT_HOT_SIMD && SNOVA_l == 4 && SNOVA_r <= 7
    SNOVA_CLEAR_OBJ(gaussp);
#elif RCT_SJ_M4
    SNOVA_CLEAR_OBJ(gaussp16);
#elif RCT_SIGN_JOG
    SNOVA_CLEAR_OBJ(gaussp);
#endif
#endif
}
#endif

#if !RCT_SIGN_JOG && !RCT_OQDF && !(RCT_USE_SIMD && SNOVA_r <= 16) && !RCT_OQWV
static void rct_sign_whipbuild_scalar(rct_sign_ctx *c) {
    const gf_t *signature_in_GF = c->signature_in_GF;
    gf_t *whipped_sig = c->whipped_sig;
    for (int ab = 0; ab < SNOVA_l; ++ab)
        for (int ni = 0; ni < SNOVA_v; ++ni)
            for (int i1 = 0; i1 < SNOVA_l; i1++)
                for (int j1 = 0; j1 < SNOVA_r; j1++)
                    for (int k1 = 0; k1 < SNOVA_l; k1++)
                        gf_set_add(&whipped_sig[(ab * SNOVA_v + ni) * SNOVA_lr + i1 * SNOVA_r + j1],
                                   gf_mult_sec(rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1],
                                               signature_in_GF[ni * SNOVA_lr + k1 * SNOVA_r + j1]));
}
#endif

#if !RCT_SIGN_JOG && !RCT_USE_SIMD && !(RCT_Q_SIMD && (SNOVA_r == SNOVA_l)) && !(RCT_Q_SIMD && RCT_Q_HAVE_MAGIC)
static void rct_sign_sumt_scalar(rct_sign_ctx *c) {
    const gf_t *P11 = c->P11, *whipped_sig = c->whipped_sig;
    gf_t *sum_t1 = c->sum_t1;
    RCT_SCRATCH _Alignas(32) gf_t sum_t0[SNOVA_m1 * SNOVA_l * SNOVA_v * SNOVA_lr];
    memset(sum_t0, 0, sizeof(sum_t0));
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        for (int ni = 0; ni < SNOVA_v; ++ni)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int nj = 0; nj < SNOVA_v; ++nj)
                    gf_mat_mul_add_lr_sec(&sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr],
                                      &P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2],
                                      &whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_lr], SNOVA_l, SNOVA_l, SNOVA_r);

        for (int a1 = 0; a1 < SNOVA_l; ++a1)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int ni = 0; ni < SNOVA_v; ++ni)
                    for (int i1 = 0; i1 < SNOVA_r; i1++)
                        for (int j1 = 0; j1 < SNOVA_r; j1++)
                            for (int k1 = 0; k1 < SNOVA_l; k1++)
                                gf_set_add(&sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1],
                                           gf_mult_sec(whipped_sig[(a1 * SNOVA_v + ni) * SNOVA_lr + k1 * SNOVA_r + i1],
                                                   sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_lr + k1 * SNOVA_r + j1]));
    }
    SNOVA_CLEAR_OBJ(sum_t0);
}
#endif

#if !RCT_SIGN_JOG && !RCT_Q_SIMD && !(RCT_USE_GFNI && (SNOVA_l == 4)) && !(RCT_HOT_QRP16 && (SNOVA_l == 4)) \
    && !(RCT_USE_PSHUFB && (SNOVA_l == 4) && defined(RCT_MULLO) && (RCT_MULLO + 0))
static void rct_skx_fold_F_scalar(gf_t *F21, gf_t *F12, const gf_t *T12, const gf_t *P11) {
    for (int i1 = 0; i1 < SNOVA_m1; i1++)
        for (int j1 = 0; j1 < SNOVA_v; j1++)
            for (int j2 = 0; j2 < SNOVA_v; j2++)
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    RCT_MATMUL_ADD_ASEC(&F21[((i1 * SNOVA_o + k1) * SNOVA_v + j1) * SNOVA_l2],
                                   &T12[(j2 * SNOVA_o + k1) * SNOVA_l2],
                                   &P11[((i1 * SNOVA_v + j2) * SNOVA_v + j1) * SNOVA_l2]);

    for (int i1 = 0; i1 < SNOVA_m1; i1++)
        for (int j1 = 0; j1 < SNOVA_v; j1++)
            for (int j2 = 0; j2 < SNOVA_v; j2++)
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    RCT_MATMUL_ADD_BSEC(&F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
                                   &P11[((i1 * SNOVA_v + j1) * SNOVA_v + j2) * SNOVA_l2],
                                   &T12[(j2 * SNOVA_o + k1) * SNOVA_l2]);
}
#endif

#endif
