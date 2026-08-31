#ifndef RCT_PKX_H
#define RCT_PKX_H

#ifndef RCT_PKX_JOGSIMD
#define RCT_PKX_JOGSIMD 1
#endif
#if RCT_JOG_PKXJOG && RCT_PKX_JOGSIMD && (SNOVA_l >= 4) && (SNOVA_l <= 7)
#define RCT_PKX_JOG_SIMD 1
static inline void rct_pkx_jog_scatter(gf_t *P, size_t mi_row_base, int ni,
                                       size_t col0, const gf_t *src, int K) {
    gf_t *row[SNOVA_l];
    for (int ei = 0; ei < SNOVA_l; ++ei)
        row[ei] = P + (mi_row_base + (size_t)ni * SNOVA_l + ei) * RCT_JOG_L32 + col0;
    for (int nj = 0; nj + 1 < K; ++nj) {
        const gf_t *s = src + (size_t)nj * SNOVA_l2;
        for (int ei = 0; ei < SNOVA_l; ++ei)
            memcpy(row[ei] + (size_t)nj * SNOVA_l, s + (size_t)ei * SNOVA_l, 8);
    }
    const gf_t *s = src + (size_t)(K - 1) * SNOVA_l2;
    for (int ei = 0; ei < SNOVA_l; ++ei)
        memcpy(row[ei] + (size_t)(K - 1) * SNOVA_l, s + (size_t)ei * SNOVA_l, SNOVA_l);
}
#else
#define RCT_PKX_JOG_SIMD 0
#endif

static int rct_pkx_expand(rct_pkx_ctx *c) {
    rct_pk_t *pkx = c->pkx;
    const uint8_t *pk = c->pk;
#if RCT_VERIFY_STREAM
    memset(pkx, 0, sizeof(*pkx));
    memcpy(pkx->pk_seed, pk, SEED_LENGTH_PUBLIC);
#if HASH_PK
    shake256(pkx->pk_hash, BYTES_PK_HASH, pk, BYTES_PK);
#endif
    memcpy(pkx->Am, rct_fixed_Am, sizeof(pkx->Am));
    memcpy(pkx->Bm, rct_fixed_Bm, sizeof(pkx->Bm));
    memcpy(pkx->q1, rct_fixed_abq + (size_t)SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr),
           SNOVA_o * SNOVA_alpha * SNOVA_l);
    memcpy(pkx->q2, rct_fixed_abq + (size_t)SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + SNOVA_l),
           SNOVA_o * SNOVA_alpha * SNOVA_l);
#if SNOVA_q != 16
    if (expand_pk(pkx->P22, pk + SEED_LENGTH_PUBLIC)) return -1;
#else
    (void)expand_pk;
    rct_unpack_nib_seg(pkx->P22, pk + SEED_LENGTH_PUBLIC, NUMGF_PK / 2);
#endif
    return 0;
#elif RCT_PKX_FUSED
    memset((uint8_t *)pkx + offsetof(rct_pk_t, Am), 0, sizeof(*pkx) - offsetof(rct_pk_t, Am));
    memcpy(pkx->pk_seed, pk, SEED_LENGTH_PUBLIC);
#if HASH_PK
    shake256(pkx->pk_hash, BYTES_PK_HASH, pk, BYTES_PK);
#endif

    memcpy(pkx->Am, rct_fixed_Am, sizeof(pkx->Am));
    memcpy(pkx->Bm, rct_fixed_Bm, sizeof(pkx->Bm));
    memcpy(pkx->q1, rct_fixed_abq + (size_t)SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr),
           SNOVA_o * SNOVA_alpha * SNOVA_l);
    memcpy(pkx->q2, rct_fixed_abq + (size_t)SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + SNOVA_l),
           SNOVA_o * SNOVA_alpha * SNOVA_l);

#if !RCT_PKX_PGEN
#if SNOVA_WRAPPER_STACK
    _Alignas(32) uint8_t pub_bytes[((NUM_GEN_PUB_BYTES + 15) & ~(size_t)7u)];
#else
    _Alignas(8) static uint8_t pub_bytes[((NUM_GEN_PUB_BYTES + 15) & ~(size_t)7u)];
#endif
    rct_public_xof(pk, pub_bytes, NUM_GEN_PUB_BYTES);
#endif

#if SNOVA_q != 16
#if SNOVA_WRAPPER_STACK
    _Alignas(32) gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
#else
    static gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
#endif
    if (expand_pk(P22, pk + SEED_LENGTH_PUBLIC)) return -1;
#else
    (void)expand_pk;
#endif

    enum {
        RCT_OFF_P11 = 0,
        RCT_OFF_P12 = SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2,
        RCT_OFF_P21 = SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2,
    };
#if SNOVA_q == 16
#if !RCT_PKX_PGEN
#define RCT_PKX_SEG(dst, gfoff, ngf) \
    rct_unpack_nib_seg((dst), pub_bytes + ((size_t)(gfoff) >> 1), (size_t)(ngf) >> 1)
#endif
#define RCT_PKX_P22SEG(dst, gfoff, ngf) \
    rct_unpack_nib_seg((dst), pk + SEED_LENGTH_PUBLIC + ((size_t)(gfoff) >> 1), (size_t)(ngf) >> 1)
#else
#if !RCT_PKX_PGEN
#define RCT_PKX_SEG(dst, gfoff, ngf) rct_unpack_modq_seg((dst), pub_bytes + (size_t)(gfoff), (size_t)(ngf))
#endif
#define RCT_PKX_P22SEG(dst, gfoff, ngf) memcpy((dst), P22 + (size_t)(gfoff), (size_t)(ngf))
#endif
#if RCT_PKX_PGEN
    {
        snova_pgen_t pg;
        snova_pgen_init(&pg, pk);
        snova_prow_t rw;
        int blk, bmi, bni, bnc;
        (void)bnc;
        while (snova_pgen_peek(&pg, &blk, &bmi, &bni, &bnc)) {
            gf_t *dst;
            if (blk == 0)
                dst = &pkx->P[((size_t)(bmi * SNOVA_n + bni) * SNOVA_n) * SNOVA_l2];
            else if (blk == 1)
                dst = &pkx->P[(((size_t)(bmi * SNOVA_n + bni) * SNOVA_n) + SNOVA_v) * SNOVA_l2];
            else
                dst = &pkx->P[((size_t)(bmi * SNOVA_n + SNOVA_v + bni) * SNOVA_n) * SNOVA_l2];
            (void)snova_pgen_next_row_into(&pg, &rw, dst);
        }
    }
    for (int mi = 0; mi < SNOVA_m1; ++mi)
        for (int ni = 0; ni < SNOVA_o; ++ni)
            RCT_PKX_P22SEG(&pkx->P[(((size_t)(mi * SNOVA_n + SNOVA_v + ni) * SNOVA_n) + SNOVA_v) * SNOVA_l2],
                           (size_t)(mi * SNOVA_o + ni) * SNOVA_o * SNOVA_l2, SNOVA_o * SNOVA_l2);
#else
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        for (int ni = 0; ni < SNOVA_v; ++ni) {
            gf_t *row = &pkx->P[((size_t)(mi * SNOVA_n + ni) * SNOVA_n) * SNOVA_l2];
            RCT_PKX_SEG(row, RCT_OFF_P11 + (size_t)(mi * SNOVA_v + ni) * SNOVA_v * SNOVA_l2,
                        SNOVA_v * SNOVA_l2);
            RCT_PKX_SEG(row + (size_t)SNOVA_v * SNOVA_l2,
                        RCT_OFF_P12 + (size_t)(mi * SNOVA_v + ni) * SNOVA_o * SNOVA_l2,
                        SNOVA_o * SNOVA_l2);
        }
        for (int ni = 0; ni < SNOVA_o; ++ni) {
            gf_t *row = &pkx->P[((size_t)(mi * SNOVA_n + SNOVA_v + ni) * SNOVA_n) * SNOVA_l2];
            RCT_PKX_SEG(row, RCT_OFF_P21 + (size_t)(mi * SNOVA_o + ni) * SNOVA_v * SNOVA_l2,
                        SNOVA_v * SNOVA_l2);
            RCT_PKX_P22SEG(row + (size_t)SNOVA_v * SNOVA_l2,
                           (size_t)(mi * SNOVA_o + ni) * SNOVA_o * SNOVA_l2, SNOVA_o * SNOVA_l2);
        }
    }
#undef RCT_PKX_SEG
#endif
#undef RCT_PKX_P22SEG
    return 0;
#else
    memset(pkx, 0, sizeof(*pkx));
    memcpy(pkx->pk_seed, pk, SEED_LENGTH_PUBLIC);
#if HASH_PK
    shake256(pkx->pk_hash, BYTES_PK_HASH, pk, BYTES_PK);
#endif

#if SNOVA_WRAPPER_STACK
    _Alignas(32) gf_t P_matrix[NUM_PUB_GF];
    _Alignas(32) gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
#else
    gf_t *P_matrix = rct_pub_Pmatrix;
    static gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
#endif
    gf_t *P11 = P_matrix;
    gf_t *P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
    gf_t *P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;

    if (expand_pk(P22, pk + SEED_LENGTH_PUBLIC)) return -1;
    expand_public(P_matrix, pk);

#if RCT_PKX_JOG_SIMD
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        const size_t mrb = (size_t)mi * RCT_JOG_NL;
        for (int ni = 0; ni < SNOVA_v; ++ni) {
            rct_pkx_jog_scatter(pkx->P, mrb, ni, 0,
                                P11 + ((size_t)(mi * SNOVA_v + ni) * SNOVA_v) * SNOVA_l2, SNOVA_v);
            rct_pkx_jog_scatter(pkx->P, mrb, ni, (size_t)SNOVA_v * SNOVA_l,
                                P12 + ((size_t)(mi * SNOVA_v + ni) * SNOVA_o) * SNOVA_l2, SNOVA_o);
        }
        for (int ni = SNOVA_v; ni < SNOVA_n; ++ni) {
            const int nio = ni - SNOVA_v;
            rct_pkx_jog_scatter(pkx->P, mrb, ni, 0,
                                P21 + ((size_t)(mi * SNOVA_o + nio) * SNOVA_v) * SNOVA_l2, SNOVA_v);
            rct_pkx_jog_scatter(pkx->P, mrb, ni, (size_t)SNOVA_v * SNOVA_l,
                                P22 + ((size_t)(mi * SNOVA_o + nio) * SNOVA_o) * SNOVA_l2, SNOVA_o);
        }
    }
#else
#if RCT_JOG_PKXJOG
#define RCT_PKX_DST(mi, ni, nj, idx)                                              \
    pkx->P[(((size_t)(mi) * RCT_JOG_NL + (size_t)(ni) * SNOVA_l + (idx) / SNOVA_l) \
            * RCT_JOG_L32) + (size_t)(nj) * SNOVA_l + (idx) % SNOVA_l]
#else
#define RCT_PKX_DST(mi, ni, nj, idx) \
    pkx->P[(((mi) * SNOVA_n + (ni)) * SNOVA_n + (nj)) * SNOVA_l2 + (idx)]
#endif
    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        for (int ni = 0; ni < SNOVA_v; ++ni) {
            for (int nj = 0; nj < SNOVA_v; ++nj)
                for (int idx = 0; idx < SNOVA_l2; idx++)
                    RCT_PKX_DST(mi, ni, nj, idx) =
                        P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + idx];
            for (int nj = SNOVA_v; nj < SNOVA_n; ++nj)
                for (int idx = 0; idx < SNOVA_l2; idx++)
                    RCT_PKX_DST(mi, ni, nj, idx) =
                        P12[((mi * SNOVA_v + ni) * SNOVA_o + (nj - SNOVA_v)) * SNOVA_l2 + idx];
        }
        for (int ni = SNOVA_v; ni < SNOVA_n; ++ni) {
            for (int nj = 0; nj < SNOVA_v; ++nj)
                for (int idx = 0; idx < SNOVA_l2; idx++)
                    RCT_PKX_DST(mi, ni, nj, idx) =
                        P21[((mi * SNOVA_o + (ni - SNOVA_v)) * SNOVA_v + nj) * SNOVA_l2 + idx];
            for (int nj = SNOVA_v; nj < SNOVA_n; ++nj)
                for (int idx = 0; idx < SNOVA_l2; idx++)
                    RCT_PKX_DST(mi, ni, nj, idx) =
                        P22[((mi * SNOVA_o + (ni - SNOVA_v)) * SNOVA_o + nj - SNOVA_v) * SNOVA_l2 + idx];
        }
    }
#undef RCT_PKX_DST
#endif

    gf_t *A = P_matrix + (SNOVA_m1 * (SNOVA_n * SNOVA_n - SNOVA_o * SNOVA_o)) * SNOVA_l2;
    gf_t *B = A + SNOVA_o * SNOVA_alpha * SNOVA_r2;
    gf_t *q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
    gf_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;
#if FIXED_ABQ
    memcpy(A, rct_fixed_abq, sizeof(rct_fixed_abq));
#endif
    for (size_t idx = 0; idx < (size_t)SNOVA_o * SNOVA_alpha; idx++) {
        be_invertible_by_add_aS(&pkx->Am[idx * SNOVA_r2], &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
        be_invertible_by_add_aS(&pkx->Bm[idx * SNOVA_lr], &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);
#if ROUND2_T12
        if (!q1[idx * SNOVA_l + SNOVA_l - 1])
            q1[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q1[idx * SNOVA_l] + (q1[idx * SNOVA_l] == 0));
        if (!q2[idx * SNOVA_l + SNOVA_l - 1])
            q2[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q2[idx * SNOVA_l] + (q2[idx * SNOVA_l] == 0));
#endif
    }
    memcpy(pkx->q1, q1, SNOVA_o * SNOVA_alpha * SNOVA_l);
    memcpy(pkx->q2, q2, SNOVA_o * SNOVA_alpha * SNOVA_l);
    return 0;
#endif
}

#endif
