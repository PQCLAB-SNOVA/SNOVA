#ifndef RCT_VERIFY_SCALAR_H
#define RCT_VERIFY_SCALAR_H

#if !RCT_USE_SIMD && !RCT_Q_SIMD && !RCT_VF_JOG

static void rct_vf_contract_ref(rct_vf_ctx *c) {
    const rct_pk_t *pkx = c->pkx;
    const gf_t *signature_in_GF = c->sig_gf;
    gf_t *sum_t1 = c->sum_t1;
    static gf_t whipped_sig[SNOVA_l * SNOVA_n * SNOVA_lr];
    memset(whipped_sig, 0, sizeof(whipped_sig));
    for (int ab = 0; ab < SNOVA_l; ++ab)
        for (int idx = 0; idx < SNOVA_n; ++idx)
            for (int i1 = 0; i1 < SNOVA_l; i1++)
                for (int j1 = 0; j1 < SNOVA_r; j1++)
                    for (int k1 = 0; k1 < SNOVA_l; k1++)
                        gf_set_add(&whipped_sig[(ab * SNOVA_n + idx) * SNOVA_lr + i1 * SNOVA_r + j1],
                                   gf_mult(rct_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1],
                                           signature_in_GF[idx * SNOVA_lr + k1 * SNOVA_r + j1]));

    static gf_t sum_t0[SNOVA_m1 * SNOVA_l * SNOVA_n * SNOVA_lr];
    memset(sum_t0, 0, sizeof(sum_t0));

    for (int mi = 0; mi < SNOVA_m1; ++mi) {
        for (int ni = 0; ni < SNOVA_n; ++ni)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int nj = 0; nj < SNOVA_n; ++nj)
                    gf_mat_mul_add_lr(&sum_t0[((mi * SNOVA_l + b1) * SNOVA_n + ni) * SNOVA_lr],
                                      &pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2],
                                      &whipped_sig[(b1 * SNOVA_n + nj) * SNOVA_lr], SNOVA_l, SNOVA_l, SNOVA_r);

        for (int a1 = 0; a1 < SNOVA_l; ++a1)
            for (int b1 = 0; b1 < SNOVA_l; ++b1)
                for (int ni = 0; ni < SNOVA_n; ++ni)
                    for (int i1 = 0; i1 < SNOVA_r; i1++)
                        for (int j1 = 0; j1 < SNOVA_r; j1++)
                            for (int k1 = 0; k1 < SNOVA_l; k1++)
                                gf_set_add(&sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_r2 + i1 * SNOVA_r + j1],
                                           gf_mult(whipped_sig[(a1 * SNOVA_n + ni) * SNOVA_lr + k1 * SNOVA_r + i1],
                                                   sum_t0[((mi * SNOVA_l + b1) * SNOVA_n + ni) * SNOVA_lr + k1 * SNOVA_r + j1]));
    }
}

#endif

#endif
