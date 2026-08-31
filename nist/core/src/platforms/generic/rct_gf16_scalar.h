#ifndef RCT_GF16_SCALAR_H
#define RCT_GF16_SCALAR_H

static inline void gf_mat_mul(gf_t *a, const gf_t *b, const gf_t *c) {
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int j1 = 0; j1 < SNOVA_l; j1++) {
            gf_t sum = 0;
            for (int k1 = 0; k1 < SNOVA_l; k1++)
                gf_set_add(&sum, gf_mult(b[i1 * SNOVA_l + k1], c[k1 * SNOVA_l + j1]));
            a[i1 * SNOVA_l + j1] = sum;
        }
}

static inline void gf_mat_mul_add(gf_t *a, const gf_t *b, const gf_t *c) {
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int j1 = 0; j1 < SNOVA_l; j1++) {
            gf_t sum = 0;
            for (int k1 = 0; k1 < SNOVA_l; k1++)
                gf_set_add(&sum, gf_mult(b[i1 * SNOVA_l + k1], c[k1 * SNOVA_l + j1]));
            gf_set_add(&a[i1 * SNOVA_l + j1], sum);
        }
}

static inline void gf_mat_mul_add_lr(gf_t *a, const gf_t *b, const gf_t *c, int ad, int bd, int cd) {
    for (int i1 = 0; i1 < ad; i1++)
        for (int j1 = 0; j1 < cd; j1++) {
            gf_t sum = 0;
            for (int k1 = 0; k1 < bd; k1++)
                gf_set_add(&sum, gf_mult(b[i1 * bd + k1], c[k1 * cd + j1]));
            gf_set_add(&a[i1 * cd + j1], sum);
        }
}

static inline void gf_mat_mul_add_sec(gf_t *a, const gf_t *b, const gf_t *c) {
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int j1 = 0; j1 < SNOVA_l; j1++) {
            gf_t sum = 0;
            for (int k1 = 0; k1 < SNOVA_l; k1++)
                gf_set_add(&sum, gf_mult_sec(b[i1 * SNOVA_l + k1], c[k1 * SNOVA_l + j1]));
            gf_set_add(&a[i1 * SNOVA_l + j1], sum);
        }
}
static inline void gf_mat_mul_add_lr_sec(gf_t *a, const gf_t *b, const gf_t *c, int ad, int bd, int cd) {
    for (int i1 = 0; i1 < ad; i1++)
        for (int j1 = 0; j1 < cd; j1++) {
            gf_t sum = 0;
            for (int k1 = 0; k1 < bd; k1++)
                gf_set_add(&sum, gf_mult_sec(b[i1 * bd + k1], c[k1 * cd + j1]));
            gf_set_add(&a[i1 * cd + j1], sum);
        }
}

#endif
