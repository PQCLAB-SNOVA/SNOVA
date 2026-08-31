/**
 * @file gf16m.c
 */
#include "gf16m.h"
#include "../gf16_spec.h"

#include <string.h>

static gf16m_t S_mat;
static gf16m_t S_pow[SNOVA_L];
static int ring_initialised = 0;

void gf16m_zero(gf16m_t a) { memset(a, 0, SNOVA_SQ_RANK); }

void gf16m_copy(gf16m_t dst, const gf16m_t src) {
    memcpy(dst, src, SNOVA_SQ_RANK);
}

int gf16m_eq(const gf16m_t a, const gf16m_t b) {
    return memcmp(a, b, SNOVA_SQ_RANK) == 0;
}

void gf16m_identity(gf16m_t a) {
    gf16m_zero(a);
    for (int i = 0; i < SNOVA_RANK; ++i) GF16M_AT(a, i, i) = 1;
}

void gf16m_add(const gf16m_t a, const gf16m_t b, gf16m_t c) {
    for (int i = 0; i < SNOVA_SQ_RANK; ++i) c[i] = gf16_add(a[i], b[i]);
}

void gf16m_scale(const gf16m_t a, gf16_t k, gf16m_t c) {
    for (int i = 0; i < SNOVA_SQ_RANK; ++i) c[i] = gf16_mul(a[i], k);
}

void gf16m_transpose(const gf16m_t a, gf16m_t at) {
    gf16m_t t;
    for (int i = 0; i < SNOVA_RANK; ++i)
        for (int j = 0; j < SNOVA_RANK; ++j)
            GF16M_AT(t, i, j) = GF16M_AT(a, j, i);
    gf16m_copy(at, t);
}

__attribute__((unused))
static gf16_t det_cofactor(const gf16_t *M, int n) {
    if (n == 1) return M[0];
    if (n == 2) return gf16_add(gf16_mul(M[0], M[3]), gf16_mul(M[1], M[2]));
    gf16_t acc = 0;
    gf16_t minor[SNOVA_SQ_RANK];
    for (int col = 0; col < n; ++col) {
        int mi = 0;
        for (int i = 1; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (j == col) continue;
                minor[mi++] = M[i * n + j];
            }
        acc = gf16_add(acc, gf16_mul(M[col], det_cofactor(minor, n - 1)));
    }
    return acc;
}

#include "xgf16.h"

static uint8_t gf16_mul_tab[256];
static int gf16_mul_tab_done = 0;
static void init_gf16_mul_tab(void) {
    if (gf16_mul_tab_done) return;
    for (int a = 0; a < 16; ++a) {
        for (int b = 0; b < 16; ++b) {
            uint8_t p = 0;
            for (int i = 0; i < 4; ++i)
                if (a & (1 << i)) p ^= (b << i);
            for (int i = 6; i >= 4; --i)
                if (p & (1 << i)) p ^= 0x13 << (i - 4);
            gf16_mul_tab[(a << 4) | b] = p & 0xF;
        }
    }
    gf16_mul_tab_done = 1;
}

static inline gf16_t det_gauss_mul(gf16_t a, gf16_t b) {
    return gf16_mul_tab[((a & 0xF) << 4) | (b & 0xF)];
}

static const gf16_t gf16_inv_tab[16] = {
    0,  1,  9, 14, 13, 11,  7,  6,
    15, 2, 12,  5, 10,  4,  3,  8
};
static inline gf16_t det_gauss_inv(gf16_t a) {
    return gf16_inv_tab[a & 0x0F];
}

__attribute__((unused))
static gf16_t det_gauss(const gf16_t *M, int n) {
    gf16_t a[SNOVA_SQ_RANK];
    for (int i = 0; i < n * n; ++i) a[i] = M[i];
    gf16_t det = 1;
    for (int i = 0; i < n; ++i) {
        if (a[i * n + i] == 0) {
            int j;
            for (j = i + 1; j < n; ++j) if (a[j * n + i] != 0) break;
            if (j == n) return 0;
            for (int k = 0; k < n; ++k) {
                gf16_t t = a[i * n + k]; a[i * n + k] = a[j * n + k]; a[j * n + k] = t;
            }
        }
        gf16_t piv = a[i * n + i];
        det = det_gauss_mul(det, piv);
        gf16_t inv = det_gauss_inv(piv);
        for (int j = i + 1; j < n; ++j) {
            gf16_t f = det_gauss_mul(a[j * n + i], inv);
            if (f == 0) continue;
            for (int k = i; k < n; ++k)
                a[j * n + k] = gf16_add(a[j * n + k], det_gauss_mul(f, a[i * n + k]));
        }
    }
    return det;
}

#include "xgf16.h"

#if SNOVA_RANK >= 2 && SNOVA_RANK <= 4
static inline gf16_t gf16_mul_inline(gf16_t a, gf16_t b) {
    return gf16_mul_tab[((a & 0xF) << 4) | (b & 0xF)];
}
#endif

#if SNOVA_RANK == 2
static inline gf16_t det_inline(const gf16m_t a) {
    return gf16_add(
        gf16_mul_inline(GF16M_AT(a, 0, 0), GF16M_AT(a, 1, 1)),
        gf16_mul_inline(GF16M_AT(a, 0, 1), GF16M_AT(a, 1, 0)));
}
#elif SNOVA_RANK == 3
static inline gf16_t det_inline(const gf16m_t a) {
    gf16_t d0 = gf16_mul_inline(GF16M_AT(a, 0, 0),
        gf16_add(gf16_mul_inline(GF16M_AT(a, 1, 1), GF16M_AT(a, 2, 2)),
                 gf16_mul_inline(GF16M_AT(a, 1, 2), GF16M_AT(a, 2, 1))));
    gf16_t d1 = gf16_mul_inline(GF16M_AT(a, 0, 1),
        gf16_add(gf16_mul_inline(GF16M_AT(a, 1, 0), GF16M_AT(a, 2, 2)),
                 gf16_mul_inline(GF16M_AT(a, 1, 2), GF16M_AT(a, 2, 0))));
    gf16_t d2 = gf16_mul_inline(GF16M_AT(a, 0, 2),
        gf16_add(gf16_mul_inline(GF16M_AT(a, 1, 0), GF16M_AT(a, 2, 1)),
                 gf16_mul_inline(GF16M_AT(a, 1, 1), GF16M_AT(a, 2, 0))));
    return gf16_add(gf16_add(d0, d1), d2);
}
#elif SNOVA_RANK == 4
#define POD4(a, i0, j0, i1, j1, i2, j2, i3, j3, i4, j4) \
    gf16_mul_inline(GF16M_AT(a, i0, j0), \
        gf16_add(gf16_mul_inline(GF16M_AT(a, i1, j1), GF16M_AT(a, i2, j2)), \
                 gf16_mul_inline(GF16M_AT(a, i3, j3), GF16M_AT(a, i4, j4))))

static inline gf16_t det_inline(const gf16m_t a) {
    gf16_t d0 = gf16_mul_inline(GF16M_AT(a, 0, 0),
        gf16_add(gf16_add(POD4(a, 1, 1, 2, 2, 3, 3, 2, 3, 3, 2),
                          POD4(a, 1, 2, 2, 1, 3, 3, 2, 3, 3, 1)),
                 POD4(a, 1, 3, 2, 1, 3, 2, 2, 2, 3, 1)));
    gf16_t d1 = gf16_mul_inline(GF16M_AT(a, 0, 1),
        gf16_add(gf16_add(POD4(a, 1, 0, 2, 2, 3, 3, 2, 3, 3, 2),
                          POD4(a, 1, 2, 2, 0, 3, 3, 2, 3, 3, 0)),
                 POD4(a, 1, 3, 2, 0, 3, 2, 2, 2, 3, 0)));
    gf16_t d2 = gf16_mul_inline(GF16M_AT(a, 0, 2),
        gf16_add(gf16_add(POD4(a, 1, 0, 2, 1, 3, 3, 2, 3, 3, 1),
                          POD4(a, 1, 1, 2, 0, 3, 3, 2, 3, 3, 0)),
                 POD4(a, 1, 3, 2, 0, 3, 1, 2, 1, 3, 0)));
    gf16_t d3 = gf16_mul_inline(GF16M_AT(a, 0, 3),
        gf16_add(gf16_add(POD4(a, 1, 0, 2, 1, 3, 2, 2, 2, 3, 1),
                          POD4(a, 1, 1, 2, 0, 3, 2, 2, 2, 3, 0)),
                 POD4(a, 1, 2, 2, 0, 3, 1, 2, 1, 3, 0)));
    return gf16_add(gf16_add(d0, d1), gf16_add(d2, d3));
}
#elif SNOVA_RANK == 5 && 0
static inline gf16_t det3_n5(const gf16m_t a, int j0, int j1, int j2) {
    gf16_t t0 = gf16_mul_inline(GF16M_AT(a, 0, j0),
        gf16_add(gf16_mul_inline(GF16M_AT(a, 1, j1), GF16M_AT(a, 2, j2)),
                 gf16_mul_inline(GF16M_AT(a, 1, j2), GF16M_AT(a, 2, j1))));
    gf16_t t1 = gf16_mul_inline(GF16M_AT(a, 0, j1),
        gf16_add(gf16_mul_inline(GF16M_AT(a, 1, j0), GF16M_AT(a, 2, j2)),
                 gf16_mul_inline(GF16M_AT(a, 1, j2), GF16M_AT(a, 2, j0))));
    gf16_t t2 = gf16_mul_inline(GF16M_AT(a, 0, j2),
        gf16_add(gf16_mul_inline(GF16M_AT(a, 1, j0), GF16M_AT(a, 2, j1)),
                 gf16_mul_inline(GF16M_AT(a, 1, j1), GF16M_AT(a, 2, j0))));
    return gf16_add(gf16_add(t0, t1), t2);
}

static inline gf16_t det2_n5(const gf16m_t a, int j0, int j1) {
    return gf16_add(gf16_mul_inline(GF16M_AT(a, 3, j0), GF16M_AT(a, 4, j1)),
                    gf16_mul_inline(GF16M_AT(a, 3, j1), GF16M_AT(a, 4, j0)));
}

static inline gf16_t det_inline(const gf16m_t a) {
    gf16_t d012 = gf16_mul_inline(det3_n5(a, 0, 1, 2), det2_n5(a, 3, 4));
    gf16_t d013 = gf16_mul_inline(det3_n5(a, 0, 1, 3), det2_n5(a, 2, 4));
    gf16_t d014 = gf16_mul_inline(det3_n5(a, 0, 1, 4), det2_n5(a, 2, 3));
    gf16_t d023 = gf16_mul_inline(det3_n5(a, 0, 2, 3), det2_n5(a, 1, 4));
    gf16_t d024 = gf16_mul_inline(det3_n5(a, 0, 2, 4), det2_n5(a, 1, 3));
    gf16_t d034 = gf16_mul_inline(det3_n5(a, 0, 3, 4), det2_n5(a, 1, 2));
    gf16_t d123 = gf16_mul_inline(det3_n5(a, 1, 2, 3), det2_n5(a, 0, 4));
    gf16_t d124 = gf16_mul_inline(det3_n5(a, 1, 2, 4), det2_n5(a, 0, 3));
    gf16_t d134 = gf16_mul_inline(det3_n5(a, 1, 3, 4), det2_n5(a, 0, 2));
    gf16_t d234 = gf16_mul_inline(det3_n5(a, 2, 3, 4), det2_n5(a, 0, 1));
    return gf16_add(
        gf16_add(gf16_add(gf16_add(d012, d013), gf16_add(d014, d023)),
                 gf16_add(gf16_add(d024, d034), gf16_add(d123, d124))),
        gf16_add(d134, d234));
}
#endif

gf16_t gf16m_det(const gf16m_t a) {
#if SNOVA_RANK == 5
    return det_gauss(a, SNOVA_RANK);
#elif SNOVA_RANK >= 2 && SNOVA_RANK <= 4
    return det_inline(a);
#else
    return det_cofactor(a, SNOVA_RANK);
#endif
}

static void build_S(void) {
    for (int i = 0; i < SNOVA_RANK; ++i)
        for (int j = 0; j < SNOVA_RANK; ++j)
            GF16M_AT(S_mat, i, j) = GF16_S_ENTRY(i, j);
#if SNOVA_L == 5
    GF16M_AT(S_mat, 4, 4) = 9;
#endif
}

void gf16m_ring_init(void) {
    if (ring_initialised) return;
    gf16_field_init();
    init_gf16_mul_tab();
    build_S();
    gf16m_identity(S_pow[0]);
    for (int k = 1; k < SNOVA_L; ++k)
        gf16m_mul(S_pow[k - 1], S_mat, S_pow[k]);
    ring_initialised = 1;
}

const gf16_t *gf16m_S(void) { return S_mat; }
const gf16_t *gf16m_Spow(int k) { return S_pow[k]; }

void gf16m_make_invertible(gf16m_t a) {
    if (gf16m_det(a) != 0) return;
    for (gf16_t t = 1; t < 16; ++t) {
        gf16m_t tS, cand;
        gf16m_scale(S_mat, t, tS);
        gf16m_add(a, tS, cand);
        if (gf16m_det(cand) != 0) { gf16m_copy(a, cand); return; }
    }
}
