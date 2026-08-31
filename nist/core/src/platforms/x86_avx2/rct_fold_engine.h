#ifndef SNOVA_RCT_FOLD_ENGINE_H
#define SNOVA_RCT_FOLD_ENGINE_H

#ifndef RCT_L4G_FOLD
#define RCT_L4G_FOLD (!RCT_Q_SIMD && RCT_USE_GFNI && (SNOVA_l == 4))
#endif
#if RCT_L4G_FOLD

#define RCT_L4G_PATS \
    const __m128i sA0 = _mm_setr_epi8(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12); \
    const __m128i sA1 = _mm_setr_epi8(1,1,1,1,5,5,5,5,9,9,9,9,13,13,13,13); \
    const __m128i sA2 = _mm_setr_epi8(2,2,2,2,6,6,6,6,10,10,10,10,14,14,14,14); \
    const __m128i sA3 = _mm_setr_epi8(3,3,3,3,7,7,7,7,11,11,11,11,15,15,15,15); \
    const __m128i sB0 = _mm_setr_epi8(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3); \
    const __m128i sB1 = _mm_setr_epi8(4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7); \
    const __m128i sB2 = _mm_setr_epi8(8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11); \
    const __m128i sB3 = _mm_setr_epi8(12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15)

static void rct_l4g_fold_bsec(gf_t *C, const gf_t *A, const gf_t *T12, int nrows, int rmw) {
    RCT_L4G_PATS;
    for (int i1 = 0; i1 < SNOVA_m1; i1++)
        for (int j1 = 0; j1 < nrows; j1++) {
            __m128i acc[SNOVA_o];
            for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm_setzero_si128();
            for (int j2 = 0; j2 < SNOVA_v; j2++) {
                __m128i av = _mm_loadu_si128((const __m128i *)&A[((i1 * nrows + j1) * SNOVA_v + j2) * SNOVA_l2]);
                __m128i ak0 = _mm_shuffle_epi8(av, sA0), ak1 = _mm_shuffle_epi8(av, sA1),
                        ak2 = _mm_shuffle_epi8(av, sA2), ak3 = _mm_shuffle_epi8(av, sA3);
                for (int k1 = 0; k1 < SNOVA_o; k1++) {
                    __m128i bv = _mm_loadu_si128((const __m128i *)&T12[(j2 * SNOVA_o + k1) * SNOVA_l2]);
                    __m128i p = RCT_GFMUL128(ak0, _mm_shuffle_epi8(bv, sB0));
                    p = _mm_xor_si128(p, RCT_GFMUL128(ak1, _mm_shuffle_epi8(bv, sB1)));
                    p = _mm_xor_si128(p, RCT_GFMUL128(ak2, _mm_shuffle_epi8(bv, sB2)));
                    p = _mm_xor_si128(p, RCT_GFMUL128(ak3, _mm_shuffle_epi8(bv, sB3)));
                    acc[k1] = _mm_xor_si128(acc[k1], p);
                }
            }
            for (int k1 = 0; k1 < SNOVA_o; k1++) {
                gf_t *c = &C[((i1 * nrows + j1) * SNOVA_o + k1) * SNOVA_l2];
                __m128i r = rct_gfni_cleanup128(acc[k1]);
                if (rmw) r = _mm_xor_si128(r, _mm_loadu_si128((const __m128i *)c));
                _mm_storeu_si128((__m128i *)c, r);
            }
        }
}

static void rct_l4g_fold_asec(gf_t *P22, const gf_t *T12, const gf_t *F12) {
    RCT_L4G_PATS;
    for (int i1 = 0; i1 < SNOVA_m1; i1++)
        for (int j1 = 0; j1 < SNOVA_o; j1++) {
            __m128i acc[SNOVA_o];
            for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm_setzero_si128();
            for (int idx = 0; idx < SNOVA_v; idx++) {
                __m128i av = _mm_loadu_si128((const __m128i *)&T12[(idx * SNOVA_o + j1) * SNOVA_l2]);
                __m128i ak0 = _mm_shuffle_epi8(av, sA0), ak1 = _mm_shuffle_epi8(av, sA1),
                        ak2 = _mm_shuffle_epi8(av, sA2), ak3 = _mm_shuffle_epi8(av, sA3);
                for (int k1 = 0; k1 < SNOVA_o; k1++) {
                    __m128i bv = _mm_loadu_si128((const __m128i *)&F12[((i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]);
                    __m128i p = RCT_GFMUL128(ak0, _mm_shuffle_epi8(bv, sB0));
                    p = _mm_xor_si128(p, RCT_GFMUL128(ak1, _mm_shuffle_epi8(bv, sB1)));
                    p = _mm_xor_si128(p, RCT_GFMUL128(ak2, _mm_shuffle_epi8(bv, sB2)));
                    p = _mm_xor_si128(p, RCT_GFMUL128(ak3, _mm_shuffle_epi8(bv, sB3)));
                    acc[k1] = _mm_xor_si128(acc[k1], p);
                }
            }
            for (int k1 = 0; k1 < SNOVA_o; k1++) {
                gf_t *c = &P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2];
                _mm_storeu_si128((__m128i *)c,
                    _mm_xor_si128(rct_gfni_cleanup128(acc[k1]),
                                  _mm_loadu_si128((const __m128i *)c)));
            }
        }
}

#endif

#if RCT_SIGN_JOG
#ifndef RCT_FOLD_SEG
#define RCT_FOLD_SEG 1
#endif
#ifndef RCT_FOLD_SEG_M4
#define RCT_FOLD_SEG_M4 1
#endif
#define RCT_F5_A4 (RCT_SIGN_JOG && RCT_HAVE_GFNI  && SNOVA_l == 5 && RCT_FOLD_SEG)
#define RCT_F5_M4 (RCT_SIGN_JOG && !RCT_HAVE_GFNI && SNOVA_l == 5 && RCT_FOLD_SEG && RCT_FOLD_SEG_M4)
#ifndef RCT_FOLD_WIDE
#define RCT_FOLD_WIDE SNOVA_WRAPPER_STACK
#endif
#define RCT_F5_WIDE (RCT_SIGN_JOG && !RCT_HAVE_GFNI && SNOVA_l == 5 && RCT_FOLD_WIDE)

#if RCT_F5_WIDE
#define RCT_WIDE_MVL       (SNOVA_m1 * SNOVA_v * SNOVA_l)
#define RCT_WIDE_P11AW_LEN (RCT_WIDE_MVL * SNOVA_v * SNOVA_l)
#define RCT_WIDE_FW_LEN    (RCT_WIDE_MVL * SNOVA_o * SNOVA_l + 16)
#if SNOVA_WRAPPER_STACK
#define RCT_WIDE_SCRATCH_DECL \
    _Alignas(32) uint16_t P11aw[RCT_WIDE_P11AW_LEN]; \
    _Alignas(32) uint16_t Fw[RCT_WIDE_FW_LEN]
#else
static _Alignas(32) uint16_t rct_wide_p11aw_s[RCT_WIDE_P11AW_LEN];
static _Alignas(32) uint16_t rct_wide_fw_s[RCT_WIDE_FW_LEN];
#define RCT_WIDE_SCRATCH_DECL \
    uint16_t *const P11aw = rct_wide_p11aw_s; \
    uint16_t *const Fw = rct_wide_fw_s
#endif
#endif

#if RCT_F5_WIDE
#include "gf16_core/gf16_mullo16.h"
static inline __attribute__((always_inline)) void rct_wide_fold_bsec(
    gf_t *C, const gf_t *A, const gf_t *T12, int nrows, int rmw,
    uint16_t *const P11aw, uint16_t *const Fw) {
    const int mvl = SNOVA_m1 * nrows * SNOVA_l;
    memset(Fw, 0, ((size_t)mvl * SNOVA_o * SNOVA_l + 16) * sizeof(uint16_t));
    for (int ni = 0; ni < SNOVA_v; ++ni)
        for (int k1 = 0; k1 < SNOVA_l; ++k1)
            for (int mi = 0; mi < SNOVA_m1; ++mi)
                for (int nj = 0; nj < nrows; ++nj)
                    for (int i1 = 0; i1 < SNOVA_l; ++i1)
                        P11aw[(ni * SNOVA_l + k1) * mvl + (mi * nrows + nj) * SNOVA_l + i1] =
                            A[((mi * nrows + nj) * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + k1];
    for (int nk = 0; nk < SNOVA_o; ++nk)
        for (int j1 = 0; j1 < SNOVA_l; ++j1)
            for (int ni = 0; ni < SNOVA_v; ++ni)
                for (int k1 = 0; k1 < SNOVA_l; ++k1) {
                    uint16_t s = cl_expand_scalar16(T12[(ni * SNOVA_o + nk) * SNOVA_l2 + k1 * SNOVA_l + j1]);
                    uint16_t *Fr = &Fw[(nk * SNOVA_l + j1) * mvl];
                    const uint16_t *Pr = &P11aw[(ni * SNOVA_l + k1) * mvl];
                    for (int mi = 0; mi < mvl; ++mi) Fr[mi] ^= (uint16_t)(s * Pr[mi]);
                }
    for (int i = 0; i < mvl * SNOVA_o * SNOVA_l; i += 16)
        _mm256_storeu_si256((__m256i *)&Fw[i], cl_gf16_compress_u16x16(_mm256_loadu_si256((const __m256i *)&Fw[i])));
    for (int mi = 0; mi < SNOVA_m1; ++mi)
        for (int nj = 0; nj < nrows; ++nj)
            for (int nk = 0; nk < SNOVA_o; ++nk)
                for (int i1 = 0; i1 < SNOVA_l; ++i1)
                    for (int j1 = 0; j1 < SNOVA_l; ++j1) {
                        gf_t r = (gf_t)Fw[(nk * SNOVA_l + j1) * mvl + (mi * nrows + nj) * SNOVA_l + i1];
                        gf_t *c = &C[((mi * nrows + nj) * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + j1];
                        *c = rmw ? (gf_t)(*c ^ r) : r;
                    }
    SNOVA_CLEAR(Fw, ((size_t)mvl * SNOVA_o * SNOVA_l + 16) * sizeof(uint16_t));
}
static void rct_wide_fold_F12(gf_t *F12, const gf_t *P11, const gf_t *T12,
                              uint16_t *const P11aw, uint16_t *const Fw) {
    rct_wide_fold_bsec(F12, P11, T12, SNOVA_v, 0, P11aw, Fw);
}
#endif

#if RCT_F5_A4
_Static_assert(SNOVA_l == 5, "RCT_F5_A4 fold f-clone requires l == 5");
typedef struct { __m128i b012, b34; } rct_a4f_bc_t;
typedef struct { __m128i a, b; } rct_a4f_acc_t;
static inline rct_a4f_acc_t rct_a4f_zero(void) {
    rct_a4f_acc_t z; z.a = z.b = _mm_setzero_si128(); return z;
}
static inline rct_a4f_bc_t rct_a4f_bc(const gf_t *sc) {
    const __m128i P012 = _mm_setr_epi8(0,0,0,0,0, 1,1,1,1,1, 2,2,2,2,2, -1);
    const __m128i P34  = _mm_setr_epi8(3,3,3,3,3, 4,4,4,4,4, -1,-1,-1,-1,-1,-1);
    __m128i v = _mm_loadu_si128((const __m128i *)sc);
    rct_a4f_bc_t b;
    b.b012 = _mm_shuffle_epi8(v, P012);
    b.b34  = _mm_shuffle_epi8(v, P34);
    return b;
}
static inline void rct_a4f_mac(rct_a4f_acc_t *A, rct_a4f_bc_t b, const gf_t *w) {
    A->a = _mm_xor_si128(A->a, _mm_gf2p8mul_epi8(b.b012,
               _mm_loadu_si128((const __m128i *)w)));
    A->b = _mm_xor_si128(A->b, _mm_gf2p8mul_epi8(b.b34,
               _mm_loadu_si128((const __m128i *)&w[15])));
}
static inline __m128i rct_a4f_fold(rct_a4f_acc_t A) {
    __m128i f = _mm_xor_si128(A.a, _mm_srli_si128(A.a, 5));
    f = _mm_xor_si128(f, _mm_srli_si128(A.a, 10));
    f = _mm_xor_si128(f, A.b);
    return _mm_xor_si128(f, _mm_srli_si128(A.b, 5));
}
static inline void rct_a4f_xor5(gf_t *C, __m128i acc_folded) {
    _Alignas(16) uint8_t pb[16];
    _mm_store_si128((__m128i *)pb, RCT_SJ_CLEAN(acc_folded));
    for (int j = 0; j < 5; ++j) C[j] = (gf_t)(C[j] ^ pb[j]);
}
static inline void rct_a4f_cell_add(gf_t *C, const rct_a4f_bc_t bc[5], const gf_t *B) {
    for (int i = 0; i < 5; ++i) {
        rct_a4f_acc_t a = rct_a4f_zero();
        rct_a4f_mac(&a, bc[i], B);
        rct_a4f_xor5(&C[i * 5], rct_a4f_fold(a));
    }
}
#endif

#if RCT_F5_M4
_Static_assert(SNOVA_l == 5, "RCT_F5_M4 fold f-clone requires l == 5");
#include "gf16_core/gf16_mullo16.h"
#define RCT_M4F_PAT0 _mm256_setr_epi8(0,1,0,1,0,1,0,1,0,1,2,3,2,3,2,3,2,3,2,3,4,5,4,5,4,5,4,5,4,5,6,7)
#define RCT_M4F_PAT1 _mm256_setr_epi8(6,7,6,7,6,7,6,7,8,9,8,9,8,9,8,9,8,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
typedef struct { __m256i b0, b1; } rct_m4f_bc_t;
typedef struct { __m256i a0, a1; } rct_m4f_acc_t;
static inline rct_m4f_acc_t rct_m4f_zero(void) {
    rct_m4f_acc_t z; z.a0 = z.a1 = _mm256_setzero_si256(); return z;
}
static inline rct_m4f_bc_t rct_m4f_bc(const gf_t *sc) {
    _Alignas(16) uint16_t ev[8] = {0};
    for (int k = 0; k < 5; ++k) ev[k] = cl_expand_scalar16(sc[k]);
    __m256i evb = _mm256_broadcastsi128_si256(_mm_load_si128((const __m128i *)ev));
    rct_m4f_bc_t b;
    b.b0 = _mm256_shuffle_epi8(evb, RCT_M4F_PAT0);
    b.b1 = _mm256_shuffle_epi8(evb, RCT_M4F_PAT1);
    return b;
}
static inline void rct_m4f_mac(rct_m4f_acc_t *A, rct_m4f_bc_t b, const gf_t *w) {
    A->a0 = _mm256_xor_si256(A->a0, _mm256_mullo_epi16(b.b0,
        _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)w))));
    A->a1 = _mm256_xor_si256(A->a1, _mm256_mullo_epi16(b.b1,
        _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)&w[16]))));
}
static inline __m128i rct_m4f_fold(rct_m4f_acc_t A) {
    _Alignas(32) uint16_t tb[48];
    _mm256_store_si256((__m256i *)tb, A.a0);
    _mm256_store_si256((__m256i *)&tb[16], A.a1);
    _mm256_store_si256((__m256i *)&tb[32], _mm256_setzero_si256());
    __m128i f = _mm_loadu_si128((const __m128i *)tb);
    f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[5]));
    f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[10]));
    f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[15]));
    f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[20]));
    __m256i c = cl_gf16_compress_u16x16(_mm256_castsi128_si256(f));
    return cl_gf16_pack_u16_to_bytes(c);
}
static inline void rct_m4f_xor5(gf_t *C, __m128i acc_folded) {
    _Alignas(16) uint8_t pb[16];
    _mm_store_si128((__m128i *)pb, acc_folded);
    for (int j = 0; j < 5; ++j) C[j] = (gf_t)(C[j] ^ pb[j]);
}
static inline void rct_m4f_cell_add(gf_t *C, const rct_m4f_bc_t bc[5], const gf_t *B) {
    for (int i = 0; i < 5; ++i) {
        rct_m4f_acc_t a = rct_m4f_zero();
        rct_m4f_mac(&a, bc[i], B);
        rct_m4f_xor5(&C[i * 5], rct_m4f_fold(a));
    }
}

static inline rct_m4f_bc_t rct_m4f2_bc_mirror(const uint16_t *ev5) {
    __m256i evb = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)ev5));
    rct_m4f_bc_t b;
    b.b0 = _mm256_shuffle_epi8(evb, RCT_M4F_PAT0);
    b.b1 = _mm256_shuffle_epi8(evb, RCT_M4F_PAT1);
    return b;
}
static inline void rct_m4f2_build_bc(rct_m4f_bc_t (*bc)[5], const gf_t *src, int n_cells) {
    for (int c = 0; c < n_cells; ++c) {
        _Alignas(32) uint16_t ev[32];
        rct_m4_expand_buf(ev, &src[c * SNOVA_l2], SNOVA_l2);
        for (int i = 0; i < 5; ++i) bc[c][i] = rct_m4f2_bc_mirror(&ev[i * 5]);
    }
}
static inline void rct_m4f2_mac5(__m256i *a0, __m256i *a1, const rct_m4f_bc_t *bc, const gf_t *w) {
    const __m256i w0 = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)w));
    const __m256i w1 = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)&w[16]));
    for (int i = 0; i < 5; ++i) {
        a0[i] = _mm256_xor_si256(a0[i], _mm256_mullo_epi16(bc[i].b0, w0));
        a1[i] = _mm256_xor_si256(a1[i], _mm256_mullo_epi16(bc[i].b1, w1));
    }
}
static inline void rct_m4f2_fold5_xor(gf_t *C, const __m256i *a0, const __m256i *a1) {
    _Alignas(32) uint16_t tb[5][48];
    for (int i = 0; i < 5; ++i) {
        _mm256_store_si256((__m256i *)tb[i], a0[i]);
        _mm256_store_si256((__m256i *)&tb[i][16], a1[i]);
        _mm256_store_si256((__m256i *)&tb[i][32], _mm256_setzero_si256());
    }
    for (int i = 0; i < 5; ++i) {
        __m128i f = _mm_loadu_si128((const __m128i *)tb[i]);
        f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[i][5]));
        f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[i][10]));
        f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[i][15]));
        f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[i][20]));
        __m256i c = cl_gf16_compress_u16x16(_mm256_castsi128_si256(f));
        _Alignas(16) uint8_t pb[16];
        _mm_store_si128((__m128i *)pb, cl_gf16_pack_u16_to_bytes(c));
        for (int j = 0; j < 5; ++j) C[i * 5 + j] = (gf_t)(C[i * 5 + j] ^ pb[j]);
    }
}
static inline void rct_m4f2_fold_F12(gf_t *F12, const gf_t *P11, const gf_t *T12) {
    for (int i1 = 0; i1 < SNOVA_m1; ++i1)
        for (int j1 = 0; j1 < SNOVA_v; ++j1) {
            rct_m4f_bc_t bcp[SNOVA_v][5];
            rct_m4f2_build_bc(bcp, &P11[(i1 * SNOVA_v + j1) * SNOVA_v * SNOVA_l2], SNOVA_v);
            for (int k1 = 0; k1 < SNOVA_o; ++k1) {
                gf_t *C = &F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2];
                __m256i a0[5], a1[5];
                for (int i = 0; i < 5; ++i) { a0[i] = _mm256_setzero_si256(); a1[i] = _mm256_setzero_si256(); }
                for (int j2 = 0; j2 < SNOVA_v; ++j2)
                    rct_m4f2_mac5(a0, a1, bcp[j2], &T12[(j2 * SNOVA_o + k1) * SNOVA_l2]);
                rct_m4f2_fold5_xor(C, a0, a1);
            }
        }
}
static inline void rct_m4f2_fold_P22_pass1(gf_t *P22, const gf_t *T12, const gf_t *F12) {
    for (int j1 = 0; j1 < SNOVA_o; ++j1) {
        rct_m4f_bc_t bc[SNOVA_v][5];
        for (int idx = 0; idx < SNOVA_v; ++idx) {
            _Alignas(32) uint16_t ev[32];
            rct_m4_expand_buf(ev, &T12[(idx * SNOVA_o + j1) * SNOVA_l2], SNOVA_l2);
            for (int i = 0; i < 5; ++i) bc[idx][i] = rct_m4f2_bc_mirror(&ev[i * 5]);
        }
        for (int i1 = 0; i1 < SNOVA_m1; ++i1)
            for (int k1 = 0; k1 < SNOVA_o; ++k1) {
                gf_t *C = &P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2];
                __m256i a0[5], a1[5];
                for (int i = 0; i < 5; ++i) { a0[i] = _mm256_setzero_si256(); a1[i] = _mm256_setzero_si256(); }
                for (int idx = 0; idx < SNOVA_v; ++idx)
                    rct_m4f2_mac5(a0, a1, bc[idx], &F12[((i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]);
                rct_m4f2_fold5_xor(C, a0, a1);
            }
        SNOVA_CLEAR_OBJ(bc);
    }
}
static inline void rct_m4f2_fold_P22_pass2(gf_t *P22, const gf_t *P21, const gf_t *T12) {
    for (int i1 = 0; i1 < SNOVA_m1; ++i1)
        for (int j1 = 0; j1 < SNOVA_o; ++j1) {
            rct_m4f_bc_t bcp[SNOVA_v][5];
            rct_m4f2_build_bc(bcp, &P21[(i1 * SNOVA_o + j1) * SNOVA_v * SNOVA_l2], SNOVA_v);
            for (int k1 = 0; k1 < SNOVA_o; ++k1) {
                gf_t *C = &P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2];
                __m256i a0[5], a1[5];
                for (int i = 0; i < 5; ++i) { a0[i] = _mm256_setzero_si256(); a1[i] = _mm256_setzero_si256(); }
                for (int idx = 0; idx < SNOVA_v; ++idx)
                    rct_m4f2_mac5(a0, a1, bcp[idx], &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
                rct_m4f2_fold5_xor(C, a0, a1);
            }
        }
}
static inline void rct_m4f2_fold_row_bsec(gf_t *Crow, const gf_t *Arow, const gf_t *T12) {
    rct_m4f_bc_t bcp[SNOVA_v][5];
    rct_m4f2_build_bc(bcp, Arow, SNOVA_v);
    for (int k1 = 0; k1 < SNOVA_o; ++k1) {
        gf_t *C = &Crow[k1 * SNOVA_l2];
        __m256i a0[5], a1[5];
        for (int i = 0; i < 5; ++i) { a0[i] = _mm256_setzero_si256(); a1[i] = _mm256_setzero_si256(); }
        for (int idx = 0; idx < SNOVA_v; ++idx)
            rct_m4f2_mac5(a0, a1, bcp[idx], &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
        rct_m4f2_fold5_xor(C, a0, a1);
    }
}
#endif

#if RCT_F5_A4 || RCT_F5_M4 || RCT_F5_WIDE
#if RCT_F5_A4 || RCT_F5_M4
#if RCT_F5_A4
#define RCT_F5_BC_T     rct_a4f_bc_t
#define RCT_F5_ACC_T    rct_a4f_acc_t
#define RCT_F5_ZERO     rct_a4f_zero
#define RCT_F5_BC       rct_a4f_bc
#define RCT_F5_MAC      rct_a4f_mac
#define RCT_F5_FOLD     rct_a4f_fold
#define RCT_F5_XOR5     rct_a4f_xor5
#define RCT_F5_CELL_ADD rct_a4f_cell_add
#else
#define RCT_F5_BC_T     rct_m4f_bc_t
#define RCT_F5_ACC_T    rct_m4f_acc_t
#define RCT_F5_ZERO     rct_m4f_zero
#define RCT_F5_BC       rct_m4f_bc
#define RCT_F5_MAC      rct_m4f_mac
#define RCT_F5_FOLD     rct_m4f_fold
#define RCT_F5_XOR5     rct_m4f_xor5
#define RCT_F5_CELL_ADD rct_m4f_cell_add
#endif
#endif

static inline void rct_f5_fold_F12(gf_t *F12, const gf_t *P11, const gf_t *T12) {
#if RCT_F5_WIDE
    { RCT_WIDE_SCRATCH_DECL; rct_wide_fold_F12(F12, P11, T12, P11aw, Fw); }
#elif RCT_F5_M4
    rct_m4f2_fold_F12(F12, P11, T12);
#else
    for (int i1 = 0; i1 < SNOVA_m1; ++i1)
        for (int j1 = 0; j1 < SNOVA_v; ++j1) {
            RCT_F5_BC_T bcp[SNOVA_v][5];
            for (int j2 = 0; j2 < SNOVA_v; ++j2)
                for (int i = 0; i < 5; ++i)
                    bcp[j2][i] = RCT_F5_BC(&P11[((i1 * SNOVA_v + j1) * SNOVA_v + j2) * SNOVA_l2 + i * 5]);
            for (int k1 = 0; k1 < SNOVA_o; ++k1) {
                gf_t *C = &F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2];
                RCT_F5_ACC_T a[5];
                for (int i = 0; i < 5; ++i) a[i] = RCT_F5_ZERO();
                for (int j2 = 0; j2 < SNOVA_v; ++j2) {
                    const gf_t *w = &T12[(j2 * SNOVA_o + k1) * SNOVA_l2];
                    for (int i = 0; i < 5; ++i) RCT_F5_MAC(&a[i], bcp[j2][i], w);
                }
                for (int i = 0; i < 5; ++i) RCT_F5_XOR5(&C[i * 5], RCT_F5_FOLD(a[i]));
            }
        }
#endif
}

#if RCT_F5_A4 || RCT_F5_M4
static inline void rct_f5_fold_P22_pass1(gf_t *P22, const gf_t *T12, const gf_t *F12) {
#if RCT_F5_M4
    rct_m4f2_fold_P22_pass1(P22, T12, F12);
#else
    for (int j1 = 0; j1 < SNOVA_o; ++j1) {
        RCT_F5_BC_T bc[SNOVA_v][5];
        for (int idx = 0; idx < SNOVA_v; ++idx)
            for (int i = 0; i < 5; ++i)
                bc[idx][i] = RCT_F5_BC(&T12[(idx * SNOVA_o + j1) * SNOVA_l2 + i * 5]);
        for (int i1 = 0; i1 < SNOVA_m1; ++i1)
            for (int k1 = 0; k1 < SNOVA_o; ++k1) {
                gf_t *C = &P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2];
                RCT_F5_ACC_T a[5];
                for (int i = 0; i < 5; ++i) a[i] = RCT_F5_ZERO();
                for (int idx = 0; idx < SNOVA_v; ++idx) {
                    const gf_t *w = &F12[((i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2];
                    for (int i = 0; i < 5; ++i) RCT_F5_MAC(&a[i], bc[idx][i], w);
                }
                for (int i = 0; i < 5; ++i) RCT_F5_XOR5(&C[i * 5], RCT_F5_FOLD(a[i]));
            }
        SNOVA_CLEAR_OBJ(bc);
    }
#endif
}

static inline void rct_f5_fold_P22_pass2(gf_t *P22, const gf_t *P21, const gf_t *T12) {
#if RCT_F5_WIDE
    { RCT_WIDE_SCRATCH_DECL; rct_wide_fold_bsec(P22, P21, T12, SNOVA_o, 1, P11aw, Fw); }
#elif RCT_F5_M4
    rct_m4f2_fold_P22_pass2(P22, P21, T12);
#else
    for (int i1 = 0; i1 < SNOVA_m1; ++i1)
        for (int j1 = 0; j1 < SNOVA_o; ++j1) {
            RCT_F5_BC_T bcp[SNOVA_v][5];
            for (int idx = 0; idx < SNOVA_v; ++idx)
                for (int i = 0; i < 5; ++i)
                    bcp[idx][i] = RCT_F5_BC(&P21[((i1 * SNOVA_o + j1) * SNOVA_v + idx) * SNOVA_l2 + i * 5]);
            for (int k1 = 0; k1 < SNOVA_o; ++k1) {
                gf_t *C = &P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2];
                RCT_F5_ACC_T a[5];
                for (int i = 0; i < 5; ++i) a[i] = RCT_F5_ZERO();
                for (int idx = 0; idx < SNOVA_v; ++idx) {
                    const gf_t *w = &T12[(idx * SNOVA_o + k1) * SNOVA_l2];
                    for (int i = 0; i < 5; ++i) RCT_F5_MAC(&a[i], bcp[idx][i], w);
                }
                for (int i = 0; i < 5; ++i) RCT_F5_XOR5(&C[i * 5], RCT_F5_FOLD(a[i]));
            }
        }
#endif
}

static inline void rct_f5_fold_row_bsec(gf_t *Crow, const gf_t *Arow, const gf_t *T12) {
#if RCT_F5_M4
    rct_m4f2_fold_row_bsec(Crow, Arow, T12);
#else
    RCT_F5_BC_T bcp[SNOVA_v][5];
    for (int idx = 0; idx < SNOVA_v; ++idx)
        for (int i = 0; i < 5; ++i)
            bcp[idx][i] = RCT_F5_BC(&Arow[idx * SNOVA_l2 + i * 5]);
    for (int k1 = 0; k1 < SNOVA_o; ++k1) {
        gf_t *C = &Crow[k1 * SNOVA_l2];
        RCT_F5_ACC_T a[5];
        for (int i = 0; i < 5; ++i) a[i] = RCT_F5_ZERO();
        for (int idx = 0; idx < SNOVA_v; ++idx) {
            const gf_t *w = &T12[(idx * SNOVA_o + k1) * SNOVA_l2];
            for (int i = 0; i < 5; ++i) RCT_F5_MAC(&a[i], bcp[idx][i], w);
        }
        for (int i = 0; i < 5; ++i) RCT_F5_XOR5(&C[i * 5], RCT_F5_FOLD(a[i]));
    }
#endif
}
#endif
#endif
#endif
#endif
