#ifndef RCT_SIGN_ENGINE_H
#define RCT_SIGN_ENGINE_H
#if RCT_SIGN_JOG
static inline __m128i rct_sj_cleanup128(__m128i v) {
    const __m128i m0f = _mm_set1_epi8(0x0f);
    __m128i vhi = _mm_and_si128(v, _mm_set1_epi8((char)0xf0));
    __m128i a = _mm_srli_epi16(vhi, 3);
    __m128i b = _mm_srli_epi16(v, 4);
    return _mm_and_si128(_mm_xor_si128(_mm_xor_si128(v, a), b), m0f);
}
static inline __m256i rct_sj_cleanup256(__m256i v) {
    const __m256i m0f = _mm256_set1_epi8(0x0f);
    __m256i vhi = _mm256_and_si256(v, _mm256_set1_epi8((char)0xf0));
    __m256i a = _mm256_srli_epi16(vhi, 3);
    __m256i b = _mm256_srli_epi16(v, 4);
    return _mm256_and_si256(_mm256_xor_si256(_mm256_xor_si256(v, a), b), m0f);
}

#if RCT_HAVE_GFNI
static inline __m128i rct_sj_bc128(gf_t s) { return _mm_set1_epi8((char)s); }
static inline __m128i rct_sj_sv128(__m128i bc, __m128i v) { return _mm_gf2p8mul_epi8(bc, v); }
#define RCT_SJ_CLEAN(v) rct_sj_cleanup128(v)
static inline __m256i rct_sj_bc256(gf_t s) { return _mm256_set1_epi8((char)s); }
static inline __m256i rct_sj_sv256(__m256i bc, __m256i v) { return _mm256_gf2p8mul_epi8(bc, v); }
#define RCT_SJ_CLEAN256(v) rct_sj_cleanup256(v)
static inline void rct_sj_ensure(void) {}
#else
static __m128i rct_sj_vtl[16];
static int rct_sj_vtl_done = 0;
static inline void rct_sj_ensure(void) {
    if (rct_sj_vtl_done) return;
    for (int k = 0; k < 16; ++k) {
        _Alignas(16) uint8_t t[16];
        for (int x = 0; x < 16; ++x) t[x] = gf_mult((gf_t)k, (gf_t)x);
        rct_sj_vtl[k] = _mm_load_si128((const __m128i *)t);
    }
    rct_sj_vtl_done = 1;
}
static inline __m128i rct_sj_bc128(gf_t s) {
    const __m128i sv = _mm_set1_epi8((char)(s & 0x0F));
    const __m128i b1 = _mm_set1_epi8(1), b2 = _mm_set1_epi8(2),
                  b4 = _mm_set1_epi8(4), b8 = _mm_set1_epi8(8);
    __m128i r = _mm_and_si128(rct_sj_vtl[1], _mm_cmpeq_epi8(_mm_and_si128(sv, b1), b1));
    r = _mm_xor_si128(r, _mm_and_si128(rct_sj_vtl[2], _mm_cmpeq_epi8(_mm_and_si128(sv, b2), b2)));
    r = _mm_xor_si128(r, _mm_and_si128(rct_sj_vtl[4], _mm_cmpeq_epi8(_mm_and_si128(sv, b4), b4)));
    r = _mm_xor_si128(r, _mm_and_si128(rct_sj_vtl[8], _mm_cmpeq_epi8(_mm_and_si128(sv, b8), b8)));
    return r;
}
static inline __m128i rct_sj_sv128(__m128i bc, __m128i v) { return _mm_shuffle_epi8(bc, v); }
#define RCT_SJ_CLEAN(v) (v)
static inline __m256i rct_sj_bc256(gf_t s) { __m128i b = rct_sj_bc128(s); return _mm256_set_m128i(b, b); }
static inline __m256i rct_sj_sv256(__m256i bc, __m256i v) { return _mm256_shuffle_epi8(bc, v); }
#define RCT_SJ_CLEAN256(v) (v)
#endif

static inline __m128i rct_sj_bc128_pub(gf_t s) {
    SNOVA_CT_ASSERT_PUBLIC_MEM(&s, sizeof s);
#if RCT_HAVE_GFNI
    return _mm_set1_epi8((char)s);
#else
    return rct_sj_vtl[s & 0x0F];
#endif
}

static inline __m256i rct_sj_bc256_pub(gf_t s) {
    SNOVA_CT_ASSERT_PUBLIC_MEM(&s, sizeof s);
#if RCT_HAVE_GFNI
    return _mm256_set1_epi8((char)s);
#else
    return _mm256_broadcastsi128_si256(rct_sj_vtl[s & 0x0F]);
#endif
}

static inline void rct_sj_store_r(gf_t *dst, __m128i acc) {
    _Alignas(16) uint8_t tmp[16];
    _mm_store_si128((__m128i *)tmp, RCT_SJ_CLEAN(acc));
    memcpy(dst, tmp, SNOVA_r);
}

static inline void rct_sj_mm_add(gf_t *C, const gf_t *A, const gf_t *B, int ad, int bd, int cd) {
    for (int i = 0; i < ad; ++i) {
        __m128i acc = _mm_setzero_si128();
        for (int k = 0; k < bd; ++k)
            acc = _mm_xor_si128(acc, rct_sj_sv128(rct_sj_bc128(A[i * bd + k]),
                _mm_loadu_si128((const __m128i *)&B[k * cd])));
        _Alignas(16) uint8_t pb[16];
        _mm_store_si128((__m128i *)pb, RCT_SJ_CLEAN(acc));
        for (int j = 0; j < cd; ++j) C[i * cd + j] = (gf_t)(C[i * cd + j] ^ pb[j]);
    }
}

static inline void rct_sj_mm_add_pub(gf_t *C, const gf_t *A, const gf_t *B, int ad, int bd, int cd) {
    for (int i = 0; i < ad; ++i) {
        __m128i acc = _mm_setzero_si128();
        for (int k = 0; k < bd; ++k)
            acc = _mm_xor_si128(acc, rct_sj_sv128(rct_sj_bc128_pub(A[i * bd + k]),
                _mm_loadu_si128((const __m128i *)&B[k * cd])));
        _Alignas(16) uint8_t pb[16];
        _mm_store_si128((__m128i *)pb, RCT_SJ_CLEAN(acc));
        for (int j = 0; j < cd; ++j) C[i * cd + j] = (gf_t)(C[i * cd + j] ^ pb[j]);
    }
}

#if RCT_HAVE_GFNI && SNOVA_l == 5 && (SNOVA_r == 5 || SNOVA_r == 6 || SNOVA_r == 8)
#define RCT_SJ_A4 1
#if SNOVA_r == 8
typedef struct { __m256i bv; __m128i bt; } rct_a4_bc_t;
typedef struct { __m256i a; __m128i t; } rct_a4_acc_t;
static inline rct_a4_acc_t rct_a4_zero(void) {
    rct_a4_acc_t z; z.a = _mm256_setzero_si256(); z.t = _mm_setzero_si128(); return z;
}
static inline rct_a4_bc_t rct_a4_bc(const gf_t *sc) {
    const __m256i PAT = _mm256_setr_epi8(0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1,
                                         2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3);
    rct_a4_bc_t b;
    b.bv = _mm256_shuffle_epi8(_mm256_broadcastsi128_si256(
               _mm_loadu_si128((const __m128i *)sc)), PAT);
    b.bt = _mm_set1_epi8((char)sc[4]);
    return b;
}
static inline void rct_a4_mac(rct_a4_acc_t *A, rct_a4_bc_t b, const gf_t *w) {
    A->a = _mm256_xor_si256(A->a, _mm256_gf2p8mul_epi8(b.bv,
               _mm256_loadu_si256((const __m256i *)w)));
    A->t = _mm_xor_si128(A->t, _mm_gf2p8mul_epi8(b.bt,
               _mm_loadu_si128((const __m128i *)&w[32])));
}
static inline __m128i rct_a4_fold(rct_a4_acc_t A) {
    __m128i f = _mm_xor_si128(_mm256_castsi256_si128(A.a),
                              _mm256_extracti128_si256(A.a, 1));
    f = _mm_xor_si128(f, _mm_unpackhi_epi64(f, f));
    return _mm_xor_si128(f, A.t);
}
#elif SNOVA_r == 6
typedef struct { __m128i b01, b23, bt; } rct_a4_bc_t;
typedef struct { __m128i a01, a23, t; } rct_a4_acc_t;
static inline rct_a4_acc_t rct_a4_zero(void) {
    rct_a4_acc_t z; z.a01 = z.a23 = z.t = _mm_setzero_si128(); return z;
}
static inline rct_a4_bc_t rct_a4_bc(const gf_t *sc) {
    const __m128i P01 = _mm_setr_epi8(0,0,0,0,0,0, 1,1,1,1,1,1, -1,-1,-1,-1);
    const __m128i P23 = _mm_setr_epi8(2,2,2,2,2,2, 3,3,3,3,3,3, -1,-1,-1,-1);
    __m128i v = _mm_loadu_si128((const __m128i *)sc);
    rct_a4_bc_t b;
    b.b01 = _mm_shuffle_epi8(v, P01);
    b.b23 = _mm_shuffle_epi8(v, P23);
    b.bt  = _mm_set1_epi8((char)sc[4]);
    return b;
}
static inline void rct_a4_mac(rct_a4_acc_t *A, rct_a4_bc_t b, const gf_t *w) {
    A->a01 = _mm_xor_si128(A->a01, _mm_gf2p8mul_epi8(b.b01,
                 _mm_loadu_si128((const __m128i *)w)));
    A->a23 = _mm_xor_si128(A->a23, _mm_gf2p8mul_epi8(b.b23,
                 _mm_loadu_si128((const __m128i *)&w[12])));
    A->t   = _mm_xor_si128(A->t, _mm_gf2p8mul_epi8(b.bt,
                 _mm_loadu_si128((const __m128i *)&w[24])));
}
static inline __m128i rct_a4_fold(rct_a4_acc_t A) {
    __m128i f = _mm_xor_si128(A.a01, A.a23);
    f = _mm_xor_si128(f, _mm_srli_si128(f, 6));
    return _mm_xor_si128(f, A.t);
}
#else
typedef struct { __m128i b012, b34; } rct_a4_bc_t;
typedef struct { __m128i a, b; } rct_a4_acc_t;
static inline rct_a4_acc_t rct_a4_zero(void) {
    rct_a4_acc_t z; z.a = z.b = _mm_setzero_si128(); return z;
}
static inline rct_a4_bc_t rct_a4_bc(const gf_t *sc) {
    const __m128i P012 = _mm_setr_epi8(0,0,0,0,0, 1,1,1,1,1, 2,2,2,2,2, -1);
    const __m128i P34  = _mm_setr_epi8(3,3,3,3,3, 4,4,4,4,4, -1,-1,-1,-1,-1,-1);
    __m128i v = _mm_loadu_si128((const __m128i *)sc);
    rct_a4_bc_t b;
    b.b012 = _mm_shuffle_epi8(v, P012);
    b.b34  = _mm_shuffle_epi8(v, P34);
    return b;
}
static inline void rct_a4_mac(rct_a4_acc_t *A, rct_a4_bc_t b, const gf_t *w) {
    A->a = _mm_xor_si128(A->a, _mm_gf2p8mul_epi8(b.b012,
               _mm_loadu_si128((const __m128i *)w)));
    A->b = _mm_xor_si128(A->b, _mm_gf2p8mul_epi8(b.b34,
               _mm_loadu_si128((const __m128i *)&w[15])));
}
static inline __m128i rct_a4_fold(rct_a4_acc_t A) {
    __m128i f = _mm_xor_si128(A.a, _mm_srli_si128(A.a, 5));
    f = _mm_xor_si128(f, _mm_srli_si128(A.a, 10));
    f = _mm_xor_si128(f, A.b);
    return _mm_xor_si128(f, _mm_srli_si128(A.b, 5));
}
#endif
#else
#define RCT_SJ_A4 0
#endif

#if RCT_SIGN_JOG && !RCT_HAVE_GFNI && SNOVA_l == 5 && (SNOVA_r == 5 || SNOVA_r == 6 || SNOVA_r == 8)
#define RCT_SJ_M4 1
#include "gf16_core/gf16_mullo16.h"
#if SNOVA_r == 5
#define RCT_M4_PAT0 _mm256_setr_epi8(0,1,0,1,0,1,0,1,0,1,2,3,2,3,2,3,2,3,2,3,4,5,4,5,4,5,4,5,4,5,6,7)
#define RCT_M4_PAT1 _mm256_setr_epi8(6,7,6,7,6,7,6,7,8,9,8,9,8,9,8,9,8,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
#endif
#if SNOVA_r == 6
#define RCT_M4_PAT0 _mm256_setr_epi8(0,1,0,1,0,1,0,1,0,1,0,1,2,3,2,3,2,3,2,3,2,3,2,3,4,5,4,5,4,5,4,5)
#define RCT_M4_PAT1 _mm256_setr_epi8(4,5,4,5,6,7,6,7,6,7,6,7,6,7,6,7,8,9,8,9,8,9,8,9,8,9,8,9,-1,-1,-1,-1)
#endif
#if SNOVA_r == 8
#define RCT_M4_PAT0 _mm256_setr_epi8(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3)
#define RCT_M4_PAT1 _mm256_setr_epi8(4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,6,7,6,7,6,7,6,7,6,7,6,7,6,7,6,7)
#endif

typedef struct { __m256i b0, b1;
#if SNOVA_r == 8
    __m256i b2;
#endif
} rct_m4_bc_t;
typedef struct { __m256i a0, a1;
#if SNOVA_r == 8
    __m256i a2;
#endif
} rct_m4_acc_t;
static inline rct_m4_acc_t rct_m4_zero(void) {
    rct_m4_acc_t z; z.a0 = z.a1 = _mm256_setzero_si256();
#if SNOVA_r == 8
    z.a2 = _mm256_setzero_si256();
#endif
    return z;
}
static inline rct_m4_bc_t rct_m4_bc(const gf_t *sc) {
    _Alignas(16) uint16_t ev[8];
    for (int k = 0; k < 5; ++k) ev[k] = cl_expand_scalar16(sc[k]);
    __m256i evb = _mm256_broadcastsi128_si256(_mm_load_si128((const __m128i *)ev));
    rct_m4_bc_t b;
    b.b0 = _mm256_shuffle_epi8(evb, RCT_M4_PAT0);
    b.b1 = _mm256_shuffle_epi8(evb, RCT_M4_PAT1);
#if SNOVA_r == 8
    b.b2 = _mm256_set1_epi16((short)ev[4]);
#endif
    return b;
}
static inline void rct_m4_mac(rct_m4_acc_t *A, rct_m4_bc_t b, const gf_t *w) {
    A->a0 = _mm256_xor_si256(A->a0, _mm256_mullo_epi16(b.b0,
        _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)w))));
    A->a1 = _mm256_xor_si256(A->a1, _mm256_mullo_epi16(b.b1,
        _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)&w[16]))));
#if SNOVA_r == 8
    A->a2 = _mm256_xor_si256(A->a2, _mm256_mullo_epi16(b.b2,
        _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)&w[32]))));
#endif
}
static inline __m128i rct_m4_fold(rct_m4_acc_t A) {
    _Alignas(32) uint16_t tb[48];
    _mm256_store_si256((__m256i *)tb, A.a0);
    _mm256_store_si256((__m256i *)&tb[16], A.a1);
#if SNOVA_r == 8
    _mm256_store_si256((__m256i *)&tb[32], A.a2);
#else
    _mm256_store_si256((__m256i *)&tb[32], _mm256_setzero_si256());
#endif
    __m128i f = _mm_loadu_si128((const __m128i *)tb);
    f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[SNOVA_r]));
    f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[2 * SNOVA_r]));
    f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[3 * SNOVA_r]));
    f = _mm_xor_si128(f, _mm_loadu_si128((const __m128i *)&tb[4 * SNOVA_r]));
    __m256i c = cl_gf16_compress_u16x16(_mm256_castsi128_si256(f));
    return cl_gf16_pack_u16_to_bytes(c);
}
#else
#define RCT_SJ_M4 0
#endif

#if RCT_SIGN_JOG && !RCT_HAVE_GFNI && SNOVA_l == 5
#include "gf16_core/gf16_mullo16.h"
static inline void rct_m4_expand_buf(uint16_t *dst, const gf_t *src, int n) {
    for (int i = 0; i < n; i += 16)
        _mm256_storeu_si256((__m256i *)&dst[i], cl_gf16_expand_u16x16(
            _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)&src[i]))));
}
static inline void rct_m4_raw_buf(uint16_t *dst, const gf_t *src, int n) {
    for (int i = 0; i < n; i += 16)
        _mm256_storeu_si256((__m256i *)&dst[i],
            _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)&src[i])));
}
static inline __m128i rct_m4_compress128(__m128i a) {
    const __m128i m0f = _mm_set1_epi16(0x000f);
    __m128i val = _mm_xor_si128(
        _mm_xor_si128(_mm_and_si128(a, m0f),
                      _mm_srli_epi16(_mm_and_si128(a, _mm_set1_epi16(0x00f0)), 3)),
        _mm_xor_si128(_mm_srli_epi16(_mm_and_si128(a, _mm_set1_epi16(0x0f00)), 6),
                      _mm_srli_epi16(_mm_and_si128(a, _mm_set1_epi16((short)0xf000)), 9)));
    val = _mm_xor_si128(_mm_xor_si128(val,
              _mm_srli_epi16(_mm_and_si128(val, _mm_set1_epi16(0x00f0)), 3)),
              _mm_srli_epi16(val, 4));
    return _mm_and_si128(val, m0f);
}
#endif

#if RCT_SJ_M4
static inline void rct_m4_mm_add(gf_t *C, const gf_t *A, const gf_t *B, int ad, int bd, int cd) {
    _Alignas(32) uint16_t xe[64];
    rct_m4_expand_buf(xe, A, ad * bd);
    for (int i = 0; i < ad; ++i) {
        __m128i acc = _mm_setzero_si128();
        for (int k = 0; k < bd; ++k)
            acc = _mm_xor_si128(acc, _mm_mullo_epi16(_mm_set1_epi16((short)xe[i * bd + k]),
                _mm256_castsi256_si128(_mm256_cvtepu8_epi16(
                    _mm_loadu_si128((const __m128i *)&B[k * cd])))));
        _Alignas(16) uint8_t pb[16];
        _mm_store_si128((__m128i *)pb, _mm_packus_epi16(rct_m4_compress128(acc), _mm_setzero_si128()));
        for (int j = 0; j < cd; ++j) C[i * cd + j] = (gf_t)(C[i * cd + j] ^ pb[j]);
    }
}
#endif

_Static_assert(SNOVA_r <= 16, "RCT_SIGN_JOG engine requires r <= 16 (16-byte slot/lane)");
_Static_assert(SNOVA_l <= 16, "RCT_SIGN_JOG engine requires l <= 16 (16-byte slot/lane)");

static inline void rct_gauss_row_scale_sec(gf_t *row, gf_t s, int k0, int kend) {
    __m256i bs = rct_sj_bc256(s);
    for (int k = k0; k < kend; k += 32)
        _mm256_storeu_si256((__m256i *)&row[k],
            RCT_SJ_CLEAN256(rct_sj_sv256(bs, _mm256_loadu_si256((const __m256i *)&row[k]))));
}
static inline void rct_gauss_row_axpy_sec(gf_t *dst, const gf_t *src, gf_t s, int k0, int kend) {
    __m256i bs = rct_sj_bc256(s);
    for (int k = k0; k < kend; k += 32)
        _mm256_storeu_si256((__m256i *)&dst[k], _mm256_xor_si256(
            _mm256_loadu_si256((const __m256i *)&dst[k]),
            RCT_SJ_CLEAN256(rct_sj_sv256(bs, _mm256_loadu_si256((const __m256i *)&src[k])))));
}
#endif
#endif
