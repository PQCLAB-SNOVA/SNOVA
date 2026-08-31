/**
 * @file snova_rect.h
 */
#ifndef SNOVA_RECT_H
#define SNOVA_RECT_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "../snova_params.h"
#include "../primitives/sym_shake.h"
#include "../primitives/sym_aes.h"
#include "drbg.h"
#include "ct_poison.h"
#include "secure_clear.h"

#if !FIXED_ABQ
#error "rct core requires FIXED_ABQ=1 (R3): expand_public never writes the ABQ tail of P_matrix"
#endif

#if SNOVA_WRAPPER_STACK
#define RCT_SCRATCH
#else
#define RCT_SCRATCH static
#endif

typedef uint8_t gf_t;

#define SNOVA_l   SNOVA_L
#define SNOVA_r   SNOVA_R
#define SNOVA_v   SNOVA_V
#define SNOVA_o   SNOVA_O
#define SNOVA_n   SNOVA_N
#define SNOVA_q   SNOVA_Q
#define SNOVA_l2  (SNOVA_L * SNOVA_L)
#define SNOVA_r2  (SNOVA_R * SNOVA_R)
#define SNOVA_lr  (SNOVA_L * SNOVA_R)
#define SNOVA_m1  SNOVA_M1
#define SNOVA_alpha SNOVA_ALPHA

#define SEED_LENGTH_PUBLIC  16
#define SEED_LENGTH_PRIVATE 32
#define BYTES_SALT   16
#define BYTES_DIGEST 64
#define BYTES_PK_HASH SNOVA_BYTES_PK_HASH

#define PACK_GF    SNOVA_PACK_GF
#define PACK_BYTES SNOVA_PACK_BYTES
#define BYTES_GF(x) ((PACK_BYTES * (x) + PACK_GF - 1) / PACK_GF)

#define NUMGF_PK        (SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2)
#define NUMGF_SIGNATURE (SNOVA_n * SNOVA_lr)
#define GF16_HASH       (SNOVA_o * SNOVA_l * SNOVA_r)
#define BYTES_HASH      (BYTES_GF(GF16_HASH))
#define BYTES_SIGNATURE (BYTES_GF(NUMGF_SIGNATURE) + BYTES_SALT)
#define BYTES_PK        (BYTES_GF(NUMGF_PK) + SEED_LENGTH_PUBLIC)
#define BYTES_SK        (SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE + BYTES_PK_HASH)
#define CRYPTO_BYTES_R  (BYTES_SIGNATURE)

#define NUM_GEN_PUB_GF   (SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2)
#if SNOVA_q != 16
#define NUM_GEN_PUB_BYTES (NUM_GEN_PUB_GF)
#else
#define NUM_GEN_PUB_BYTES ((NUM_GEN_PUB_GF + 1) / 2)
#endif
#define NUM_PUB_GF                                                             \
    (SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2 +       \
     SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr) + 2 * SNOVA_o * SNOVA_alpha * SNOVA_l)
#define NUM_GEN_SEC_BYTES (BYTES_GF(SNOVA_v * SNOVA_lr))

#define i_prime(mi, alpha) (((alpha) + (mi)) % SNOVA_m1)

static gf_t rct_multtab[SNOVA_q * SNOVA_q];
static gf_t rct_addtab[SNOVA_q * SNOVA_q];
static gf_t rct_S[SNOVA_l * SNOVA_l2];

static inline gf_t gf_mult(const gf_t a, const gf_t b) { return rct_multtab[a * SNOVA_q + b]; }

static inline gf_t gf_add(const gf_t a, const gf_t b) {
#if SNOVA_q != 16
    unsigned t = (unsigned)a + b;
    return (gf_t)(t - (SNOVA_q & (0u - (unsigned)(t >= SNOVA_q))));
#else
    return (gf_t)(a ^ b);
#endif
}
static inline void gf_set_add(gf_t *a, const gf_t b) { *a = gf_add(*a, b); }
static inline gf_t gf_sub(const gf_t a, const gf_t b) {
#if SNOVA_q != 16
    unsigned t = (unsigned)a + SNOVA_q - b;
    return (gf_t)(t - (SNOVA_q & (0u - (unsigned)(t >= SNOVA_q))));
#else
    return (gf_t)(a ^ b);
#endif
}

#if SNOVA_q == 16
#include "../gf16_core/xgf16.h"
#endif
static inline gf_t gf_mult_sec(const gf_t a, const gf_t b) {
#if SNOVA_q != 16
    return (gf_t)(((unsigned)a * b) % SNOVA_q);
#else
    return (gf_t)xgf16_unspread(xgf16_reduce(xgf16_spread(a) * xgf16_spread(b)));
#endif
}
static inline gf_t gf_inv_sec(const gf_t a) {
#if SNOVA_q == 16
    uint32_t a0 = a & 1u, a1 = (a >> 1) & 1u, a2 = (a >> 2) & 1u, a3 = (a >> 3) & 1u;
    gf_t s2  = (gf_t)((a0 ^ a2) | (a2 << 1) | ((a1 ^ a3) << 2) | (a3 << 3));
    gf_t s3  = gf_mult_sec(s2, a);
    uint32_t b0 = s3 & 1u, b1 = (s3 >> 1) & 1u, b2 = (s3 >> 2) & 1u, b3 = (s3 >> 3) & 1u;
    gf_t s6  = (gf_t)((b0 ^ b2) | (b2 << 1) | ((b1 ^ b3) << 2) | (b3 << 3));
    uint32_t c0 = s6 & 1u, c1 = (s6 >> 1) & 1u, c2 = (s6 >> 2) & 1u, c3 = (s6 >> 3) & 1u;
    gf_t s12 = (gf_t)((c0 ^ c2) | (c2 << 1) | ((c1 ^ c3) << 2) | (c3 << 3));
    return gf_mult_sec(s12, s2);
#else
    gf_t val = a;
    for (int j1 = 3; j1 < SNOVA_q; j1++) val = gf_mult_sec(val, a);
    return val;
#endif
}
static inline uint32_t ct_gf_nz(const uint32_t v) { return (0u - v) >> 31; }
static inline gf_t ct_gf_sel(const uint32_t cond, const gf_t a, const gf_t b) {
    uint32_t m = 0u - cond;
    return (gf_t)((a & m) | (b & ~m));
}

#if defined(__AVX2__) && (SNOVA_Q == 16) && !defined(RCT_FORCE_SCALAR)
#define RCT_AVX2_Q16 1
#include <immintrin.h>
#include "gf16_core/gf16_qrp16.h"
#define SNOVA_lr16 (((SNOVA_lr) + 31) / 32)
#define SNOVA_lr32 (SNOVA_lr16 * 32)
#if defined(__GFNI__) && !defined(RCT_FORCE_PSHUFB)
#define RCT_HAVE_GFNI 1
#else
#define RCT_HAVE_GFNI 0
#endif
#else
#define RCT_AVX2_Q16 0
#define RCT_HAVE_GFNI 0
#endif

#if RCT_AVX2_Q16 && (SNOVA_l == 4)
#if RCT_HAVE_GFNI
#define RCT_USE_GFNI 1
#define RCT_USE_PSHUFB 0
#else
#define RCT_USE_GFNI 0
#define RCT_USE_PSHUFB 1
#endif
#else
#define RCT_USE_GFNI 0
#define RCT_USE_PSHUFB 0
#endif
#define RCT_USE_SIMD (RCT_USE_GFNI || RCT_USE_PSHUFB)

#if RCT_AVX2_Q16 && ((SNOVA_l != 4) || defined(RCT_FORCE_JOG))
#define RCT_VF_JOG 1
#else
#define RCT_VF_JOG 0
#endif

#if RCT_AVX2_Q16 && (SNOVA_l != 4)
#define RCT_SIGN_JOG 1
#else
#define RCT_SIGN_JOG 0
#endif

#include "platforms/x86_avx2/rct_sign_engine.h"

#if RCT_VF_JOG
#define RCT_JOG_NL (SNOVA_n * SNOVA_l)
#define RCT_JOG_VTL ((RCT_JOG_NL + 31) / 32)
#define RCT_JOG_L32 (RCT_JOG_VTL * 32)
#define RCT_JOG_RP ((SNOVA_r + 1) / 2)
#if SNOVA_l != 4
#define RCT_JOG_PKXJOG 1
#else
#define RCT_JOG_PKXJOG 0
#endif
#else
#define RCT_JOG_PKXJOG 0
#endif

#if !defined(RCT_TILE4)
#if RCT_VF_JOG && RCT_HAVE_GFNI && (SNOVA_l == 5) && (SNOVA_q == 16) && (SNOVA_lr <= 32)
#define RCT_TILE4 1
#else
#define RCT_TILE4 0
#endif
#endif
#if (RCT_TILE4 + 0) && RCT_VF_JOG && RCT_HAVE_GFNI && (SNOVA_l == 5) && (SNOVA_q == 16)
#define RCT_VF_TILE4 1
#else
#define RCT_VF_TILE4 0
#endif

#if defined(RCT_OQ_WHIPVEC) && (RCT_OQ_WHIPVEC + 0) && (SNOVA_q != 16)
#define RCT_OQWV 1
#else
#define RCT_OQWV 0
#endif

#if RCT_USE_PSHUFB && defined(RCT_QRP16) && (RCT_QRP16 + 0) && (SNOVA_l == 4)
#define RCT_HOT_QRP16 1
#else
#define RCT_HOT_QRP16 0
#endif
#if RCT_USE_GFNI
#define RCT_GFMUL256(a, b) _mm256_gf2p8mul_epi8((a), (b))
#define RCT_GFMUL128(a, b) _mm_gf2p8mul_epi8((a), (b))
#elif RCT_HOT_QRP16
#define RCT_GFMUL256(a, b) gf16_qrp16_256_byte_mul((a), (b))
#define RCT_GFMUL128(a, b) gf16_qrp16_128_byte_mul((a), (b))
#elif RCT_SIGN_JOG && RCT_HAVE_GFNI
#define RCT_GFMUL256(a, b) _mm256_gf2p8mul_epi8((a), (b))
#define RCT_GFMUL128(a, b) _mm_gf2p8mul_epi8((a), (b))
#elif RCT_SIGN_JOG && defined(RCT_QRP16) && (RCT_QRP16 + 0)
#define RCT_GFMUL256(a, b) gf16_qrp16_256_byte_mul((a), (b))
#define RCT_GFMUL128(a, b) gf16_qrp16_128_byte_mul((a), (b))
#endif
#ifdef RCT_GFMUL256
#define RCT_GFMUL_ANY 1
#else
#define RCT_GFMUL_ANY 0
#endif

#if defined(SNOVA_SIGN_STREAM) && (SNOVA_SIGN_STREAM + 0)
#if SNOVA_q != 16
#error "SNOVA_SIGN_STREAM is q16-only: the T12 block-transpose fold identities assume S-symmetry (q16); odd-q needs a cell-transpose re-verify first — build without SIGN_STREAM=1"
#endif
#if !SNOVA_PK_EXPAND_SHAKE && !defined(SNOVA_ARCH_X86_AVX2)
#error "SNOVA_SIGN_STREAM AES path needs aes128_ctr_zero_at (AVX2); build with AES=0 (SHAKE XOF) or ARCH=x86_avx2"
#endif
#define RCT_SIGN_STREAM 1
#else
#define RCT_SIGN_STREAM 0
#endif

#if defined(SNOVA_KEYGEN_STREAM) && (SNOVA_KEYGEN_STREAM + 0)
#if !SNOVA_PK_EXPAND_SHAKE && !defined(SNOVA_ARCH_X86_AVX2)
#error "SNOVA_KEYGEN_STREAM AES path needs aes128_ctr_zero_at (AVX2); build with AES=0 (SHAKE XOF) or ARCH=x86_avx2"
#endif
#define RCT_KG_STREAM 1
#else
#define RCT_KG_STREAM 0
#endif

#ifdef RCT_BASELINE_SCALAR_HOT
#define RCT_HOT_SIMD RCT_USE_GFNI
#else
#define RCT_HOT_SIMD RCT_USE_SIMD
#endif

#include "platforms/x86_avx2/rct_gf16_cell_avx2.h"

#if !defined(RCT_CM_ACTIVE) && defined(RCT_SQ_CM_LEFT) && (RCT_SQ_CM_LEFT + 0) \
    && RCT_USE_SIMD && (SNOVA_q == 16) && (SNOVA_r == SNOVA_l) \
    && (SNOVA_L == 4)
#define RCT_CML_ONLY 1
#endif

#if defined(RCT_SQ_CM_S2) && (RCT_SQ_CM_S2 + 0) && RCT_USE_SIMD \
    && (SNOVA_q == 16) && (SNOVA_r == SNOVA_l) && (SNOVA_L == 4)
#define RCT_CMS2_ONLY 1
static inline __m128i rct_s2_pcol(int k) {
    return _mm_add_epi8(_mm_setr_epi8(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12),
                        _mm_set1_epi8((char)k));
}
static inline __m128i rct_s2_prow(int k) {
    return _mm_add_epi8(_mm_setr_epi8(0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3),
                        _mm_set1_epi8((char)(4 * k)));
}
#endif

#if defined(RCT_SQ_RIGHT_DUAL) && (RCT_SQ_RIGHT_DUAL + 0) && RCT_USE_SIMD \
    && (SNOVA_q == 16) && (SNOVA_r == SNOVA_l) && (SNOVA_L == 4)
#define RCT_RD_ACTIVE 1
#endif

#if defined(RCT_KG_MULLO) && (RCT_KG_MULLO + 0) && RCT_USE_PSHUFB \
    && (SNOVA_q == 16) && (SNOVA_l == 4)
#define RCT_KGM_ACTIVE 1
#else
#define RCT_KGM_ACTIVE 0
#endif

#if RCT_USE_SIMD && (SNOVA_l == 4) && (SNOVA_lr16 == 1)
#define RCT_VF_MTK2 1
#else
#define RCT_VF_MTK2 0
#endif

#if RCT_USE_GFNI && RCT_VF_MTK2 && (SNOVA_r2 <= 64) && (SNOVA_r <= 7) && !RCT_VF_JOG
#define RCT_VF_EMM 1
#else
#define RCT_VF_EMM 0
#endif
#if RCT_USE_GFNI && RCT_VF_MTK2 && !defined(RCT_AQ_OFF) && !RCT_VF_JOG
#define RCT_VF_AQ 1
#else
#define RCT_VF_AQ 0
#endif

#if RCT_USE_GFNI
#include "platforms/x86_avx2/rct_gf16_gfni_avx2.h"
#elif RCT_USE_PSHUFB
#include "platforms/x86_avx2/rct_gf16_pshufb_avx2.h"
#endif

#if RCT_SIGN_JOG
#define RCT_MATMUL_ADD(a, b, c)      rct_sj_mm_add((a), (b), (c), SNOVA_l, SNOVA_l, SNOVA_l)
#define RCT_MATMUL_ADD_ASEC(a, b, c) rct_sj_mm_add((a), (b), (c), SNOVA_l, SNOVA_l, SNOVA_l)
#define RCT_MATMUL_ADD_BSEC(a, b, c) rct_sj_mm_add_pub((a), (b), (c), SNOVA_l, SNOVA_l, SNOVA_l)
#elif RCT_USE_SIMD && SNOVA_l == 4
#define RCT_MATMUL_ADD(a, b, c) rct_gf4_matmul_add((a), (b), (c))
#define RCT_MATMUL_ADD_ASEC(a, b, c) rct_gf4_matmul_add_sec((a), (b), (c))
#define RCT_MATMUL_ADD_BSEC(a, b, c) rct_gf4_matmul_add((a), (b), (c))
#else
#define RCT_MATMUL_ADD(a, b, c) gf_mat_mul_add((a), (b), (c))
#define RCT_MATMUL_ADD_ASEC(a, b, c) gf_mat_mul_add_sec((a), (b), (c))
#define RCT_MATMUL_ADD_BSEC(a, b, c) gf_mat_mul_add_sec((a), (b), (c))
#endif

#include "platforms/x86_avx2/rct_oddq_avx2.h"

#include "platforms/x86_avx2/rct_verify_kernels_avx2.h"

#include "platforms/generic/rct_gf16_scalar.h"

#if defined(SNOVA_VERIFY_STREAM) && (SNOVA_VERIFY_STREAM + 0) && !RCT_VF_JOG && \
    (((SNOVA_q == 16) && RCT_VF_AQ) || ((SNOVA_q != 16) && RCT_Q_SIMD && RCT_Q_MADD))
#if !SNOVA_PK_EXPAND_SHAKE && !defined(SNOVA_ARCH_X86_AVX2)
#error "SNOVA_VERIFY_STREAM (rct) AES path needs aes128_ctr_zero_at (AVX2); build with AES=0 (SHAKE XOF) or ARCH=x86_avx2"
#endif
#define RCT_VERIFY_STREAM 1
#else
#define RCT_VERIFY_STREAM 0
#endif

#if defined(SNOVA_PKX_PGEN) && (SNOVA_PKX_PGEN + 0)
#if !SNOVA_PK_EXPAND_SHAKE && !defined(SNOVA_ARCH_X86_AVX2)
#error "SNOVA_PKX_PGEN AES path needs aes128_ctr_zero_at (AVX2); build with AES=0 (SHAKE XOF) or ARCH=x86_avx2"
#endif
#if ((SNOVA_q == 16) && ((SNOVA_l2 % 2) != 0)) || RCT_JOG_PKXJOG
#error "SNOVA_PKX_PGEN requires the fused pk_expand gate (q16 needs even l2; JOG l!=4 pkx uses the jogress-layout legacy path) - build without PKX_PGEN=1"
#endif
#define RCT_PKX_PGEN 1
#else
#define RCT_PKX_PGEN 0
#endif

#if RCT_VERIFY_STREAM && (SNOVA_q != 16)
static inline void rct_vf_qs_rowseg(uint16_t *accm, const gf_t *cells, int col0, int ncols) {
    __m256i acc[SNOVA_l][RCT_Q_LR16];
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int g = 0; g < RCT_Q_LR16; g++)
            acc[i1][g] = _mm256_load_si256((const __m256i *)&accm[(i1 * RCT_Q_LR16 + g) * 16]);
    for (int j = 0; j < ncols; ++j) {
        const uint8_t *wb = &rct_qv_wpair[(size_t)(col0 + j) * 2 * 2 * RCT_Q_LRP];
        const gf_t *pcell = cells + (size_t)j * SNOVA_l2;
        for (int h = 0; h < 2; ++h) {
            __m256i w0 = _mm256_load_si256((const __m256i *)(wb + h * 2 * RCT_Q_LRP));
#if RCT_Q_LR16 == 2
            __m256i w1 = _mm256_load_si256((const __m256i *)(wb + h * 2 * RCT_Q_LRP + 32));
#endif
            for (int i1 = 0; i1 < SNOVA_l; i1++) {
                uint16_t pw;
                memcpy(&pw, pcell + i1 * SNOVA_l + 2 * h, 2);
                __m256i pb = _mm256_set1_epi16((short)pw);
                acc[i1][0] = _mm256_add_epi16(acc[i1][0], _mm256_maddubs_epi16(w0, pb));
#if RCT_Q_LR16 == 2
                acc[i1][1] = _mm256_add_epi16(acc[i1][1], _mm256_maddubs_epi16(w1, pb));
#endif
            }
        }
    }
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int g = 0; g < RCT_Q_LR16; g++)
            _mm256_store_si256((__m256i *)&accm[(i1 * RCT_Q_LR16 + g) * 16], acc[i1][g]);
}
#endif

static inline gf_t gf_mat_det(gf_t *a) {
    gf_t det = 0;
#if SNOVA_l == 1
    det = a[0];
#elif SNOVA_l == 2
    det = gf_sub(gf_mult(a[0], a[3]), gf_mult(a[1], a[2]));
#elif SNOVA_l == 3
    det = gf_mult(a[0], gf_sub(gf_mult(a[4], a[8]), gf_mult(a[5], a[7])));
    gf_set_add(&det, gf_mult(a[1], gf_sub(gf_mult(a[5], a[6]), gf_mult(a[3], a[8]))));
    gf_set_add(&det, gf_mult(a[2], gf_sub(gf_mult(a[3], a[7]), gf_mult(a[4], a[6]))));
#elif SNOVA_l == 4
    gf_t det_l, det_r;
#define DET_L(x, y) det_l = gf_sub(gf_mult(a[x], a[4 + y]), gf_mult(a[y], a[4 + x]))
#define DET_R(x, y) det_r = gf_sub(gf_mult(a[8 + x], a[12 + y]), gf_mult(a[8 + y], a[12 + x]))
#define DET22(x1, y1, x2, y2) DET_L(x1, y1); DET_R(x2, y2); gf_set_add(&det, gf_mult(det_l, det_r))
    DET22(0, 1, 2, 3); DET22(0, 2, 3, 1); DET22(0, 3, 1, 2);
    DET22(1, 2, 0, 3); DET22(1, 3, 2, 0); DET22(2, 3, 0, 1);
#undef DET_R
#undef DET22
#undef DET_L
#elif SNOVA_l == 5
    gf_t det_l, det_r;
#define DET_L(x, y) det_l = gf_sub(gf_mult(a[x], a[5 + y]), gf_mult(a[y], a[5 + x]))
#define DET_R2(x, y, z) gf_mult(gf_sub(gf_mult(a[10 + x], a[15 + y]), gf_mult(a[10 + y], a[15 + x])), a[20 + z])
#define DET_R3(x, y, z) det_r = gf_add(DET_R2(x, y, z), gf_add(DET_R2(y, z, x), DET_R2(z, x, y)))
#define DET23(x1, y1, x2, y2, z2) DET_L(x1, y1); DET_R3(x2, y2, z2); gf_set_add(&det, gf_mult(det_l, det_r))
    DET23(0, 1, 2, 3, 4); DET23(0, 2, 3, 1, 4); DET23(0, 3, 1, 2, 4); DET23(0, 4, 1, 3, 2);
    DET23(1, 2, 0, 3, 4); DET23(1, 3, 2, 0, 4); DET23(1, 4, 2, 3, 0);
    DET23(2, 3, 0, 1, 4); DET23(2, 4, 0, 3, 1); DET23(3, 4, 2, 0, 1);
#undef DET_R2
#undef DET_R3
#undef DET23
#undef DET_L
#else
#error "Unsupported rank"
#endif
    return det;
}

static void rct_init_gf_tables(void) {
#if SNOVA_q == 16
    uint8_t F_star[15] = {1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9};
    for (int i1 = 0; i1 < 16; i1++) {
        rct_multtab[i1] = 0;
        rct_multtab[i1 * SNOVA_q] = 0;
    }
    for (int i1 = 0; i1 < SNOVA_q - 1; i1++)
        for (int j1 = 0; j1 < SNOVA_q - 1; j1++)
            rct_multtab[F_star[i1] * SNOVA_q + F_star[j1]] = F_star[(i1 + j1) % (SNOVA_q - 1)];
    for (int i1 = 0; i1 < SNOVA_q; i1++)
        for (int j1 = 0; j1 < SNOVA_q; j1++)
            rct_addtab[i1 * SNOVA_q + j1] = (i1 ^ j1);
#else
    for (int i1 = 0; i1 < SNOVA_q; i1++)
        for (int j1 = 0; j1 < SNOVA_q; j1++) {
            rct_multtab[i1 * SNOVA_q + j1] = (i1 * j1) % SNOVA_q;
            rct_addtab[i1 * SNOVA_q + j1] = (i1 + j1) % SNOVA_q;
        }
#endif
}

static void rct_set_S(gf_t *gf_S1) {
#if SNOVA_q == 16
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int j1 = 0; j1 < SNOVA_l; j1++)
            gf_S1[i1 * SNOVA_l + j1] = 8 - (i1 + j1);
#if SNOVA_l == 5
    gf_S1[SNOVA_l2 - 1] = 9;
#endif
#else
    for (int i1 = 0; i1 < SNOVA_l; i1++)
        for (int j1 = 0; j1 < SNOVA_l; j1++)
            gf_S1[i1 * SNOVA_l + j1] = (SNOVA_Q_A + i1 + j1) & SNOVA_Q_B;
    gf_S1[SNOVA_l2 - 1] = SNOVA_Q_C;
#endif
}

static void rct_gen_S_array(void) {
    memset(rct_S, 0, sizeof(rct_S));
    for (int i1 = 0; i1 < SNOVA_l; i1++) rct_S[i1 * SNOVA_l + i1] = 1;
#if SNOVA_l > 1
    rct_set_S(&rct_S[1 * SNOVA_l2]);
    for (int i1 = 2; i1 < SNOVA_l; i1++)
        gf_mat_mul(&rct_S[i1 * SNOVA_l2], &rct_S[1 * SNOVA_l2], &rct_S[(i1 - 1) * SNOVA_l2]);
#endif
}

static void convert_bytes_to_GF(gf_t *gf_array, const uint8_t *byte_array, size_t num) {
#if SNOVA_q != 16
    for (size_t idx = 0; idx < num; idx++)
        gf_array[idx] = byte_array[idx] % SNOVA_q;
#else
    for (size_t idx = 0; idx < num / 2; idx++) {
        gf_array[2 * idx] = (byte_array[idx] & 0xf) % SNOVA_q;
        gf_array[2 * idx + 1] = (byte_array[idx] >> 4) % SNOVA_q;
    }
    if (num & 1) gf_array[num - 1] = (byte_array[num / 2] & 0xf) % SNOVA_q;
#endif
}

#if SNOVA_q == 16
static inline void rct_unpack_nib_seg(gf_t *dst, const uint8_t *src, size_t nb) {
    size_t i = 0;
#if RCT_USE_SIMD
    const __m256i m0f = _mm256_set1_epi8(0x0f);
    for (; i + 32 <= nb; i += 32) {
        __m256i b = _mm256_loadu_si256((const __m256i *)(src + i));
        __m256i lo = _mm256_and_si256(b, m0f);
        __m256i hi = _mm256_and_si256(_mm256_srli_epi16(b, 4), m0f);
        __m256i u0 = _mm256_unpacklo_epi8(lo, hi);
        __m256i u1 = _mm256_unpackhi_epi8(lo, hi);
        _mm256_storeu_si256((__m256i *)(dst + 2 * i), _mm256_permute2x128_si256(u0, u1, 0x20));
        _mm256_storeu_si256((__m256i *)(dst + 2 * i + 32), _mm256_permute2x128_si256(u0, u1, 0x31));
    }
    if (i + 16 <= nb) {
        __m128i b = _mm_loadu_si128((const __m128i *)(src + i));
        __m128i lo = _mm_and_si128(b, _mm256_castsi256_si128(m0f));
        __m128i hi = _mm_and_si128(_mm_srli_epi16(b, 4), _mm256_castsi256_si128(m0f));
        _mm_storeu_si128((__m128i *)(dst + 2 * i), _mm_unpacklo_epi8(lo, hi));
        _mm_storeu_si128((__m128i *)(dst + 2 * i + 16), _mm_unpackhi_epi8(lo, hi));
        i += 16;
    }
#endif
    for (; i < nb; ++i) {
        dst[2 * i] = (gf_t)(src[i] & 0x0f);
        dst[2 * i + 1] = (gf_t)(src[i] >> 4);
    }
}
#else
static inline void rct_unpack_modq_seg(gf_t *dst, const uint8_t *src, size_t ngf) {
    for (size_t i = 0; i < ngf; ++i) dst[i] = src[i] % SNOVA_q;
}
#endif

static void compress_gf(uint8_t *byte_array, const gf_t *gf_array, size_t num) {
    size_t idx = 0, out_idx = 0;
    size_t num_bytes = BYTES_GF(num);
    do {
        uint64_t val = 0, fact = 1;
        int i1 = 0;
        while (i1 < PACK_GF && idx < num) {
            val += fact * (gf_array[idx] % SNOVA_q);
            idx++; i1++; fact *= SNOVA_q;
        }
        i1 = (i1 + 1) / 2;
        int j1 = 0;
        while (j1 < PACK_BYTES && out_idx < num_bytes) {
            byte_array[out_idx] = val & 0xff;
            out_idx++; val = val >> 8; j1++;
        }
    } while (idx < num);
}

static int expand_gf(gf_t *gf_array, const uint8_t *byte_array, size_t num) {
    size_t num_bytes = BYTES_GF(num);
    size_t idx = 0, out_idx = 0;
    uint64_t val;
    uint64_t res = 0;
    do {
        val = 0;
        int i1 = 0;
        while (i1 < PACK_BYTES && idx < num_bytes) {
            val = val ^ ((uint64_t)(byte_array[idx]) << (8 * i1));
            idx++; i1++;
        }
        int j1 = 0;
        while (j1 < PACK_GF && out_idx < num) {
            gf_array[out_idx] = val % SNOVA_q;
            val = val / SNOVA_q;
            out_idx++; j1++;
        }
        res |= val;
    } while (out_idx < num);
#if SNOVA_q == 16
    if (num & 1) return byte_array[num / 2] & 0xF0;
#endif
    return res != 0;
}

static void compress_pk(uint8_t *pk, gf_t *P22) { compress_gf(pk, P22, NUMGF_PK); }
static int expand_pk(gf_t *P22, const uint8_t *pk) { return expand_gf(P22, pk, NUMGF_PK); }

static void rct_public_xof(const uint8_t seed[16], uint8_t *out, size_t outlen) {
#if SNOVA_PK_EXPAND_SHAKE
    size_t padded = (outlen + 7u) & ~(size_t)7u;
#if defined(SNOVA_ARCH_X86_AVX2)
    snova_shake_opt(seed, 16, (uint64_t *)out, padded);
#else
    snova_shake_ref(seed, 16, (uint64_t *)out, padded);
#endif
#else
    aes128_ctr_zero(out, outlen, seed);
#endif
}

static void expand_public(gf_t *P_matrix, const uint8_t *seed) {
#if SNOVA_WRAPPER_STACK
    _Alignas(32) uint8_t pk_bytes[((NUM_GEN_PUB_BYTES + 15) & ~(size_t)7u)];
#else
    _Alignas(8) static uint8_t pk_bytes[((NUM_GEN_PUB_BYTES + 15) & ~(size_t)7u)];
#endif
    rct_public_xof(seed, pk_bytes, NUM_GEN_PUB_BYTES);
    convert_bytes_to_GF(P_matrix, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);
}

#if RCT_SIGN_STREAM || RCT_KG_STREAM || RCT_VERIFY_STREAM || RCT_PKX_PGEN
#include "pgen.h"
#endif

static void hash_combined(uint8_t *hash_out, const uint8_t *digest, const size_t len_digest,
                          const uint8_t *pk_seed, const uint8_t *salt) {
    keccak_ctx state;
    shake_init(&state, 256);
#if HASH_PK
    shake_absorb(&state, pk_seed, BYTES_PK_HASH);
#else
    shake_absorb(&state, pk_seed, SEED_LENGTH_PUBLIC);
#endif
    shake_absorb(&state, digest, len_digest);
    shake_absorb(&state, salt, BYTES_SALT);
    shake_finalize(&state);
    shake_squeeze(&state, hash_out, BYTES_HASH);
}

static inline void gen_a_FqS(gf_t *Qm, gf_t *q) {
#if ROUND2_T12
    if (!q[SNOVA_l - 1]) q[SNOVA_l - 1] = SNOVA_q - (q[0] + (q[0] == 0));
#endif
    for (int i1 = 0; i1 < SNOVA_l2; i1++) {
        gf_t sum = 0;
        for (int j1 = 0; j1 < SNOVA_l; j1++)
            gf_set_add(&sum, gf_mult(q[j1], rct_S[j1 * SNOVA_l2 + i1]));
        Qm[i1] = sum;
    }
}

static inline void gen_a_FqS_sec(gf_t *Qm, gf_t *q) {
#if ROUND2_T12
    uint32_t nz_last = ct_gf_nz(q[SNOVA_l - 1]);
    gf_t fallback = (gf_t)(SNOVA_q - (q[0] + (1u - ct_gf_nz(q[0]))));
    q[SNOVA_l - 1] = ct_gf_sel(nz_last, q[SNOVA_l - 1], fallback);
#endif
#if (RCT_USE_GFNI || RCT_HOT_QRP16) && (SNOVA_l == 4)
    __m128i acc = RCT_GFMUL128(_mm_set1_epi8((char)q[0]), _mm_loadu_si128((const __m128i *)&rct_S[0]));
    acc = _mm_xor_si128(acc, RCT_GFMUL128(_mm_set1_epi8((char)q[1]), _mm_loadu_si128((const __m128i *)&rct_S[SNOVA_l2])));
    acc = _mm_xor_si128(acc, RCT_GFMUL128(_mm_set1_epi8((char)q[2]), _mm_loadu_si128((const __m128i *)&rct_S[2 * SNOVA_l2])));
    acc = _mm_xor_si128(acc, RCT_GFMUL128(_mm_set1_epi8((char)q[3]), _mm_loadu_si128((const __m128i *)&rct_S[3 * SNOVA_l2])));
    _mm_storeu_si128((__m128i *)Qm, rct_gfni_cleanup128(acc));
#elif RCT_USE_PSHUFB && (SNOVA_l == 4)
    __m128i acc = rct_gf16_mul128_sec(_mm_set1_epi8((char)q[0]), _mm_loadu_si128((const __m128i *)&rct_S[0]));
    acc = _mm_xor_si128(acc, rct_gf16_mul128_sec(_mm_set1_epi8((char)q[1]), _mm_loadu_si128((const __m128i *)&rct_S[SNOVA_l2])));
    acc = _mm_xor_si128(acc, rct_gf16_mul128_sec(_mm_set1_epi8((char)q[2]), _mm_loadu_si128((const __m128i *)&rct_S[2 * SNOVA_l2])));
    acc = _mm_xor_si128(acc, rct_gf16_mul128_sec(_mm_set1_epi8((char)q[3]), _mm_loadu_si128((const __m128i *)&rct_S[3 * SNOVA_l2])));
    _mm_storeu_si128((__m128i *)Qm, acc);
#else
    for (int i1 = 0; i1 < SNOVA_l2; i1++) {
        gf_t sum = 0;
        for (int j1 = 0; j1 < SNOVA_l; j1++)
            gf_set_add(&sum, gf_mult_sec(q[j1], rct_S[j1 * SNOVA_l2 + i1]));
        Qm[i1] = sum;
    }
#endif
}

#define SK_BLOCK_SIZE 32
static void expand_T12(gf_t *T12, const uint8_t *seed) {
    gf_t T12coef[SNOVA_o * SNOVA_v * SNOVA_l];
    gf_t sk_data[SK_BLOCK_SIZE];
    keccak_ctx state;
    shake_init(&state, 256);
    shake_absorb(&state, seed, SEED_LENGTH_PRIVATE);
    shake_finalize(&state);

    size_t idx = SK_BLOCK_SIZE, t_idx = 0;
    while (t_idx < (size_t)SNOVA_o * SNOVA_v * SNOVA_l) {
        if (idx >= SK_BLOCK_SIZE) {
            shake_squeeze(&state, sk_data, SK_BLOCK_SIZE);
            idx = 0;
        }
#if SNOVA_q != 16
        {
            int accept = (sk_data[idx] < SNOVA_REJECTION_LIMIT);
            SNOVA_CT_DECLASSIFY(&accept, sizeof accept);
            if (accept) {
                T12coef[t_idx] = sk_data[idx] % SNOVA_q;
                t_idx++;
            }
        }
#else
        T12coef[t_idx] = sk_data[idx] & 0xf; t_idx++;
        T12coef[t_idx] = sk_data[idx] >> 4; t_idx++;
#endif
        idx++;
    }
    for (size_t i1 = 0; i1 < (size_t)SNOVA_o * SNOVA_v; i1++)
        gen_a_FqS_sec(&T12[i1 * SNOVA_l2], &T12coef[i1 * SNOVA_l]);
    SNOVA_CLEAR_OBJ(T12coef);
    SNOVA_CLEAR_OBJ(sk_data);
}

static inline void be_invertible_by_add_aS(gf_t *mat, const gf_t *orig, const int l1, const int l2) {
    memcpy(mat, orig, l1 * l2);
#if ABQ_ALG2
    if ((l1 == SNOVA_l) && (l2 == SNOVA_l))
        if (gf_mat_det(mat) == 0)
            for (gf_t f1 = 1; f1 < SNOVA_q; f1++) {
#if SNOVA_l > 1
                for (int i1 = 0; i1 < SNOVA_l2; i1++)
                    gf_set_add(&mat[i1], gf_mult(f1, rct_S[SNOVA_l2 + i1]));
#else
                mat[0] = 1;
#endif
                if (gf_mat_det(mat) != 0) break;
            }
#endif
}

static void gen_ABQ(gf_t *A, gf_t *Am, gf_t *Bm, gf_t *Q1m, gf_t *Q2m) {
    gf_t *B = A + SNOVA_o * SNOVA_alpha * SNOVA_r2;
    gf_t *q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
    gf_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;
    for (size_t idx = 0; idx < (size_t)SNOVA_o * SNOVA_alpha; idx++) {
        be_invertible_by_add_aS(&Am[idx * SNOVA_r2], &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
        be_invertible_by_add_aS(&Bm[idx * SNOVA_lr], &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);
        gen_a_FqS(&Q1m[idx * SNOVA_l2], &q1[idx * SNOVA_l]);
        gen_a_FqS(&Q2m[idx * SNOVA_l2], &q2[idx * SNOVA_l]);
    }
}

#define ABQ_RAW_N (SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l))
static gf_t rct_fixed_abq[ABQ_RAW_N];
#if FIXED_ABQ
static gf_t rct_fixed_Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
static gf_t rct_fixed_Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
static gf_t rct_fixed_Q1[SNOVA_o * SNOVA_alpha * SNOVA_l2];
static gf_t rct_fixed_Q2[SNOVA_o * SNOVA_alpha * SNOVA_l2];
#endif
#if !SNOVA_WRAPPER_STACK && !(RCT_KG_STREAM && !RCT_JOG_PKXJOG)
static gf_t rct_pub_Pmatrix[NUM_PUB_GF];
#endif
static int rct_inited = 0;

static void gen_fixed_ABQ(const char *abq_seed) {
    uint8_t rng_out[ABQ_RAW_N];
    shake256(rng_out, ABQ_RAW_N, (const uint8_t *)abq_seed, (size_t)strlen(abq_seed));
    convert_bytes_to_GF(rct_fixed_abq, rng_out, ABQ_RAW_N);
}

#if RCT_VF_TILE4
static void rct_build_t4(void);
#endif
static void rct_init(void) {
    if (rct_inited) return;
    rct_init_gf_tables();
#if RCT_USE_PSHUFB
    rct_build_vtl();
#endif
    rct_gen_S_array();
#if RCT_VF_MTK2 || RCT_VF_JOG
    rct_build_mtk2();
#endif
#if RCT_VF_MTK2 && RCT_USE_GFNI
    rct_build_vf_gfni();
#endif
#if RCT_VF_AQ
    rct_build_vf_aq();
#endif
#if RCT_VF_TILE4
    rct_build_t4();
#endif
    gen_fixed_ABQ("SNOVA_ABQ");
#if FIXED_ABQ
    gen_ABQ(rct_fixed_abq, rct_fixed_Am, rct_fixed_Bm, rct_fixed_Q1, rct_fixed_Q2);
#endif
    rct_inited = 1;
}

#include "platforms/x86_avx2/rct_fold_engine.h"

#include "snova_rect_keygen.h"

#include "snova_rect_sign.h"

#include "snova_rect_verify.h"

#include "snova_rect_api.h"

#endif
