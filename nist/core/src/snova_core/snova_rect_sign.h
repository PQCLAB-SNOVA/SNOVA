#ifndef SNOVA_RECT_SIGN_H
#define SNOVA_RECT_SIGN_H

#ifndef RCT_SKX_SLIM
#define RCT_SKX_SLIM 0
#endif
typedef struct {
    uint8_t sk[BYTES_SK];
    gf_t T12[SNOVA_o * SNOVA_v * SNOVA_l2 + 16];
#if RCT_SIGN_STREAM
    gf_t abq[ABQ_RAW_N];
#elif RCT_SKX_SLIM
    gf_t P11[SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2 + 32];
    gf_t abq[ABQ_RAW_N];
#else
    gf_t P_matrix[NUM_PUB_GF];
#endif
    gf_t F21[SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2];
    gf_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2];
    gf_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2 + 16];
    gf_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr + 16];
    gf_t Q1[SNOVA_o * SNOVA_alpha * SNOVA_l2];
    gf_t Q2[SNOVA_o * SNOVA_alpha * SNOVA_l2 + 16];
} rct_skx_t;

#if defined(RCT_OQ_DEVFLOW) && (RCT_OQ_DEVFLOW + 0) && defined(RCT_Q_SIMD) && \
    (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC
#define RCT_OQDF 1
#else
#define RCT_OQDF 0
#endif
#if RCT_OQDF
static _Alignas(32) uint16_t rct_oq_whip_w[SNOVA_l * SNOVA_v * RCT_Q_LRP];
static _Alignas(32) uint16_t rct_oq_sum_t1u[SNOVA_m1 * SNOVA_l2 * SNOVA_r2 + 64];
static _Alignas(32) uint16_t rct_oq_P11u[SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2];
static _Alignas(32) uint16_t rct_oq_F21u[SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2];
static _Alignas(32) uint16_t rct_oq_F12u[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2];
static inline void rct_oq_expand_u16(uint16_t *dst, const gf_t *src, int n) {
    int i = 0;
    for (; i + 16 <= n; i += 16)
        _mm256_storeu_si256((__m256i *)(dst + i),
            _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(src + i))));
    for (; i < n; ++i) dst[i] = src[i];
}
#endif

#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC
enum { RCT_GN = SNOVA_o * SNOVA_lr, RCT_GPAD = (SNOVA_o * SNOVA_lr / 16 + 1) * 16 };
#endif
#if SNOVA_Q == 16 && RCT_GFMUL_ANY && !defined(RCT_GAUSS_SCALAR)
enum { RCT_SNB = (SNOVA_o * SNOVA_lr + 1 + 31) / 32 };
#endif
#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
enum { RCT_CM_LR16 = SNOVA_lr32 / 16,
       RCT_CM_OLR16 = (SNOVA_o * SNOVA_lr + 15) / 16,
       RCT_CM_OLR = RCT_CM_OLR16 * 16,
       RCT_CM_OLR32 = (RCT_CM_OLR + 31) / 32 * 32,
       RCT_CM_OLR32N = RCT_CM_OLR32 / 32 };
#if (SNOVA_r == SNOVA_l) && defined(RCT_SQ_CM_TIGHT) && (RCT_SQ_CM_TIGHT + 0)
enum { RCT_CMW = SNOVA_lr, RCT_CMW16 = SNOVA_lr / 16 };
#else
enum { RCT_CMW = SNOVA_lr32, RCT_CMW16 = SNOVA_lr32 / 16 };
#endif
#endif

typedef struct rct_cm_ctx {
    uint16_t *cm_whip; uint8_t *cm_whipb;
    const uint16_t *Amx, *Bmx, *Q1x, *Q2x, *q1x, *q2x;
} rct_cm_ctx;

typedef struct rct_sign_ctx {
    rct_skx_t *skx;
    const gf_t *T12, *P11, *aptr, *F21, *F12, *Am, *Bm, *Q1, *Q2, *q1, *q2;
    gf_t *signature_in_GF; gf_t *hash_in_GF16; gf_t *Fvv;
    gf_t (*gauss)[SNOVA_o * SNOVA_lr + 1 + 64];
    gf_t *sum_t1; gf_t *whipped_sig; gf_t *whipped_F21, *whipped_F12;
    rct_cm_ctx *cm; uint8_t num_sign; int flag_redo; uint64_t _pa, _pb;
} rct_sign_ctx;

#include "platforms/x86_avx2/rct_sign_gfni.h"
#include "platforms/x86_avx2/rct_sign_pshufb.h"
#include "platforms/x86_avx2/rct_sign_mullo.h"
#include "platforms/x86_avx2/rct_sign_cm.h"
#include "platforms/x86_avx2/rct_sign_oddq.h"
#include "platforms/x86_avx2/rct_sign_jog.h"
#include "platforms/generic/rct_sign_scalar.h"
static void rct_sk_expand(const uint8_t *sk, rct_skx_t *skx) {
    rct_init();
    uint64_t _pa,_pb; (void)_pa; (void)_pb;
    const uint8_t *seed = sk;
    memcpy(skx->sk, sk, BYTES_SK);
    memset(skx->T12 + SNOVA_o * SNOVA_v * SNOVA_l2, 0, 16);
    memset(skx->Am + SNOVA_o * SNOVA_alpha * SNOVA_r2, 0, 16);
    memset(skx->Bm + SNOVA_o * SNOVA_alpha * SNOVA_lr, 0, 16);
    memset(skx->Q2 + SNOVA_o * SNOVA_alpha * SNOVA_l2, 0, 16);

    expand_T12(skx->T12, seed + SEED_LENGTH_PUBLIC);
#if RCT_SIGN_STREAM || RCT_SKX_SLIM
    gf_t *P_stream = (gf_t *)malloc((size_t)NUM_GEN_PUB_GF * sizeof(gf_t));
    if (!P_stream) { memset(skx, 0, sizeof(*skx)); return; }
#if RCT_SIGN_STREAM
    snova_pgen_fill_pblocks(seed, P_stream);
#else
    expand_public(P_stream, seed);
#endif
    gf_t *T12 = skx->T12;
    gf_t *P11 = P_stream;
    gf_t *P12 = P_stream + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
    gf_t *P21 = P_stream + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;
#else
    expand_public(skx->P_matrix, seed);

    gf_t *T12 = skx->T12;
    gf_t *P11 = skx->P_matrix;
    gf_t *P12 = skx->P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
    gf_t *P21 = skx->P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;
#endif
    gf_t *F21 = skx->F21;
    gf_t *F12 = skx->F12;
    memset(F21, 0, sizeof(skx->F21));
    memset(F12, 0, sizeof(skx->F12));

    RCT_PT(_pa);
#if RCT_Q_SIMD
    rct_skx_fold_F_oddq(F21, F12, T12, P11, P12, P21);
#else
#if RCT_SIGN_JOG
    rct_skx_fold_F_jog(F21, F12, T12, P11);
#elif RCT_USE_GFNI && (SNOVA_l == 4)
    rct_skx_fold_F_gfni(F21, F12, T12, P11);
#elif RCT_USE_PSHUFB && (SNOVA_l == 4) && defined(RCT_MULLO) && (RCT_MULLO + 0)
    rct_skx_fold_F_mullo(F21, F12, T12, P11);
#elif RCT_HOT_QRP16 && (SNOVA_l == 4)
    rct_skx_fold_F_qrp16(F21, F12, T12, P11);
#else
    rct_skx_fold_F_scalar(F21, F12, T12, P11);
#endif

    for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) gf_set_add(&F12[i1], P12[i1]);
    for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) gf_set_add(&F21[i1], P21[i1]);
#endif
    RCT_PT(_pb); RCT_PACC(3,_pa,_pb);

#if RCT_SIGN_STREAM || RCT_SKX_SLIM
#if RCT_SKX_SLIM && !RCT_SIGN_STREAM
    memcpy(skx->P11, P_stream, (size_t)SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2 * sizeof(gf_t));
    memset(skx->P11 + (size_t)SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2, 0, 32 * sizeof(gf_t));
#endif
    free(P_stream);
    gf_t *aptr = skx->abq;
#else
    gf_t *aptr = skx->P_matrix + SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2;
#endif
#if FIXED_ABQ
    memcpy(aptr, rct_fixed_abq, sizeof(rct_fixed_abq));
    memcpy(skx->Am, rct_fixed_Am, sizeof(rct_fixed_Am));
    memcpy(skx->Bm, rct_fixed_Bm, sizeof(rct_fixed_Bm));
    memcpy(skx->Q1, rct_fixed_Q1, sizeof(rct_fixed_Q1));
    memcpy(skx->Q2, rct_fixed_Q2, sizeof(rct_fixed_Q2));
#else
    gen_ABQ(aptr, skx->Am, skx->Bm, skx->Q1, skx->Q2);
#endif
}

static int rct_sign_expanded(rct_skx_t *skx, uint8_t *sig, const uint8_t *digest,
                             const size_t len_digest, const uint8_t *salt) {
    rct_init();
    uint64_t _pa,_pb; (void)_pa; (void)_pb;
    const uint8_t *seed = skx->sk;

    gf_t *T12 = skx->T12;
#if RCT_SIGN_STREAM
    size_t p11_gf = (size_t)SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
    gf_t *P11 = (gf_t *)malloc((p11_gf + 32) * sizeof(gf_t));
    if (!P11) { memset(sig, 0, BYTES_SIGNATURE); return -1; }
    memset(P11 + p11_gf, 0, 32 * sizeof(gf_t));
    snova_pgen_fill_p11(seed, P11);
    gf_t *aptr = skx->abq;
#elif RCT_SKX_SLIM
    gf_t *P11 = skx->P11;
    gf_t *aptr = skx->abq;
#else
    gf_t *P11 = skx->P_matrix;
    gf_t *aptr = skx->P_matrix + SNOVA_m1 * (SNOVA_v * SNOVA_v + 2 * SNOVA_v * SNOVA_o) * SNOVA_l2;
#endif
    gf_t *F21 = skx->F21;
    gf_t *F12 = skx->F12;
    gf_t *Am = skx->Am;
    gf_t *Bm = skx->Bm;
    gf_t *Q1 = skx->Q1;
    gf_t *Q2 = skx->Q2;
    gf_t *q1 = aptr + SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr);
    gf_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

    gf_t hash_in_GF16[GF16_HASH];
    uint8_t sign_hashb[BYTES_HASH];
#if HASH_PK
    hash_combined(sign_hashb, digest, len_digest, seed + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE, salt);
#else
    hash_combined(sign_hashb, digest, len_digest, seed, salt);
#endif
    expand_gf(hash_in_GF16, sign_hashb, GF16_HASH);

    RCT_SCRATCH _Alignas(32) gf_t gauss[SNOVA_o * SNOVA_lr][SNOVA_o * SNOVA_lr + 1 + 64];
    gf_t solution[SNOVA_o * SNOVA_lr] = {0};
#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC
    enum { RCT_GN = SNOVA_o * SNOVA_lr, RCT_GPAD = (SNOVA_o * SNOVA_lr / 16 + 1) * 16 };
    RCT_SCRATCH _Alignas(32) uint16_t gu[RCT_GN][RCT_GPAD];
    _Alignas(32) uint16_t sol16[RCT_GPAD];
#endif
#if SNOVA_Q == 16 && RCT_GFMUL_ANY && !defined(RCT_GAUSS_SCALAR)
    enum { RCT_SNB = (SNOVA_o * SNOVA_lr + 1 + 31) / 32 };
    _Alignas(32) gf_t solpad[RCT_SNB * 32 + 32];
#endif
    RCT_SCRATCH _Alignas(32) gf_t signature_in_GF[SNOVA_n * SNOVA_lr + 32];
    memset(signature_in_GF, 0, sizeof(signature_in_GF));
#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
    enum { RCT_CM_LR16 = SNOVA_lr32 / 16,
           RCT_CM_OLR16 = (SNOVA_o * SNOVA_lr + 15) / 16,
           RCT_CM_OLR = RCT_CM_OLR16 * 16,
           RCT_CM_OLR32 = (RCT_CM_OLR + 31) / 32 * 32,
           RCT_CM_OLR32N = RCT_CM_OLR32 / 32 };
#if (SNOVA_r == SNOVA_l) && defined(RCT_SQ_CM_TIGHT) && (RCT_SQ_CM_TIGHT + 0)
    enum { RCT_CMW = SNOVA_lr, RCT_CMW16 = SNOVA_lr / 16 };
#else
    enum { RCT_CMW = SNOVA_lr32, RCT_CMW16 = SNOVA_lr32 / 16 };
#endif
    RCT_SCRATCH _Alignas(32) uint16_t rct_cm_whip[SNOVA_l * SNOVA_v * RCT_CMW];
#if RCT_USE_GFNI
    RCT_SCRATCH _Alignas(32) uint8_t rct_cm_whipb[SNOVA_l * SNOVA_v * RCT_CMW];
#endif
#endif
#if defined(RCT_CM_ACTIVE) || defined(RCT_CMS3_ONLY)
    RCT_SCRATCH _Alignas(32) uint16_t rct_cm_Amx[SNOVA_o * SNOVA_alpha * SNOVA_r2];
    RCT_SCRATCH _Alignas(32) uint16_t rct_cm_Bmx[SNOVA_o * SNOVA_alpha * SNOVA_lr];
    RCT_SCRATCH _Alignas(32) uint16_t rct_cm_Q1x[SNOVA_o * SNOVA_alpha * SNOVA_l2];
    RCT_SCRATCH _Alignas(32) uint16_t rct_cm_Q2x[SNOVA_o * SNOVA_alpha * SNOVA_l2];
    RCT_SCRATCH _Alignas(32) uint16_t rct_cm_q1x[SNOVA_o * SNOVA_alpha * SNOVA_l];
    RCT_SCRATCH _Alignas(32) uint16_t rct_cm_q2x[SNOVA_o * SNOVA_alpha * SNOVA_l];
    rct_cm_expand_arr(rct_cm_Amx, Am, SNOVA_o * SNOVA_alpha * SNOVA_r2);
    rct_cm_expand_arr(rct_cm_Bmx, Bm, SNOVA_o * SNOVA_alpha * SNOVA_lr);
    rct_cm_expand_arr(rct_cm_Q1x, Q1, SNOVA_o * SNOVA_alpha * SNOVA_l2);
    rct_cm_expand_arr(rct_cm_Q2x, Q2, SNOVA_o * SNOVA_alpha * SNOVA_l2);
    rct_cm_expand_arr(rct_cm_q1x, q1, SNOVA_o * SNOVA_alpha * SNOVA_l);
    rct_cm_expand_arr(rct_cm_q2x, q2, SNOVA_o * SNOVA_alpha * SNOVA_l);
#endif
#if RCT_OQDF
    rct_oq_expand_u16(rct_oq_P11u, P11, SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2);
    rct_oq_expand_u16(rct_oq_F21u, F21, SNOVA_m1 * SNOVA_o * SNOVA_v * SNOVA_l2);
    rct_oq_expand_u16(rct_oq_F12u, F12, SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2);
#endif

#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
    rct_cm_ctx CM;
    CM.cm_whip = rct_cm_whip;
#if RCT_USE_GFNI
    CM.cm_whipb = rct_cm_whipb;
#else
    CM.cm_whipb = (uint8_t *)0;
#endif
#if defined(RCT_CM_ACTIVE) || defined(RCT_CMS3_ONLY)
    CM.Amx = rct_cm_Amx; CM.Bmx = rct_cm_Bmx; CM.Q1x = rct_cm_Q1x;
    CM.Q2x = rct_cm_Q2x; CM.q1x = rct_cm_q1x; CM.q2x = rct_cm_q2x;
#else
    CM.Amx = CM.Bmx = CM.Q1x = CM.Q2x = CM.q1x = CM.q2x = (const uint16_t *)0;
#endif
#endif
    rct_sign_ctx C;
    rct_sign_ctx *c = &C;
    c->skx = skx;
    c->T12 = T12; c->P11 = P11; c->aptr = aptr;
    c->F21 = F21; c->F12 = F12;
    c->Am = Am; c->Bm = Bm; c->Q1 = Q1; c->Q2 = Q2; c->q1 = q1; c->q2 = q2;
    c->signature_in_GF = signature_in_GF;
    c->hash_in_GF16 = hash_in_GF16;
    c->gauss = gauss;
#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
    c->cm = &CM;
#else
    c->cm = (rct_cm_ctx *)0;
#endif
    c->_pa = c->_pb = 0; c->num_sign = 0; c->flag_redo = 0;
    (void)c;

    int flag_redo = 1;
    uint8_t num_sign = 0;

    SNOVA_DUDECT_RETRY_RESET();
    int sign_rc = 0;
    do {
        memset(gauss, 0, sizeof(gauss));
        num_sign++;
        SNOVA_DUDECT_RETRY_TICK();
        if (num_sign == 255) { memset(sig, 0, BYTES_SIGNATURE);
            sign_rc = -1; goto sign_cleanup; }
        flag_redo = 0;

        uint8_t vinegar_in_byte[NUM_GEN_SEC_BYTES];
        keccak_ctx v_instance;
        shake_init(&v_instance, 256);
        shake_absorb(&v_instance, seed + SEED_LENGTH_PUBLIC, SEED_LENGTH_PRIVATE);
        shake_absorb(&v_instance, digest, BYTES_DIGEST);
        shake_absorb(&v_instance, salt, BYTES_SALT);
        shake_absorb(&v_instance, &num_sign, 1);
        shake_finalize(&v_instance);
        shake_squeeze(&v_instance, vinegar_in_byte, NUM_GEN_SEC_BYTES);

        expand_gf(signature_in_GF, vinegar_in_byte, SNOVA_v * SNOVA_lr);

        gf_t Fvv_in_GF16Matrix[SNOVA_o * SNOVA_lr] = {0};
#if !RCT_OQDF
        RCT_SCRATCH _Alignas(32) gf_t sum_t1[SNOVA_m1 * SNOVA_l2 * SNOVA_r2 + 64];
        memset(sum_t1, 0, sizeof(sum_t1));
        RCT_SCRATCH _Alignas(32) gf_t whipped_sig[SNOVA_l * SNOVA_v * SNOVA_lr + 16];
        memset(whipped_sig, 0, sizeof(whipped_sig));
        c->sum_t1 = sum_t1; c->whipped_sig = whipped_sig;
#endif
        c->Fvv = Fvv_in_GF16Matrix;

        RCT_PT(_pa);
#if RCT_SIGN_JOG
        rct_sign_whipbuild_jog(c);
#elif RCT_OQDF
        rct_sign_whipbuild_oqdf(c);
#elif RCT_USE_SIMD && SNOVA_r <= 16
        rct_sign_whipbuild_q16(c);
#elif RCT_OQWV
        rct_sign_whipbuild_oqwv(c);
#else
        rct_sign_whipbuild_scalar(c);
#endif

#if RCT_SIGN_JOG
        rct_sign_sumt_jog(c);
#elif RCT_USE_SIMD
        rct_sign_s1whip_simd(c);
#elif RCT_Q_SIMD && (SNOVA_r == SNOVA_l)
        rct_sign_sumt_oddq_sq(c);
#elif RCT_Q_SIMD && RCT_Q_HAVE_MAGIC
        rct_sign_sumt_oddq_rect(c);
#else
        rct_sign_sumt_scalar(c);
#endif
        RCT_PT(_pb); RCT_PACC(4,_pa,_pb); RCT_PT(_pa);

#if defined(SNOVA_CT_CANARY) && (SNOVA_CT_CANARY + 0) == 3
#if !RCT_SIGN_JOG
#error "SNOVA_CT_CANARY=3 needs the RCT_SIGN_JOG arm (l=5 shapes); pick an l=5 parameter set"
#endif
        { volatile int ct_canary3_sink = 0;
          __m128i cb = rct_sj_bc128_pub(c->sum_t1[0]);
          ct_canary3_sink ^= _mm_cvtsi128_si32(cb); (void)ct_canary3_sink; }
#endif
#if RCT_SIGN_JOG
        rct_sign_fvv_jog(c);
#elif RCT_Q_SIMD && (SNOVA_r == SNOVA_l) && RCT_Q_HAVE_MAGIC && !defined(RCT_FVV_SCALAR)
        rct_sign_fvv_oddq_sq(c);
#elif RCT_Q_SIMD && (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC
        rct_sign_fvv_oddq_rect(c);
#elif defined(RCT_CMS2_ONLY)
        rct_sign_fvv_cms2(c);
#else
        rct_sign_fvv_std(c);
#endif

        RCT_PT(_pb); RCT_PACC(5,_pa,_pb); RCT_PT(_pa);
        for (int mi = 0; mi < SNOVA_o; ++mi)
            for (int i1 = 0; i1 < SNOVA_lr; i1++)
                gauss[mi * SNOVA_lr + i1][SNOVA_o * SNOVA_lr] =
                    gf_sub(hash_in_GF16[mi * SNOVA_lr + i1], Fvv_in_GF16Matrix[mi * SNOVA_lr + i1]);

#ifndef RCT_CM_ACTIVE
        RCT_SCRATCH _Alignas(32) gf_t whipped_F21[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr + 32];
        RCT_SCRATCH _Alignas(32) gf_t whipped_F12[SNOVA_m1 * SNOVA_l * SNOVA_o * SNOVA_lr + 32];
        memset(whipped_F21, 0, sizeof(whipped_F21));
        memset(whipped_F12, 0, sizeof(whipped_F12));
        c->whipped_F21 = whipped_F21; c->whipped_F12 = whipped_F12;
#endif
#ifdef RCT_PROFILE
        uint64_t _s0 = __rdtsc(), _s1, _s2;
        RCT_PACC6(0, _pa, _s0);
#endif

#if RCT_SIGN_JOG
        rct_sign_wf_jog(c);
#elif defined(RCT_CM_ACTIVE) || defined(RCT_CMS3_ONLY)
        rct_sign_wF_gauss_cm(c);
#elif RCT_Q_SIMD && (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC
        rct_sign_wF_gauss_oddq_rect(c);
#elif RCT_Q_SIMD && (SNOVA_r == SNOVA_l)
        rct_sign_wf_oddq_sq(c);
#else
        rct_sign_wf_std(c);
#endif
#ifdef RCT_PROFILE
        _s1 = __rdtsc();
        RCT_PACC6(1, _s0, _s1);
#endif

#if defined(RCT_CM_ACTIVE) || defined(RCT_CMS3_ONLY)
#elif RCT_Q_SIMD && (SNOVA_r != SNOVA_l) && RCT_Q_HAVE_MAGIC
#elif RCT_Q_SIMD && (SNOVA_r == SNOVA_l)
        rct_sign_gauss_scatter_oddq_sq(c);
#else
        rct_sign_gauss_scatter_std(c);
#endif
#ifdef RCT_PROFILE
        _s2 = __rdtsc();
        RCT_PACC6(2, _s1, _s2);
#endif

        RCT_PT(_pb); RCT_PACC(6,_pa,_pb); RCT_PT(_pa);
#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR)
#if RCT_Q_HAVE_MAGIC
        flag_redo = rct_sign_gauss_oddq_magic(c, gu, sol16, solution);
#else
        flag_redo = rct_sign_gauss_oddq_nomagic(c);
#endif
#else
        flag_redo = rct_sign_gauss_q16(c);
#endif

        if (!flag_redo) {
#if SNOVA_Q == 16
#if defined(RCT_GAUSS_SCALAR)
#define RCT_SOLPAD_HAS_GS 1
#else
#define RCT_SOLPAD_HAS_GS 0
#endif
#if defined(RCT_T12_MULLO)
#define RCT_SOLPAD_HAS_TM 1
#else
#define RCT_SOLPAD_HAS_TM 0
#endif
#define RCT_SOLPAD_BACKSUB (RCT_GFMUL_ANY && !RCT_SOLPAD_HAS_GS)
#define RCT_SOLPAD_T12 ( \
    (RCT_HOT_QRP16 && !RCT_USE_GFNI && RCT_SOLPAD_HAS_TM && (RCT_T12_MULLO + 0) && !RCT_SOLPAD_HAS_GS && (SNOVA_L == 4)) || \
    ((RCT_USE_GFNI || RCT_HOT_QRP16) && !RCT_SOLPAD_HAS_GS && (SNOVA_L == 4)) || \
    (RCT_SIGN_JOG && RCT_GFMUL_ANY && !RCT_SOLPAD_HAS_GS) )
#if (RCT_SOLPAD_BACKSUB) != (RCT_SOLPAD_T12)
#error "solpad interlock broken: backsub SIMD selection != apply_t12 solpad selection (RCT_SOLPAD_BACKSUB/RCT_SOLPAD_T12 must be flipped together)"
#endif
#undef RCT_SOLPAD_T12
#endif
#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC
#elif SNOVA_Q == 16 && RCT_SOLPAD_BACKSUB
            rct_sign_backsub_gfni(c, solpad, solution);
#else
            rct_sign_backsub_scalar(c, solution);
#endif
#if SNOVA_Q == 16
#undef RCT_SOLPAD_BACKSUB
#endif
            memcpy(signature_in_GF + SNOVA_v * SNOVA_lr, solution, SNOVA_o * SNOVA_lr);
#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC
            rct_sign_apply_t12_oddq(c, sol16);
#elif SNOVA_Q == 16 && RCT_HOT_QRP16 && !RCT_USE_GFNI && defined(RCT_T12_MULLO) \
    && (RCT_T12_MULLO + 0) && !defined(RCT_GAUSS_SCALAR) && SNOVA_L == 4
            rct_sign_apply_t12_mullo(c, solpad);
#elif SNOVA_Q == 16 && (RCT_USE_GFNI || RCT_HOT_QRP16) && !defined(RCT_GAUSS_SCALAR) && SNOVA_L == 4
            rct_sign_apply_t12_gfni(c, solpad);
#elif RCT_SIGN_JOG && SNOVA_Q == 16 && RCT_GFMUL_ANY && !defined(RCT_GAUSS_SCALAR)
            rct_sign_apply_t12_jog(c, solpad);
#else
            rct_sign_apply_t12_scalar(c, solution);
#endif
        }
        RCT_PT(_pb); RCT_PACC(7,_pa,_pb);
        SNOVA_CLEAR_OBJ(vinegar_in_byte);
#if !RCT_OQDF
        SNOVA_CLEAR_OBJ(sum_t1);
        SNOVA_CLEAR_OBJ(whipped_sig);
#endif
#ifndef RCT_CM_ACTIVE
        SNOVA_CLEAR_OBJ(whipped_F21);
        SNOVA_CLEAR_OBJ(whipped_F12);
#endif
    } while (flag_redo);
#ifdef RCT_PROFILE
    fprintf(stderr,"[SG] s0=%lu s1=%lu s2=%lu s3=%lu s4=%lu\n",
        (unsigned long)rct_pf[3],(unsigned long)rct_pf[4],(unsigned long)rct_pf[5],
        (unsigned long)rct_pf[6],(unsigned long)rct_pf[7]);
    fprintf(stderr,"[SG6] s0=%lu s1=%lu s2=%lu\n",
        (unsigned long)rct_pf6[0],(unsigned long)rct_pf6[1],(unsigned long)rct_pf6[2]);
    fprintf(stderr,"[SG6b] s0=%lu s1=%lu s2=%lu\n",
        (unsigned long)rct_pf6b[0],(unsigned long)rct_pf6b[1],(unsigned long)rct_pf6b[2]);
    rct_pf6b[0]=rct_pf6b[1]=rct_pf6b[2]=0;
    rct_pf6[0]=rct_pf6[1]=rct_pf6[2]=0;
    rct_pf[3]=rct_pf[4]=rct_pf[5]=rct_pf[6]=rct_pf[7]=0;
#endif

    compress_gf(sig, signature_in_GF, SNOVA_n * SNOVA_lr);
    memcpy(sig + BYTES_SIGNATURE - BYTES_SALT, salt, BYTES_SALT);
sign_cleanup:
    SNOVA_CLEAR_OBJ(signature_in_GF);
    SNOVA_CLEAR_OBJ(gauss);
    SNOVA_CLEAR_OBJ(solution);
#if RCT_Q_SIMD && !defined(RCT_GAUSS_SCALAR) && RCT_Q_HAVE_MAGIC
    SNOVA_CLEAR_OBJ(gu);
    SNOVA_CLEAR_OBJ(sol16);
#endif
#if SNOVA_Q == 16 && RCT_GFMUL_ANY && !defined(RCT_GAUSS_SCALAR)
    SNOVA_CLEAR_OBJ(solpad);
#endif
#if defined(RCT_CM_ACTIVE) || defined(RCT_CML_ONLY) || defined(RCT_CMS3_ONLY)
    SNOVA_CLEAR_OBJ(rct_cm_whip);
#if RCT_USE_GFNI
    SNOVA_CLEAR_OBJ(rct_cm_whipb);
#endif
#endif
#if RCT_OQDF
    SNOVA_CLEAR_OBJ(rct_oq_whip_w);
    SNOVA_CLEAR_OBJ(rct_oq_sum_t1u);
    SNOVA_CLEAR_OBJ(rct_oq_F21u);
    SNOVA_CLEAR_OBJ(rct_oq_F12u);
#endif
#if RCT_SIGN_STREAM
    free(P11);
#endif
    return sign_rc;
}

static int rct_sign(const uint8_t *sk, uint8_t *sig, const uint8_t *digest,
                    const size_t len_digest, const uint8_t *salt) {
#if SNOVA_WRAPPER_STACK
    _Alignas(64) rct_skx_t skx;
#else
    static rct_skx_t skx;
#endif
    rct_sk_expand(sk, &skx);
    int rc = rct_sign_expanded(&skx, sig, digest, len_digest, salt);
    SNOVA_CLEAR_OBJ(skx);
    return rc;
}

#endif
