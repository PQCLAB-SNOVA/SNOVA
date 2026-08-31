#ifndef SNOVA_RECT_VERIFY_H
#define SNOVA_RECT_VERIFY_H

typedef struct {
#if RCT_VERIFY_STREAM
    gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
#elif RCT_JOG_PKXJOG
    gf_t P[SNOVA_m1 * RCT_JOG_NL * RCT_JOG_L32];
#else
    gf_t P[SNOVA_m1 * SNOVA_n * SNOVA_n * SNOVA_l2];
#endif
    gf_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
    gf_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
    gf_t q1[SNOVA_o * SNOVA_alpha * SNOVA_l];
    gf_t q2[SNOVA_o * SNOVA_alpha * SNOVA_l];
    uint8_t pk_seed[SEED_LENGTH_PUBLIC];
#if HASH_PK
    uint8_t pk_hash[BYTES_PK_HASH];
#endif
} rct_pk_t;

#if FIXED_ABQ && ((SNOVA_q != 16) || ((SNOVA_l2 % 2) == 0)) && !RCT_JOG_PKXJOG
#define RCT_PKX_FUSED 1
#else
#define RCT_PKX_FUSED 0
#endif

typedef struct rct_pkx_ctx {
    rct_pk_t *pkx;
    const uint8_t *pk;
} rct_pkx_ctx;

typedef struct rct_vf_ctx {
    const rct_pk_t *pkx;
    const uint8_t *sig;
    gf_t *sig_gf;
    gf_t *hash_gf;
    gf_t *sum_t1;
    gf_t *sum_t1q;
    uint16_t *sum_t1s;
    uint8_t *whipped_sig2;
    uint8_t *sum_t1p;
#ifdef RCT_VPROF
    uint64_t *vpf; int *vpc; uint64_t va, vb;
#endif
} rct_vf_ctx;

#include "platforms/generic/rct_pkx.h"
#include "platforms/x86_avx2/rct_verify_avx2.h"
#include "platforms/x86_avx2/rct_verify_aq.h"
#include "platforms/x86_avx2/rct_verify_jog.h"
#include "platforms/x86_avx2/rct_verify_tile4.h"
#include "platforms/x86_avx2/rct_verify_oddq.h"
#include "platforms/generic/rct_verify_scalar.h"
#include "platforms/x86_avx2/rct_verify_emat.h"

static int rct_pk_expand(rct_pk_t *pkx, const uint8_t *pk) {
    rct_init();
    rct_pkx_ctx C = { .pkx = pkx, .pk = pk };
#if RCT_VF_TILE4 && RCT_TILE4_PREBAKE
    {
        const int prc = rct_pkx_expand(&C);
        if (prc == 0) rct_t4_prebake(pkx);
        return prc;
    }
#else
    return rct_pkx_expand(&C);
#endif
}

static int rct_verify(const rct_pk_t *pkx, const uint8_t *sig, const uint8_t *digest, const size_t len_digest) {
    rct_init();
#ifdef RCT_VPROF
    static uint64_t vpf[8]; static int vpc = 0; uint64_t va, vb;
#endif
    RCT_SCRATCH _Alignas(32) gf_t signature_in_GF[NUMGF_SIGNATURE + 32];
    VPT(va);
#if RCT_VF_MTK2 && (NUMGF_SIGNATURE % 2 == 0)
    rct_vf_expand_sig(signature_in_GF, sig, NUMGF_SIGNATURE);
#else
    if (expand_gf(signature_in_GF, sig, NUMGF_SIGNATURE)) return -1;
#endif
    VPT(vb); VPA(0, va, vb);

    gf_t hash_in_GF[SNOVA_o * SNOVA_lr + 32] = {0};
#if !RCT_VF_EMM
    RCT_SCRATCH _Alignas(32) gf_t sum_t1[SNOVA_m1 * SNOVA_l2 * SNOVA_r2 + 64];
    memset(sum_t1, 0, sizeof(sum_t1));
#else
    RCT_SCRATCH _Alignas(64) gf_t sum_t1q[SNOVA_m1 * SNOVA_l2 * 64];
#endif
#if RCT_Q_SIMD && !RCT_Q_EMM
    RCT_SCRATCH _Alignas(32) uint16_t sum_t1s[SNOVA_m1 * SNOVA_l2 * SNOVA_r2];
#endif

    rct_vf_ctx C;
    rct_vf_ctx *c = &C;
    c->pkx = pkx;
    c->sig = sig;
    c->sig_gf = signature_in_GF;
    c->hash_gf = hash_in_GF;
#if !RCT_VF_EMM
    c->sum_t1 = sum_t1;
#else
    c->sum_t1q = sum_t1q;
#endif
#if RCT_Q_SIMD && !RCT_Q_EMM
    c->sum_t1s = sum_t1s;
#endif
#ifdef RCT_VPROF
    c->vpf = vpf; c->vpc = &vpc;
#endif

#if RCT_VF_JOG
    {
#if RCT_USE_SIMD
        RCT_SCRATCH _Alignas(32) uint8_t jog_whipped_sig2[SNOVA_l * SNOVA_n * SNOVA_lr32];
        memset(jog_whipped_sig2, 0, sizeof(jog_whipped_sig2));
        c->whipped_sig2 = jog_whipped_sig2;
        VPT(va);
        rct_vf_whip_q16(c);
        VPT(vb); VPA(1, va, vb);
#endif
        VPT(va);
#if RCT_VF_TILE4
        rct_vf_tile4(c);
#else
        rct_vf_jog(c);
#endif
        VPT(vb); VPA(2, va, vb);
    }
#elif RCT_USE_SIMD
    {
        RCT_SCRATCH _Alignas(32) uint8_t whipped_sig2[SNOVA_l * SNOVA_n * SNOVA_lr32];
#if !(RCT_USE_GFNI && RCT_VF_MTK2)
        memset(whipped_sig2, 0, sizeof(whipped_sig2));
#endif
        c->whipped_sig2 = whipped_sig2;
        VPT(va);
        rct_vf_whip_q16(c);
        VPT(vb); VPA(1, va, vb);
        RCT_SCRATCH _Alignas(32) uint8_t sum_t1p[(SNOVA_m1 * SNOVA_l * SNOVA_r
            + ((SNOVA_r) < 7 ? (8 - (SNOVA_r)) : 1)) * SNOVA_lr32];
        memset(sum_t1p, 0, sizeof(sum_t1p));
        c->sum_t1p = sum_t1p;
        VPT(va);
#if RCT_VF_AQ
        rct_vf_contract_aq(c);
#else
        rct_vf_contract_nonaq(c);
#endif
        VPT(vb); VPA(2, va, vb); VPT(va);
        rct_vf_reindex_q16(c);
        VPT(vb); VPA(3, va, vb);
    }
#else
    memset(sum_t1, 0, sizeof(sum_t1));
#if RCT_Q_SIMD
#if (SNOVA_r == SNOVA_l)
    _Static_assert((uint32_t)SNOVA_n * SNOVA_l * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "u16 accumulation overflow guard");
#else
#endif
        {
            VPT(va); rct_vf_whip_oddq(c); VPT(vb); VPA(1, va, vb);
            VPT(va); rct_vf_rl_oddq(c); VPT(vb); VPA(2, va, vb);
            VPT(va); rct_vf_reindex_oddq(c); VPT(vb); VPA(3, va, vb);
        }
#else
    rct_vf_contract_ref(c);
#endif
#endif

    VPT(va);
    rct_vf_emat(c);

    VPT(vb); VPA(4, va, vb); VPT(va);
    uint8_t signed_bytes[BYTES_HASH];
    uint8_t signed_gf[GF16_HASH] = {0};
    const uint8_t *salt = sig + BYTES_SIGNATURE - BYTES_SALT;
#if HASH_PK
    hash_combined(signed_bytes, digest, len_digest, pkx->pk_hash, salt);
#else
    hash_combined(signed_bytes, digest, len_digest, pkx->pk_seed, salt);
#endif
    expand_gf(signed_gf, signed_bytes, GF16_HASH);

    int result = 0;
    for (int i = 0; i < GF16_HASH; ++i)
        if (hash_in_GF[i] != signed_gf[i]) { result = -1; break; }
    VPT(vb); VPA(5, va, vb);
#ifdef RCT_VPROF
    if (++vpc == VPROF_AT) {
        fprintf(stderr, "[VF-PROF n=%d] s0=%lu s1=%lu s2=%lu s3=%lu s4=%lu s5=%lu  (per-call avg)\n",
                VPROF_AT,
                (unsigned long)(vpf[0] / VPROF_AT), (unsigned long)(vpf[1] / VPROF_AT),
                (unsigned long)(vpf[2] / VPROF_AT), (unsigned long)(vpf[3] / VPROF_AT),
                (unsigned long)(vpf[4] / VPROF_AT), (unsigned long)(vpf[5] / VPROF_AT));
    }
#endif
    return result;
}

#endif
