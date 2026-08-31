/**
 * @file ct_test.c
 */
#define _POSIX_C_SOURCE 200809L
#include "../src/snova.h"
#include "../src/snova_core/drbg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/snova_core/ct_poison.h"

int main(int argc, char **argv) {
    int n = (argc > 1) ? atoi(argv[1]) : 2;
    if (n < 1) n = 2;

    uint8_t entropy[48];
    for (int i = 0; i < 48; ++i) entropy[i] = (uint8_t)i;
    randombytes_init(entropy, NULL, 0);

    uint8_t *pk = (uint8_t *)malloc(1 << 20);
    uint8_t *sk = (uint8_t *)malloc(1 << 20);
    if (!pk || !sk) { fprintf(stderr, "alloc failed\n"); return 1; }

    uint8_t msg[64];
    for (int i = 0; i < 64; ++i) msg[i] = (uint8_t)(i ^ 0xA5);
    size_t mlen = 64;
    uint8_t *sm = (uint8_t *)malloc(mlen + SNOVA_SM_OVERHEAD);
    uint8_t *mout = (uint8_t *)malloc(mlen);
    if (!sm || !mout) { fprintf(stderr, "alloc failed\n"); return 1; }

    int bad = 0;
#if defined(SNOVA_CT_CANARY) && (SNOVA_CT_CANARY + 0) >= 1
    volatile int canary_sink = 0;
#endif
    for (int i = 0; i < n; ++i) {
        int rc = snova_crypto_sign_keypair(pk, sk);
        SNOVA_CT_DECLASSIFY(&rc, sizeof rc);
        SNOVA_CT_DECLASSIFY(pk, SNOVA_PK_BYTES);
        if (rc) { bad++; continue; }

#if defined(SNOVA_CT_CANARY) && (SNOVA_CT_CANARY + 0) == 1
        if (sk[SNOVA_SEED_PUB] & 1u) canary_sink ^= 1;
#endif

        size_t smlen = 0;
        rc = snova_crypto_sign(sm, &smlen, msg, mlen, sk);
        SNOVA_CT_DECLASSIFY(&rc, sizeof rc);
        SNOVA_CT_DECLASSIFY(&smlen, sizeof smlen);
        if (rc) { bad++; continue; }
#if defined(SNOVA_CT_CANARY) && (SNOVA_CT_CANARY + 0) == 2
        if (sm[0] & 1u) canary_sink ^= 2;
#endif
        SNOVA_CT_DECLASSIFY(sm, smlen);

        size_t mlen2 = 0;
        rc = snova_crypto_sign_open(mout, &mlen2, sm, smlen, pk);
        SNOVA_CT_DECLASSIFY(&rc, sizeof rc);
        SNOVA_CT_DECLASSIFY(&mlen2, sizeof mlen2);
        if (rc != 0 || mlen2 != mlen) { bad++; continue; }
        SNOVA_CT_DECLASSIFY(mout, mlen2);
        if (memcmp(mout, msg, mlen) != 0) bad++;
    }

    {
        uint8_t *skx = (uint8_t *)malloc(snova_skx_bytes());
        uint8_t *pkx = (uint8_t *)malloc(snova_pkx_bytes());
        uint8_t *sig = (uint8_t *)malloc(SNOVA_SM_OVERHEAD);
        if (!skx || !pkx || !sig) { fprintf(stderr, "alloc failed\n"); return 1; }
        uint8_t digest[64], salt[16];
        for (int i = 0; i < 64; ++i) digest[i] = (uint8_t)(i * 13 + 7);
        for (int i = 0; i < 16; ++i) salt[i] = (uint8_t)(0x5A ^ i);

        snova_sk_expand_skx(skx, sk);
        if (snova_pk_expand(pkx, pk) != 0) bad++;
        for (int k = 0; k < 2; ++k) {
            salt[0] = (uint8_t)(0x5A + k);
            (void)snova_sign_digest_skx(sig, digest, 64, salt, skx);
            SNOVA_CT_DECLASSIFY(sig, SNOVA_SM_OVERHEAD);
            int rc = snova_verify_digest_pkx(digest, 64, sig, pkx);
            SNOVA_CT_DECLASSIFY(&rc, sizeof rc);
            if (rc != 0) bad++;
        }
        snova_skx_clear(skx);
        free(skx); free(pkx); free(sig);
    }

    printf("ct_test: engine=%s\n",
#if defined(SNOVA_CT_MSAN)
           "MSan"
#else
           "valgrind"
#endif
          );
    printf("ct_test: %d iters, functional bad=%d (CT verdict = memcheck errors)\n",
           n, bad);
    free(pk); free(sk); free(sm); free(mout);
    return bad ? 2 : 0;
}
