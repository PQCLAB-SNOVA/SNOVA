/**
 * @file memtest.c
 */
#define _DEFAULT_SOURCE 1
#include "../src/snova.h"
#include "../src/snova_core/drbg.h"
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef BYTES_DIGEST
#define BYTES_DIGEST 64
#endif

int main(int argc, char **argv) {
    const char *api = (argc > 1) ? argv[1] : "sign";
    int reps = (argc > 2) ? atoi(argv[2]) : 1;
    if (reps < 1) reps = 1;

    uint8_t entropy[48];
    for (int i = 0; i < 48; ++i) entropy[i] = (uint8_t)i;
    randombytes_init(entropy, NULL, 0);

    uint8_t *pk = (uint8_t *)malloc(1 << 20);
    uint8_t *sk = (uint8_t *)malloc(1 << 20);
    uint8_t *sm = (uint8_t *)malloc(SNOVA_SM_OVERHEAD + 64);
    uint8_t *skx = (uint8_t *)malloc(snova_skx_bytes());
    uint8_t *pkx = (uint8_t *)malloc(snova_pkx_bytes());
    if (!pk || !sk || !sm || !skx || !pkx) { fprintf(stderr, "alloc failed\n"); return 1; }
    memset(skx, 0, snova_skx_bytes());
    memset(pkx, 0, snova_pkx_bytes());

    uint8_t bseed[48];
    uint8_t salt[16];
    uint8_t digest[BYTES_DIGEST];
    for (int i = 0; i < 48; ++i) bseed[i] = (uint8_t)(i + 1);
    for (int i = 0; i < 16; ++i) salt[i] = (uint8_t)(i * 3 + 1);
    memset(digest, 0, BYTES_DIGEST);

    int rc = snova_keygen(pk, sk, bseed);

    if (!strcmp(api, "keygen")) {
        for (int i = 0; i < reps; ++i) rc |= snova_keygen(pk, sk, bseed);
    } else if (!strcmp(api, "sk_expand")) {
        for (int i = 0; i < reps; ++i) snova_sk_expand_skx(skx, sk);
    } else if (!strcmp(api, "sign")) {
        snova_sk_expand_skx(skx, sk);
        for (int i = 0; i < reps; ++i) rc |= snova_sign_digest_skx(sm, digest, BYTES_DIGEST, salt, skx);
    } else if (!strcmp(api, "pk_expand")) {
        for (int i = 0; i < reps; ++i) rc |= snova_pk_expand(pkx, pk);
    } else if (!strcmp(api, "verify")) {
        snova_sk_expand_skx(skx, sk);
        rc |= snova_sign_digest_skx(sm, digest, BYTES_DIGEST, salt, skx);
        rc |= snova_pk_expand(pkx, pk);
        for (int i = 0; i < reps; ++i) rc |= snova_verify_digest_pkx(digest, BYTES_DIGEST, sm, pkx);
    } else {
        fprintf(stderr, "unknown api '%s' (keygen|sk_expand|sign|pk_expand|verify)\n", api);
        return 2;
    }

    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    printf("MAXRSS_KB=%ld\n", (long)ru.ru_maxrss);
    if (rc != 0) fprintf(stderr, "WARN: op returned rc=%d\n", rc);

    free(skx); free(pkx); free(sm); free(pk); free(sk);
    return 0;
}
