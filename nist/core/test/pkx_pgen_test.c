/**
 * @file test/pkx_pgen_test.c
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/api.h"
#include "../src/snova_core/drbg.h"

extern size_t snova_pkx_bytes(void);
extern int snova_pk_expand(uint8_t *pkx, const uint8_t *pk);

#define N_KEYS 4

int main(void) {
    uint8_t entropy[48];
    for (int i = 0; i < 48; ++i) entropy[i] = (uint8_t)(i * 13 + 1);
    randombytes_init(entropy, NULL, 256);

    static uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    static uint8_t sk[CRYPTO_SECRETKEYBYTES];
    size_t pkx_bytes = snova_pkx_bytes();
    uint8_t *pkx = (uint8_t *)malloc(pkx_bytes);
    if (!pkx) { printf("== PKX-PGEN-TEST FAIL (oom) ==\n"); return 1; }

    uint64_t digest = 1469598103934665603ULL;

    for (int k = 0; k < N_KEYS; ++k) {
        if (crypto_sign_keypair(pk, sk) != 0) { printf("  FAIL keypair %d\n", k); return 1; }
        memset(pkx, 0, pkx_bytes);
        if (snova_pk_expand(pkx, pk) != 0) { printf("  FAIL pk_expand %d\n", k); return 1; }
        for (size_t b = 0; b < pkx_bytes; ++b) { digest ^= pkx[b]; digest *= 1099511628211ULL; }
    }
    free(pkx);

    printf("pkx_pgen_test: %s  keys=%d  pkx_bytes=%zu  pkx_digest=%016llx\n",
#if defined(SNOVA_PKX_PGEN) && SNOVA_PKX_PGEN
           "PGEN(sink C)",
#else
           "fused-unpack",
#endif
           N_KEYS, pkx_bytes, (unsigned long long)digest);
    printf("== PKX-PGEN-TEST PASS ==\n");
    return 0;
}
