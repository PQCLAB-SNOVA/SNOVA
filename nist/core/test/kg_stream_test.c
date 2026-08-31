/**
 * @file test/kg_stream_test.c
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/api.h"
#include "../src/snova_core/drbg.h"

#define N_KEYS 8

int main(void) {
    uint8_t entropy[48];
    for (int i = 0; i < 48; ++i) entropy[i] = (uint8_t)(i * 7 + 3);
    randombytes_init(entropy, NULL, 256);

    static uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    static uint8_t sk[CRYPTO_SECRETKEYBYTES];
    static uint8_t m[128];
    static uint8_t sm[CRYPTO_BYTES + 128];
    static uint8_t mout[CRYPTO_BYTES + 128];

    uint64_t digest = 1469598103934665603ULL;
    int fail = 0;

    for (int k = 0; k < N_KEYS; ++k) {
        if (crypto_sign_keypair(pk, sk) != 0) { printf("  FAIL keypair %d\n", k); return 1; }
        for (size_t b = 0; b < CRYPTO_PUBLICKEYBYTES; ++b) { digest ^= pk[b]; digest *= 1099511628211ULL; }
        for (size_t b = 0; b < CRYPTO_SECRETKEYBYTES; ++b) { digest ^= sk[b]; digest *= 1099511628211ULL; }

        unsigned long long mlen = 33 + 5 * (unsigned long long)k, smlen = 0, mlen2 = 0;
        randombytes(m, (size_t)mlen);
        if (crypto_sign(sm, &smlen, m, mlen, sk) != 0 ||
            crypto_sign_open(mout, &mlen2, sm, smlen, pk) != 0 ||
            mlen2 != mlen || memcmp(mout, m, (size_t)mlen) != 0) {
            printf("  FAIL sign/verify roundtrip (key %d)\n", k);
            fail = 1;
        }
    }

    printf("kg_stream_test: %s  keys=%d  kg_digest=%016llx\n",
#if defined(SNOVA_KEYGEN_STREAM) && SNOVA_KEYGEN_STREAM
           "STREAM(sink E)",
#else
           "materialized",
#endif
           N_KEYS, (unsigned long long)digest);
    if (fail) { printf("== KG-STREAM-TEST FAIL ==\n"); return 1; }
    printf("== KG-STREAM-TEST PASS ==\n");
    return 0;
}
