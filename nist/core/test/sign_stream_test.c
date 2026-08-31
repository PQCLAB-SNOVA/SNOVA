/**
 * @file test/sign_stream_test.c
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/api.h"
#include "../src/snova_core/drbg.h"

#define N_KEYS 5
#define N_MSGS 4
#define MSG_MAX 300

static int g_fail = 0;
static void expect(int c, const char *w, int k, int j) {
    if (!c) { printf("  FAIL: %s (key %d msg %d)\n", w, k, j); g_fail = 1; }
}

int main(void) {
    uint8_t entropy[48];
    for (int i = 0; i < 48; ++i) entropy[i] = (uint8_t)(i * 11 + 5);
    randombytes_init(entropy, NULL, 256);

    static uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    static uint8_t sk[CRYPTO_SECRETKEYBYTES];
    static uint8_t m[MSG_MAX];
    static uint8_t sm[CRYPTO_BYTES + MSG_MAX];
    static uint8_t mout[CRYPTO_BYTES + MSG_MAX];

    uint64_t digest = 1469598103934665603ULL;
    size_t total = 0;

    for (int k = 0; k < N_KEYS; ++k) {
        if (crypto_sign_keypair(pk, sk) != 0) { printf("  FAIL keypair %d\n", k); return 1; }
        for (int j = 0; j < N_MSGS; ++j) {
            unsigned long long mlen = (unsigned long long)(17 + 61 * j + 7 * k);
            randombytes(m, (size_t)mlen);
            unsigned long long smlen = 0, mlen2 = 0;
            if (crypto_sign(sm, &smlen, m, mlen, sk) != 0) { printf("  FAIL sign %d,%d\n", k, j); return 1; }

            for (unsigned long long b = 0; b < smlen; ++b) {
                digest ^= sm[b]; digest *= 1099511628211ULL;
            }
            total += (size_t)smlen;

            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) == 0 &&
                       mlen2 == mlen && memcmp(mout, m, (size_t)mlen) == 0,
                   "valid signature must verify", k, j);
            sm[(j + 1) % 16] ^= 0x40;
            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) != 0,
                   "corrupted signature must reject", k, j);
            sm[(j + 1) % 16] ^= 0x40;
        }
    }

    printf("sign_stream_test: %s  keys=%d msgs=%d  sm_total=%zu  sig_digest=%016llx\n",
#if defined(SNOVA_SIGN_STREAM) && SNOVA_SIGN_STREAM
           "STREAM(sink D)",
#else
           "materialized",
#endif
           N_KEYS, N_MSGS, total, (unsigned long long)digest);
    if (g_fail) { printf("== SIGN-STREAM-TEST FAIL ==\n"); return 1; }
    printf("== SIGN-STREAM-TEST PASS ==\n");
    return 0;
}
