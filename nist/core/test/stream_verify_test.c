/**
 * @file test/stream_verify_test.c
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/api.h"
#include "../src/snova_core/drbg.h"

#define N_KEYS 4
#define N_MSGS 3
#define MSG_MAX 256

static int g_fail = 0;

static void expect(int cond, const char *what, int k, int j) {
    if (!cond) {
        printf("  FAIL: %s (key %d, msg %d)\n", what, k, j);
        g_fail = 1;
    }
}

int main(void) {
    uint8_t entropy[48];
    for (int i = 0; i < 48; ++i) entropy[i] = (uint8_t)(i * 7 + 3);
    randombytes_init(entropy, NULL, 256);

    static uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    static uint8_t sk[CRYPTO_SECRETKEYBYTES];
    static uint8_t m[MSG_MAX];
    static uint8_t sm[CRYPTO_BYTES + MSG_MAX];
    static uint8_t mout[CRYPTO_BYTES + MSG_MAX];

    printf("stream_verify_test: %s (pk=%d sk=%d sig=%d)\n",
#if defined(SNOVA_VERIFY_STREAM) && SNOVA_VERIFY_STREAM
           "STREAM build",
#else
           "default (non-stream) build",
#endif
           (int)CRYPTO_PUBLICKEYBYTES, (int)CRYPTO_SECRETKEYBYTES,
           (int)CRYPTO_BYTES);

    for (int k = 0; k < N_KEYS; ++k) {
        if (crypto_sign_keypair(pk, sk) != 0) {
            printf("  FAIL: keypair (key %d)\n", k);
            return 1;
        }
        for (int j = 0; j < N_MSGS; ++j) {
            unsigned long long mlen = (unsigned long long)(33 + 64 * j + k);
            randombytes(m, (size_t)mlen);
            unsigned long long smlen = 0, mlen2 = 0;
            if (crypto_sign(sm, &smlen, m, mlen, sk) != 0) {
                printf("  FAIL: sign (key %d, msg %d)\n", k, j);
                return 1;
            }

            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) == 0 &&
                       mlen2 == mlen && memcmp(mout, m, (size_t)mlen) == 0,
                   "valid signature must verify", k, j);

            sm[j % 8] ^= 0x10;
            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) != 0,
                   "corrupted signature must reject", k, j);
            sm[j % 8] ^= 0x10;

            sm[CRYPTO_BYTES - 1] ^= 0x01;
            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) != 0,
                   "corrupted salt must reject", k, j);
            sm[CRYPTO_BYTES - 1] ^= 0x01;

            sm[CRYPTO_BYTES + (size_t)(mlen / 2)] ^= 0x80;
            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) != 0,
                   "corrupted message must reject", k, j);
            sm[CRYPTO_BYTES + (size_t)(mlen / 2)] ^= 0x80;

            pk[CRYPTO_PUBLICKEYBYTES - 1 - j] ^= 0x0F;
            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) != 0,
                   "corrupted pk (P22 region) must reject", k, j);
            pk[CRYPTO_PUBLICKEYBYTES - 1 - j] ^= 0x0F;

            pk[j % 16] ^= 0x21;
            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) != 0,
                   "corrupted pkseed must reject", k, j);
            pk[j % 16] ^= 0x21;

            expect(crypto_sign_open(mout, &mlen2, sm, CRYPTO_BYTES - 1, pk) != 0,
                   "truncated sm must reject", k, j);

            expect(crypto_sign_open(mout, &mlen2, sm, smlen, pk) == 0,
                   "restored signature must verify again", k, j);
        }
    }

    if (g_fail) {
        printf("== STREAM-VERIFY-TEST FAIL ==\n");
        return 1;
    }
    printf("== STREAM-VERIFY-TEST PASS (%d keys x %d msgs, %s%s) ==\n",
           N_KEYS, N_MSGS,
#if defined(SNOVA_VERIFY_STREAM) && SNOVA_VERIFY_STREAM
           "stream",
#else
           "default",
#endif
#ifdef EVAL_STREAM_SELFTEST
           " + per-verify stream==materialized hash selftest"
#else
           ""
#endif
    );
    return 0;
}
