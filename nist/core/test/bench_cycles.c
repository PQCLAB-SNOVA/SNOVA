/**
 * @file bench_cycles.c
 */
#define _POSIX_C_SOURCE 200809L
#include "../src/snova.h"
#include "../src/snova_core/drbg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

static inline uint64_t rdtsc_now(void) {
    return __rdtsc();
}

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv) {
    int n = (argc > 1) ? atoi(argv[1]) : 512;
    if (n < 1) n = 512;

    uint8_t entropy[48];
    for (int i = 0; i < 48; ++i) entropy[i] = (uint8_t)i;
    randombytes_init(entropy, NULL, 0);

    uint8_t *pk = (uint8_t *)malloc(1 << 20);
    uint8_t *sk = (uint8_t *)malloc(1 << 20);
    if (!pk || !sk) { fprintf(stderr, "alloc failed\n"); return 1; }

    uint8_t msg[64];
    for (int i = 0; i < 64; ++i) msg[i] = (uint8_t)(i ^ 0xA5);
    size_t mlen = 64;
    size_t sm_buf = mlen + SNOVA_SM_OVERHEAD;
    uint8_t *sm = (uint8_t *)malloc(sm_buf);
    uint8_t *mout = (uint8_t *)malloc(mlen);
    if (!sm || !mout) { fprintf(stderr, "alloc failed\n"); return 1; }

    uint64_t *c_kg = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *c_sg = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *c_vf = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));

    int bad = 0;
    for (int i = 0; i < n; ++i) {
        uint64_t t0 = rdtsc_now();
        snova_crypto_sign_keypair(pk, sk);
        uint64_t t1 = rdtsc_now();
        size_t smlen = 0;
        snova_crypto_sign(sm, &smlen, msg, mlen, sk);
        uint64_t t2 = rdtsc_now();
        size_t mlen2 = 0;
        int rc = snova_crypto_sign_open(mout, &mlen2, sm, smlen, pk);
        uint64_t t3 = rdtsc_now();
        if (rc != 0 || mlen2 != mlen || memcmp(mout, msg, mlen) != 0) bad++;
        c_kg[i] = t1 - t0;
        c_sg[i] = t2 - t1;
        c_vf[i] = t3 - t2;
    }

    qsort(c_kg, (size_t)n, sizeof(uint64_t), cmp_u64);
    qsort(c_sg, (size_t)n, sizeof(uint64_t), cmp_u64);
    qsort(c_vf, (size_t)n, sizeof(uint64_t), cmp_u64);

    printf("E2E rdtsc cycles (N=%d, NIST API):\n", n);
    printf("  keygen      median = %lu\n", (unsigned long)c_kg[n / 2]);
    printf("  sign        median = %lu\n", (unsigned long)c_sg[n / 2]);
    printf("  verify      median = %lu\n", (unsigned long)c_vf[n / 2]);
    if (bad) printf("  WARN: %d/%d verify failures\n", bad, n);

    free(c_kg); free(c_sg); free(c_vf);
    free(sm); free(mout); free(pk); free(sk);
    return bad ? 2 : 0;
}
