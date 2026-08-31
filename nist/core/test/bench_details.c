/**
 * @file bench_details.c
 */
#define _POSIX_C_SOURCE 200809L
#include "../src/snova.h"
#include "../src/snova_core/drbg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

#ifndef BYTES_DIGEST
#define BYTES_DIGEST 64
#endif

static inline uint64_t rdtsc_now(void) { return __rdtsc(); }

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static uint64_t median_of(uint64_t *v, int n) {
    qsort(v, (size_t)n, sizeof(uint64_t), cmp_u64);
    return v[n / 2];
}

int main(int argc, char **argv) {
    int n = (argc > 1) ? atoi(argv[1]) : 512;
    if (n < 1) n = 512;

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

    printf("SIZES pk=%lu sig=%lu sk=%lu pkx=%lu skx=%lu\n",
           (unsigned long)SNOVA_PK_BYTES,
           (unsigned long)SNOVA_SM_OVERHEAD,
           (unsigned long)(SNOVA_SEED_LEN + SNOVA_BYTES_PK_HASH),
           (unsigned long)snova_pkx_bytes(),
           (unsigned long)snova_skx_bytes());

    uint64_t *c_gk = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *c_se = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *c_sg = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *c_pe = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    uint64_t *c_vf = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    if (!c_gk || !c_se || !c_sg || !c_pe || !c_vf) { fprintf(stderr, "alloc failed\n"); return 1; }

    uint8_t bseed[48];
    uint8_t salt[16];
    uint8_t digest[BYTES_DIGEST];
    int bad = 0;

    for (int i = 0; i < n; ++i) {
        randombytes(bseed, 48);
        randombytes(salt, 16);
        memset(digest, 0, BYTES_DIGEST);
        digest[0] = (uint8_t)(i & 0xff);
        digest[1] = (uint8_t)((i >> 8) & 0xff);
        digest[2] = (uint8_t)((i >> 16) & 0xff);

        uint64_t t0 = rdtsc_now();
        int r = snova_keygen(pk, sk, bseed);
        uint64_t t1 = rdtsc_now();
        snova_sk_expand_skx(skx, sk);
        uint64_t t2 = rdtsc_now();
        r |= snova_sign_digest_skx(sm, digest, BYTES_DIGEST, salt, skx);
        uint64_t t3 = rdtsc_now();
        r |= snova_pk_expand(pkx, pk);
        uint64_t t4 = rdtsc_now();
        int rv = snova_verify_digest_pkx(digest, BYTES_DIGEST, sm, pkx);
        uint64_t t5 = rdtsc_now();

        if (r != 0 || rv != 0) bad++;
        c_gk[i] = t1 - t0;
        c_se[i] = t2 - t1;
        c_sg[i] = t3 - t2;
        c_pe[i] = t4 - t3;
        c_vf[i] = t5 - t4;
    }

    printf("DETAILS rdtsc cycles (N=%d):\n", n);
    printf("  genkeys   median = %lu\n", (unsigned long)median_of(c_gk, n));
    printf("  sk_expand median = %lu\n", (unsigned long)median_of(c_se, n));
    printf("  sign      median = %lu\n", (unsigned long)median_of(c_sg, n));
    printf("  pk_expand median = %lu\n", (unsigned long)median_of(c_pe, n));
    printf("  verify    median = %lu\n", (unsigned long)median_of(c_vf, n));
    if (bad) printf("  WARN: %d/%d sign/verify failures\n", bad, n);

    free(c_gk); free(c_se); free(c_sg); free(c_pe); free(c_vf);
    free(skx); free(pkx); free(sm); free(pk); free(sk);
    return bad ? 2 : 0;
}
