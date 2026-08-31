/**
 * @file test/fuzz_verify.c
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/api.h"
#include "../src/snova_core/drbg.h"

extern size_t snova_skx_bytes(void);
extern int snova_keygen(uint8_t *pk, uint8_t *sk, const uint8_t *seed);
extern void snova_sk_expand_skx(uint8_t *skx, const uint8_t *sk);

static uint8_t *g_pk = NULL, *g_sk = NULL, *g_valid_sm = NULL;
static size_t   g_valid_smlen = 0;
static size_t   g_skx_bytes = 0;
static int      g_inited = 0;
static long     g_accept_nontrivial = 0;

static uint64_t rng_state = 0x123456789abcdef0ULL;
static inline uint64_t xrand(void) {
    uint64_t x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng_state = x;
    return x;
}
static inline uint8_t rbyte(void) { return (uint8_t)(xrand() & 0xFF); }
static void rfill(uint8_t *p, size_t n) { for (size_t i = 0; i < n; ++i) p[i] = rbyte(); }
static size_t rrange(size_t lo, size_t hi) {
    if (hi <= lo) return lo;
    return lo + (size_t)(xrand() % (hi - lo + 1));
}

static int call_open(const uint8_t *sm, size_t smlen, const uint8_t *pk) {
    size_t msgcap = (smlen > CRYPTO_BYTES) ? (smlen - CRYPTO_BYTES) : 1;
    uint8_t *m = (uint8_t *)malloc(msgcap);
    uint8_t *sm_exact = (uint8_t *)malloc(smlen ? smlen : 1);
    memcpy(sm_exact, sm, smlen);
    unsigned long long mlen = 0;
    int rc = crypto_sign_open(m, &mlen, sm_exact, (unsigned long long)smlen, pk);
    (void)mlen;
    free(m);
    free(sm_exact);
    return rc;
}

static void note_accept(int rc, const uint8_t *pk_use, const uint8_t *sm, size_t smlen) {
    if (rc != 0) return;
    int trivial = (memcmp(pk_use, g_pk, CRYPTO_PUBLICKEYBYTES) == 0) &&
                  (smlen == g_valid_smlen) &&
                  (smlen == 0 || memcmp(sm, g_valid_sm, smlen) == 0);
    if (!trivial) {
        g_accept_nontrivial++;
        fprintf(stderr, "FATAL: non-trivial accept — possible forgery/bounds bug! "
                        "(smlen=%zu pk=%s)\n", smlen, pk_use == g_pk ? "valid" : "malformed");
#ifdef FUZZ_LIBFUZZER
        abort();
#endif
    }
}

static void harness_init(void) {
    if (g_inited) return;
    g_inited = 1;

    snova_init();

    uint8_t entropy[48];
    for (int i = 0; i < 48; ++i) entropy[i] = (uint8_t)i;
    randombytes_init(entropy, NULL, 256);

    g_pk = (uint8_t *)malloc(CRYPTO_PUBLICKEYBYTES);
    g_sk = (uint8_t *)malloc(CRYPTO_SECRETKEYBYTES);
    if (crypto_sign_keypair(g_pk, g_sk) != 0) { fprintf(stderr, "keypair fail\n"); exit(2); }

    const size_t baselen = 33;
    uint8_t basemsg[33];
    rfill(basemsg, baselen);
    size_t sm_cap = CRYPTO_BYTES + baselen;
    g_valid_sm = (uint8_t *)malloc(sm_cap);
    unsigned long long smlen = 0;
    if (crypto_sign(g_valid_sm, &smlen, basemsg, baselen, g_sk) != 0) {
        fprintf(stderr, "sign fail\n"); exit(2);
    }
    g_valid_smlen = (size_t)smlen;
    g_skx_bytes = snova_skx_bytes();

    if (call_open(g_valid_sm, g_valid_smlen, g_pk) != 0) {
        fprintf(stderr, "FATAL: valid signature rejected — harness/params mismatch\n");
        exit(2);
    }
}

static void harness_consume(const uint8_t *data, size_t size) {
    harness_init();

    uint8_t sel = size ? data[0] : 0;
    const uint8_t *payload = (size > 1) ? data + 1 : data;
    size_t plen = (size > 1) ? size - 1 : 0;

    if ((sel & 0xC0) == 0xC0) {
        uint8_t *fsk = (uint8_t *)malloc(CRYPTO_SECRETKEYBYTES);
        if (sel & 0x20) {
            memcpy(fsk, g_sk, CRYPTO_SECRETKEYBYTES);
            for (size_t i = 0; i < plen && i < (size_t)CRYPTO_SECRETKEYBYTES; ++i) fsk[i] ^= payload[i];
        } else {
            for (size_t i = 0; i < (size_t)CRYPTO_SECRETKEYBYTES; ++i) fsk[i] = plen ? payload[i % plen] : 0;
        }
        void *skx = malloc(g_skx_bytes);
        (void)snova_sk_expand_skx(skx, fsk);
        free(skx);
        free(fsk);
        return;
    }

    uint8_t *pk_own = NULL;
    const uint8_t *pk_use;
    switch (sel & 0x30) {
    case 0x00:
        pk_use = g_pk;
        break;
    case 0x10:
        pk_own = (uint8_t *)malloc(CRYPTO_PUBLICKEYBYTES);
        memcpy(pk_own, g_pk, CRYPTO_PUBLICKEYBYTES);
        for (size_t i = 0; i < plen && i < (size_t)CRYPTO_PUBLICKEYBYTES; ++i) pk_own[i] ^= payload[i];
        pk_use = pk_own;
        break;
    case 0x20:
        pk_own = (uint8_t *)malloc(CRYPTO_PUBLICKEYBYTES);
        for (size_t i = 0; i < (size_t)CRYPTO_PUBLICKEYBYTES; ++i) pk_own[i] = plen ? payload[i % plen] : 0xFF;
        pk_use = pk_own;
        break;
    default:
        pk_own = (uint8_t *)malloc(CRYPTO_PUBLICKEYBYTES);
        memset(pk_own, 0xFF, CRYPTO_PUBLICKEYBYTES);
        pk_use = pk_own;
        break;
    }

    if (sel & 0x08) {
        size_t smlen = g_valid_smlen;
        uint8_t *sm = (uint8_t *)malloc(smlen ? smlen : 1);
        memcpy(sm, g_valid_sm, smlen);
        for (size_t i = 0; i < plen && i < smlen; ++i) sm[i] ^= payload[i];
        int rc = call_open(sm, smlen, pk_use);
        note_accept(rc, pk_use, sm, smlen);
        free(sm);
    } else {
        int rc = call_open(payload, plen, pk_use);
        note_accept(rc, pk_use, payload, plen);
    }

    if (pk_own) free(pk_own);
}

#ifdef FUZZ_LIBFUZZER

int LLVMFuzzerInitialize(int *argc, char ***argv) {
    (void)argc; (void)argv;
    harness_init();
    return 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    harness_consume(data, size);
    return 0;
}

#else

static int replay_file(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "corpus: cannot open %s\n", path); return -1; }
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (n < 0) { fclose(f); return -1; }
    uint8_t *buf = (uint8_t *)malloc((size_t)n ? (size_t)n : 1);
    size_t rd = fread(buf, 1, (size_t)n, f);
    fclose(f);
    harness_consume(buf, rd);
    free(buf);
    return 0;
}

int main(int argc, char **argv) {
    long iters = (argc > 1) ? atol(argv[1]) : 100000;
    if (argc > 2) rng_state = strtoull(argv[2], NULL, 0) | 1ULL;

    harness_init();

    if (argc > 3) {
        int nrep = 0;
        for (int i = 3; i < argc; ++i) if (replay_file(argv[i]) == 0) nrep++;
        printf("fuzz_verify: replayed %d corpus file(s) through harness_consume\n", nrep);
    }

    {
        uint8_t *pk2 = (uint8_t *)malloc(CRYPTO_PUBLICKEYBYTES);
        memcpy(pk2, g_pk, CRYPTO_PUBLICKEYBYTES);
        pk2[0] ^= 1;
        if (call_open(g_valid_sm, g_valid_smlen, pk2) == 0)
            fprintf(stderr, "WARN: pk bit-flip still accepted (collision?)\n");
        free(pk2);
    }

    long accepted = 0, rejected = 0;
    long accept_by_mode[12] = {0};

    size_t max_smlen = g_valid_smlen + 64;

    for (long it = 0; it < iters; ++it) {
        int mode = (int)(xrand() % 12);
        size_t smlen;
        uint8_t *sm = NULL;
        const uint8_t *pk_use = g_pk;
        uint8_t *pk_own = NULL;

        switch (mode) {
        case 0:
            smlen = rrange(0, max_smlen);
            sm = (uint8_t *)malloc(smlen ? smlen : 1);
            rfill(sm, smlen);
            break;
        case 1:
            smlen = rrange(0, max_smlen);
            sm = (uint8_t *)malloc(smlen ? smlen : 1);
            memset(sm, 0x00, smlen);
            break;
        case 2:
            smlen = rrange(0, max_smlen);
            sm = (uint8_t *)malloc(smlen ? smlen : 1);
            memset(sm, 0xFF, smlen);
            break;
        case 3:
            smlen = g_valid_smlen;
            sm = (uint8_t *)malloc(smlen);
            memcpy(sm, g_valid_sm, smlen);
            sm[xrand() % smlen] ^= (uint8_t)(1u << (xrand() & 7));
            break;
        case 4:
        {
            size_t choices[5] = {
                (CRYPTO_BYTES > 0) ? CRYPTO_BYTES - 1 : 0,
                CRYPTO_BYTES,
                CRYPTO_BYTES + 1,
                g_valid_smlen - 1,
                rrange(0, g_valid_smlen)
            };
            smlen = choices[xrand() % 5];
            if (smlen > g_valid_smlen) smlen = g_valid_smlen;
            sm = (uint8_t *)malloc(smlen ? smlen : 1);
            memcpy(sm, g_valid_sm, smlen);
            break;
        }
        case 5:
        {
            size_t extra = rrange(1, 64);
            smlen = g_valid_smlen + extra;
            sm = (uint8_t *)malloc(smlen);
            memcpy(sm, g_valid_sm, g_valid_smlen);
            rfill(sm + g_valid_smlen, extra);
            break;
        }
        case 6:
            smlen = CRYPTO_BYTES;
            sm = (uint8_t *)malloc(smlen ? smlen : 1);
            rfill(sm, smlen);
            break;
        case 7:
            smlen = CRYPTO_BYTES + rrange(0, 48);
            sm = (uint8_t *)malloc(smlen ? smlen : 1);
            memcpy(sm, g_valid_sm, CRYPTO_BYTES);
            if (smlen > CRYPTO_BYTES) rfill(sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);
            break;
        case 8:
            smlen = g_valid_smlen;
            sm = (uint8_t *)malloc(smlen);
            memcpy(sm, g_valid_sm, smlen);
            rfill(sm, CRYPTO_BYTES);
            break;
        case 9:
            smlen = g_valid_smlen;
            sm = (uint8_t *)malloc(smlen);
            memcpy(sm, g_valid_sm, smlen);
            pk_own = (uint8_t *)malloc(CRYPTO_PUBLICKEYBYTES);
            rfill(pk_own, CRYPTO_PUBLICKEYBYTES);
            pk_use = pk_own;
            break;
        case 10:
            smlen = g_valid_smlen;
            sm = (uint8_t *)malloc(smlen);
            memcpy(sm, g_valid_sm, smlen);
            pk_own = (uint8_t *)malloc(CRYPTO_PUBLICKEYBYTES);
            memset(pk_own, 0xFF, CRYPTO_PUBLICKEYBYTES);
            pk_use = pk_own;
            break;
        default:
            smlen = rrange(0, max_smlen);
            sm = (uint8_t *)malloc(smlen ? smlen : 1);
            rfill(sm, smlen);
            pk_own = (uint8_t *)malloc(CRYPTO_PUBLICKEYBYTES);
            rfill(pk_own, CRYPTO_PUBLICKEYBYTES);
            pk_use = pk_own;
            break;
        }

        int rc = call_open(sm, smlen, pk_use);
        if (rc == 0) {
            accepted++;
            accept_by_mode[mode]++;
            note_accept(rc, pk_use, sm, smlen);
        } else rejected++;

        free(sm);
        if (pk_own) free(pk_own);
    }

    long skx_iters = iters / 4;
    void *skx = malloc(g_skx_bytes);
    for (long it = 0; it < skx_iters; ++it) {
        uint8_t *fsk = (uint8_t *)malloc(CRYPTO_SECRETKEYBYTES);
        int m = (int)(xrand() % 4);
        if (m == 0) rfill(fsk, CRYPTO_SECRETKEYBYTES);
        else if (m == 1) memset(fsk, 0x00, CRYPTO_SECRETKEYBYTES);
        else if (m == 2) memset(fsk, 0xFF, CRYPTO_SECRETKEYBYTES);
        else { memcpy(fsk, g_sk, CRYPTO_SECRETKEYBYTES); fsk[xrand() % CRYPTO_SECRETKEYBYTES] ^= 0x11; }
        (void)snova_sk_expand_skx(skx, fsk);
        free(fsk);
    }
    free(skx);

    long corpus_iters = iters / 2;
    size_t cbufcap = max_smlen + CRYPTO_PUBLICKEYBYTES + 8;
    uint8_t *cbuf = (uint8_t *)malloc(cbufcap);
    for (long it = 0; it < corpus_iters; ++it) {
        size_t clen = rrange(0, cbufcap);
        rfill(cbuf, clen);
        harness_consume(cbuf, clen);
    }
    free(cbuf);

    printf("fuzz_verify: iters=%ld  accepted=%ld  rejected=%ld  (sk_expand=%ld  corpus_replay=%ld)\n",
           iters, accepted, rejected, skx_iters, corpus_iters);
    printf("fuzz_verify: accept_by_mode = [");
    for (int i = 0; i < 12; ++i) printf("%ld%s", accept_by_mode[i], i < 11 ? "," : "");
    printf("]  nontrivial_accepts(FORGERY?)=%ld\n", g_accept_nontrivial);
    if (g_accept_nontrivial != 0) {
        fprintf(stderr, "FATAL: %ld non-trivial accepts — possible forgery/bounds bug!\n", g_accept_nontrivial);
        return 3;
    }
    printf("fuzz_verify: NO CRASH — all malformed inputs handled cleanly "
           "(ASan/UBSan clean if built with sanitizers)\n");

    free(g_pk); free(g_sk); free(g_valid_sm);
    return 0;
}

#endif
