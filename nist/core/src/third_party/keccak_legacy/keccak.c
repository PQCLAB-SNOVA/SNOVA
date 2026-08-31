/*
 * @file keccak.c
 *        FIPS-202 Keccak-f[1600] and the SHAKE sponge (standard primitive).
 * @details
 *   English: compact Keccak written from FIPS-202; classic round structure
 *   (theta, rho+pi, chi, iota). Not SNOVA; lives in third_party/. Correctness
 *   gated by FIPS-202 known-answer vectors in test/sym_selftest.c.
 */
#include "keccak.h"
#include <string.h>

#define ROL64(a, o) (((a) << (o)) | ((a) >> (64 - (o))))

static const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
    0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
    0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
    0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
    0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL};

static const int RHO[24] = {1,  3,  6,  10, 15, 21, 28, 36, 45, 55, 2,  14,
                            27, 41, 56, 8,  25, 43, 62, 18, 39, 61, 20, 44};
static const int PI[24] = {10, 7,  11, 17, 18, 3, 5,  16, 8, 21, 24, 4,
                           15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6,  1};

static void keccakf(uint64_t s[25]) {
    for (int r = 0; r < 24; ++r) {
        uint64_t bc[5], t;
        /* θ / theta */
        for (int i = 0; i < 5; ++i)
            bc[i] = s[i] ^ s[i + 5] ^ s[i + 10] ^ s[i + 15] ^ s[i + 20];
        for (int i = 0; i < 5; ++i) {
            t = bc[(i + 4) % 5] ^ ROL64(bc[(i + 1) % 5], 1);
            for (int j = 0; j < 25; j += 5) s[j + i] ^= t;
        }
        /* ρ + π / rho + pi */
        t = s[1];
        for (int i = 0; i < 24; ++i) {
            int j = PI[i];
            uint64_t tmp = s[j];
            s[j] = ROL64(t, RHO[i]);
            t = tmp;
        }
        /* χ / chi */
        for (int j = 0; j < 25; j += 5) {
            for (int i = 0; i < 5; ++i) bc[i] = s[j + i];
            for (int i = 0; i < 5; ++i)
                s[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
        }
        /* ι / iota */
        s[0] ^= RC[r];
    }
}

/** @copydoc shake_init */
void shake_init(keccak_ctx *c, int bits) {
    memset(c, 0, sizeof(*c));
    c->rate = (bits == 128) ? 168 : 136; /* 200 - 2*bits/8 / capacity rule */
    c->delim = 0x1F;
    c->pos = 0;
    c->squeezing = 0;
}

/** @copydoc shake_absorb */
void shake_absorb(keccak_ctx *c, const uint8_t *in, size_t inlen) {
    uint8_t *sb = (uint8_t *)c->st;
    for (size_t k = 0; k < inlen; ++k) {
        sb[c->pos++] ^= in[k];
        if (c->pos == c->rate) {
            keccakf(c->st);
            c->pos = 0;
        }
    }
}

/** @copydoc shake_finalize */
void shake_finalize(keccak_ctx *c) {
    uint8_t *sb = (uint8_t *)c->st;
    sb[c->pos] ^= c->delim;
    sb[c->rate - 1] ^= 0x80;
    keccakf(c->st);
    c->pos = 0;
    c->squeezing = 1;
}

/** @copydoc shake_squeeze */
void shake_squeeze(keccak_ctx *c, uint8_t *out, size_t outlen) {
    const uint8_t *sb = (const uint8_t *)c->st;
    for (size_t k = 0; k < outlen; ++k) {
        if (c->pos == c->rate) {
            keccakf(c->st);
            c->pos = 0;
        }
        out[k] = sb[c->pos++];
    }
}

static void shake_oneshot(int bits, uint8_t *out, size_t outlen,
                          const uint8_t *in, size_t inlen) {
    keccak_ctx c;
    shake_init(&c, bits);
    shake_absorb(&c, in, inlen);
    shake_finalize(&c);
    shake_squeeze(&c, out, outlen);
}

/** @copydoc shake128 */
void shake128(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen) {
    shake_oneshot(128, out, outlen, in, inlen);
}

/** @copydoc shake256 */
void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen) {
    shake_oneshot(256, out, outlen, in, inlen);
}
