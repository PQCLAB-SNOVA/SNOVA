/**
 * @file drbg.c
 */
#include "drbg.h"
#include "../primitives/sym_aes.h"
#include "secure_clear.h"
#include <string.h>

static struct {
    uint8_t Key[32];
    uint8_t V[16];
    int reseed_counter;
} S;

/**
 * @param V 16-byte counter.
 * @return none.
 */
static void inc_V(uint8_t V[16]) {
    for (int i = 15; i >= 0; --i) { if (++V[i]) break; }
}

/**
 * @param provided 48-byte provided data or NULL.
 * @return none.
 */
static void drbg_update(const uint8_t *provided) {
    uint8_t tmp[48];
    uint8_t rk[240];
    aes256_keyexp(S.Key, rk);
    for (int i = 0; i < 3; ++i) {
        inc_V(S.V);
        aes256_encrypt(rk, S.V, tmp + 16 * i);
    }
    if (provided)
        for (int i = 0; i < 48; ++i) tmp[i] ^= provided[i];
    memcpy(S.Key, tmp, 32);
    memcpy(S.V, tmp + 32, 16);
    SNOVA_CLEAR_OBJ(tmp);
    SNOVA_CLEAR_OBJ(rk);
}

void randombytes_init(const uint8_t *entropy_input,
                       const uint8_t *personalization, int security_strength) {
    (void)security_strength;
    uint8_t seed[48];
    memcpy(seed, entropy_input, 48);
    if (personalization)
        for (int i = 0; i < 48; ++i) seed[i] ^= personalization[i];
    memset(S.Key, 0, 32);
    memset(S.V, 0, 16);
    drbg_update(seed);
    S.reseed_counter = 1;
    SNOVA_CLEAR_OBJ(seed);
}

int randombytes(uint8_t *x, size_t xlen) {
    uint8_t block[16];
    uint8_t rk[240];
    size_t i = 0;
    aes256_keyexp(S.Key, rk);
    while (xlen > 0) {
        inc_V(S.V);
        aes256_encrypt(rk, S.V, block);
        size_t n = (xlen > 16) ? 16 : xlen;
        memcpy(x + i, block, n);
        i += n;
        xlen -= n;
    }
    drbg_update(NULL);
    S.reseed_counter++;
    SNOVA_CLEAR_OBJ(block);
    SNOVA_CLEAR_OBJ(rk);
    return 0;
}
