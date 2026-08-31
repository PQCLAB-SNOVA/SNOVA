/*
 * @file keccak.h
 *        FIPS-202 Keccak / SHAKE128 / SHAKE256 (standard-primitive API).
 * @details
 *   English: a compact implementation written from the FIPS-202 spec; NOT
 *   part of the SNOVA algorithm — a standard primitive kept under
 *   third_party/. Provides one-shot and incremental (multi-absorb then
 *   squeeze) APIs; SNOVA signing needs incremental SHAKE256.
 */
#ifndef TP_KECCAK_H
#define TP_KECCAK_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint64_t st[25];
    size_t   pos;
    size_t   rate;
    uint8_t  delim;
    int      squeezing;
} keccak_ctx;

void shake_init(keccak_ctx *c, int bits);

void shake_absorb(keccak_ctx *c, const uint8_t *in, size_t inlen);

void shake_finalize(keccak_ctx *c);

void shake_squeeze(keccak_ctx *c, uint8_t *out, size_t outlen);

void shake128(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

#endif /* TP_KECCAK_H */
