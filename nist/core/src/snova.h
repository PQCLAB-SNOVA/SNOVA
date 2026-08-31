/**
 * @file snova.h
 */
#ifndef M0_SNOVA_H
#define M0_SNOVA_H

#include <stdint.h>
#include <stddef.h>
#include "snova_params.h"
#include "gf16_core/gf16m.h"

#if SNOVA_Q != 16
#define SNOVA_PK_BYTES (SNOVA_SEED_PUB + SNOVA_BYTES_GF(SNOVA_M1*SNOVA_O*SNOVA_O*SNOVA_L2))
#else
#define SNOVA_PK_BYTES (SNOVA_SEED_PUB + ((SNOVA_M1*SNOVA_O*SNOVA_O*SNOVA_L2 + 1) >> 1))
#endif
#define SNOVA_SK_BYTES (SNOVA_SEED_LEN + SNOVA_BYTES_PK_HASH)
#define SNOVA_SM_OVERHEAD (SNOVA_BYTES_SIG + SNOVA_SALT_BYTES)

int snova_crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

int snova_crypto_sign(uint8_t *sm, size_t *smlen, const uint8_t *m,
                      size_t mlen, const uint8_t *sk);

int snova_crypto_sign_open(uint8_t *m, size_t *mlen, const uint8_t *sm,
                           size_t smlen, const uint8_t *pk);

int snova_keygen(uint8_t *pk, uint8_t *sk, const uint8_t *seed);

size_t snova_skx_bytes(void);

void snova_sk_expand_skx(uint8_t *skx, const uint8_t *sk);

int snova_sign_digest_skx(uint8_t *sm, const uint8_t *digest, size_t dlen,
                          const uint8_t *salt, const uint8_t *skx);

void snova_skx_clear(uint8_t *skx);

size_t snova_pkx_bytes(void);

int snova_pk_expand(uint8_t *pkx, const uint8_t *pk);

int snova_verify_digest_pkx(const uint8_t *digest, size_t dlen,
                            const uint8_t *sm, const uint8_t *pkx);

#endif
