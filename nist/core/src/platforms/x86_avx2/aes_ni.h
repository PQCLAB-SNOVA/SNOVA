/**
 * @file aes.h
 */
#ifndef TP_AES_H
#define TP_AES_H

#include <stddef.h>
#include <stdint.h>

void aes128_keyexp(const uint8_t key[16], uint8_t rk[176]);

void aes128_encrypt(const uint8_t rk[176], const uint8_t in[16], uint8_t out[16]);

void aes256_keyexp(const uint8_t key[32], uint8_t rk[240]);

void aes256_encrypt(const uint8_t rk[240], const uint8_t in[16], uint8_t out[16]);

void aes128_ctr_zero(uint8_t *out, size_t outlen, const uint8_t key[16]);

#if (defined(SNOVA_VERIFY_STREAM) && SNOVA_VERIFY_STREAM) || \
    (defined(SNOVA_SIGN_STREAM) && SNOVA_SIGN_STREAM) || \
    (defined(SNOVA_KEYGEN_STREAM) && SNOVA_KEYGEN_STREAM) || \
    (defined(SNOVA_PKX_PGEN) && SNOVA_PKX_PGEN)
void aes128_ctr_zero_at(uint8_t *out, size_t outlen, const uint8_t key[16],
                        uint64_t block_offset);
#endif

#endif
