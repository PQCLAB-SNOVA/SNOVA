/**
 * @file platforms/ref/aes_ref.h
 */
#ifndef SNOVA_PLATFORMS_REF_AES_REF_H
#define SNOVA_PLATFORMS_REF_AES_REF_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void aes128_keyexp_ref(const uint8_t key[16], uint8_t rk[176]);
void aes128_encrypt_ref(const uint8_t rk[176], const uint8_t in[16], uint8_t out[16]);
void aes256_keyexp_ref(const uint8_t key[32], uint8_t rk[240]);
void aes256_encrypt_ref(const uint8_t rk[240], const uint8_t in[16], uint8_t out[16]);
void aes128_ctr_zero_ref(uint8_t *out, size_t outlen, const uint8_t key[16]);

#define aes128_keyexp     aes128_keyexp_ref
#define aes128_encrypt    aes128_encrypt_ref
#define aes256_keyexp     aes256_keyexp_ref
#define aes256_encrypt    aes256_encrypt_ref
#define aes128_ctr_zero   aes128_ctr_zero_ref

#ifdef __cplusplus
}
#endif

#endif
