/**
 * @file primitives/snova_aes.h
 */
#ifndef SNOVA_PRIMITIVES_SNOVA_AES_H
#define SNOVA_PRIMITIVES_SNOVA_AES_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void AES_256_ECB(const unsigned char *key, const uint8_t *input, unsigned char *output);

int AES_128_CTR(unsigned char *output, size_t outputByteLen,
                const unsigned char *input, size_t inputByteLen);

#ifdef __cplusplus
}
#endif

#endif
