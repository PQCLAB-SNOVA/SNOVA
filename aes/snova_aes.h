/**
 * AES primitives used by SNOVA
 * 
 * Copyright (c) 2024 SNOVA TEAM
 */

#ifndef SNOVA_AES_H
#define SNOVA_AES_H

#include <stdint.h>

int AES_128_CTR(unsigned char *output, size_t outputByteLen, const unsigned char *input, size_t inputByteLen);
void AES_256_ECB(const unsigned char *key, const uint8_t *input, unsigned char *output);

#endif
