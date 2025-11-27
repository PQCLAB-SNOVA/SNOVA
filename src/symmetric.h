#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
	uint64_t state[28];
} shake_t;

void shake128_init(shake_t *instance);
void shake256_init(shake_t *instance);
void shake_absorb(shake_t *instance, const uint8_t *in, size_t inlen);
void shake_finalize(shake_t *instance);
void shake_squeeze(uint8_t *out, size_t outlen, shake_t *instance);
void shake_squeeze_keep(uint8_t *out, size_t outlen, shake_t *instance);
void shake_release(shake_t *instance);
void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

typedef struct {
#if __AVX512F__
	uint64_t state[406];
#elif __AVX2__
	uint64_t state[206];
#else
	uint64_t state[56];
#endif
} snova_pk_expander_t;

void snova_pk_expander_init(snova_pk_expander_t *instance, const uint8_t *seed, size_t input_bytes);
void snova_pk_expander(uint8_t *data, size_t num_bytes, snova_pk_expander_t *instance);

#endif
