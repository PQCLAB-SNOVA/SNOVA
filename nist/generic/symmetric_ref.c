// SPDX-License-Identifier: MIT

/**
 * Reference implementation of the symmetric primitives used by SNOVA
 *
 * Copyright (c) 2025 SNOVA TEAM
 */

#include <assert.h>

#include "keccak_opt64.h"
#include "snova.h"
#include "symmetric.h"

void shake128_init(shake_t* instance) {
	memset(instance, 0, sizeof(shake_t));
	instance->state[26] = 168;
}

void shake256_init(shake_t* instance) {
	memset(instance, 0, sizeof(shake_t));
	instance->state[26] = 136;
}

void shake_absorb(shake_t* instance, const uint8_t* in, size_t inlen) {
	keccak_inc_absorb(instance->state, instance->state[26], in, inlen);
}

void shake_finalize(shake_t* instance) {
	keccak_inc_finalize(instance->state, instance->state[26], 0x1F);
}

void shake_squeeze(uint8_t* out, size_t outlen, shake_t* instance) {
	keccak_inc_squeeze(out, outlen, instance->state, instance->state[26]);
}

void shake_squeeze_keep(uint8_t* out, size_t outlen, shake_t* instance) {
	keccak_inc_squeeze(out, outlen, instance->state, instance->state[26]);
}

void shake_release(shake_t* instance) {
	(void)instance;
}

void shake256(uint8_t* out, size_t outlen, const uint8_t* in, size_t inlen) {
	shake_t instance;
	shake256_init(&instance);
	shake_absorb(&instance, in, inlen);
	shake_finalize(&instance);
	shake_squeeze(out, outlen, &instance);
}

/**
 * Reference SNOVA public key expander
 */
typedef struct {
	uint64_t state[56];
} snova_pk_expander_t;

#if defined(AESCTR)

void snova_pk_expander_init(snova_pk_expander_t* instance, const uint8_t* seed, size_t input_bytes) {
	(void)input_bytes;
	memcpy(instance, seed, 16);
}

#define AES128_CTR SNOVA_NAMESPACE(AES128_CTR)
int AES128_CTR(unsigned char* output, size_t outputByteLen, const unsigned char* input, size_t inputByteLen);

void snova_pk_expander(uint8_t* data, size_t num_bytes, snova_pk_expander_t* instance) {
	AES128_CTR(data, num_bytes, (uint8_t*)instance, 16);
}

#else

void snova_pk_expander_init(snova_pk_expander_t* ref_instance, const uint8_t* seed, size_t input_bytes) {
	assert(input_bytes <= 152);

	shake_t *instance = (shake_t*)ref_instance;
	shake128_init(instance);
	shake_absorb(instance, seed, input_bytes);
}

void snova_pk_expander(uint8_t* data, size_t num_bytes, snova_pk_expander_t* instance) {
	uint64_t index = 0;
	uint64_t block = 0;
	size_t rate = ((shake_t*)instance)->state[26];

	while (index < num_bytes) {
		shake_t block_instance;
		memcpy(&block_instance, instance, sizeof(shake_t));

		// Turn SHAKE128 into SHAKE128 CTR-XOF
		for (int idx = 0; idx < 8; idx++) {
			// Little endian
			uint8_t block_i = (block >> (8 * idx)) & 0xff;
			shake_absorb(&block_instance, &block_i, 1);
		}

		shake_finalize(&block_instance);
		size_t bytes = num_bytes - index;
		if (bytes > rate) {
			bytes = rate;
		}

		shake_squeeze(data, bytes, &block_instance);

		block++;
		data += bytes;
		index += bytes;
	}
}

#endif

void snova_pk_expand(uint8_t* out, size_t outlen, const uint8_t* in, size_t inlen) {
	snova_pk_expander_t instance;
	snova_pk_expander_init(&instance, in, inlen);
	snova_pk_expander(out, outlen, &instance);
}
