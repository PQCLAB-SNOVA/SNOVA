// SPDX-License-Identifier: MIT

/**
 * Glue code between SNOVA and liboqs for the symmetric primitives used by SNOVA
 *
 * Copyright (c) 2025 SNOVA TEAM
 */

#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include <oqs/sha3.h>
#include <stdalign.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

typedef struct {
	uint64_t state[28];
} shake_t;

typedef struct {
#if __AVX512F__
	uint64_t state[406];
#elif __AVX2__
	uint64_t state[206];
#else
	uint64_t state[56];
#endif
} snova_pk_expander_t;

static inline void shake256_init(shake_t* instance) {
	OQS_SHA3_shake256_inc_ctx* oqs_instance = (OQS_SHA3_shake256_inc_ctx*)instance;
	OQS_SHA3_shake256_inc_init(oqs_instance);
}

static inline void shake_absorb(shake_t* instance, const uint8_t* in, size_t inlen) {
	OQS_SHA3_shake256_inc_ctx* oqs_instance = (OQS_SHA3_shake256_inc_ctx*)instance;
	OQS_SHA3_shake256_inc_absorb(oqs_instance, in, inlen);
}

static inline void shake_finalize(shake_t* instance) {
	OQS_SHA3_shake256_inc_ctx* oqs_instance = (OQS_SHA3_shake256_inc_ctx*)instance;
	OQS_SHA3_shake256_inc_finalize(oqs_instance);
}

static inline void shake_squeeze(uint8_t* out, size_t outlen, shake_t* instance) {
	OQS_SHA3_shake256_inc_ctx* oqs_instance = (OQS_SHA3_shake256_inc_ctx*)instance;
	OQS_SHA3_shake256_inc_squeeze(out, outlen, oqs_instance);
	OQS_SHA3_shake256_inc_ctx_release(oqs_instance);
}

static inline void shake_squeeze_keep(uint8_t* out, size_t outlen, shake_t* instance) {
	OQS_SHA3_shake256_inc_ctx* oqs_instance = (OQS_SHA3_shake256_inc_ctx*)instance;
	OQS_SHA3_shake256_inc_squeeze(out, outlen, oqs_instance);
}

static inline void shake_release(shake_t* instance) {
	OQS_SHA3_shake256_inc_ctx* oqs_instance = (OQS_SHA3_shake256_inc_ctx*)instance;
	OQS_SHA3_shake256_inc_ctx_release(oqs_instance);
}

static inline void shake256(uint8_t* out, size_t outlen, const uint8_t* in, size_t inlen) {
	shake_t instance;
	shake256_init(&instance);
	shake_absorb(&instance, in, inlen);
	shake_finalize(&instance);
	shake_squeeze(out, outlen, &instance);
}

/**
 * SNOVA-SHAKE XOF
 */

#ifdef AESCTR

#include <aes.h>

static inline void snova_pk_expander_init(snova_pk_expander_t* instance, const uint8_t* seed, size_t input_bytes) {
	(void)input_bytes;
	memcpy(instance, seed, 16);
}

static inline void snova_pk_expander(uint8_t* data, size_t num_bytes, snova_pk_expander_t* instance) {
	const unsigned char *input = (unsigned char*)instance;
	const uint8_t iv[16] = {0};

	void *state;
	OQS_AES128_CTR_inc_init(input, &state);
	OQS_AES128_CTR_inc_stream_iv(iv, 12, state, data, num_bytes);
	OQS_AES128_free_schedule(state);
}

#elif defined(OQS_ENABLE_SHA3_xkcp_low_avx2)

/**
 * liboqs x4 version
 */
#include <oqs/sha3x4.h>

static inline void snova_pk_expander_init(snova_pk_expander_t* instance, const uint8_t* seed, size_t input_bytes) {
	(void)input_bytes;
	memcpy(instance, seed, 16);
}

static inline void snova_pk_expander(uint8_t* data, size_t num_bytes, snova_pk_expander_t* instance) {
	uint8_t *pt_seed_array = (uint8_t*)instance;
	size_t input_bytes = SEED_LENGTH_PUBLIC;
	size_t output_bytes = num_bytes;
	size_t index = 0;
	uint64_t block = 0;

	while (index < output_bytes) {
		OQS_SHA3_shake128_x4_inc_ctx hashInstance;
		OQS_SHA3_shake128_x4_inc_init(&hashInstance);
		OQS_SHA3_shake128_x4_inc_absorb(&hashInstance, pt_seed_array, pt_seed_array, pt_seed_array, pt_seed_array, input_bytes);

		// Turn SHAKE128 into SHAKE128 CTR-XOF
		// Little endian
		uint64_t block_0 = block;
		block++;
		uint64_t block_1 = block;
		block++;
		uint64_t block_2 = block;
		block++;
		uint64_t block_3 = block;
		block++;
		OQS_SHA3_shake128_x4_inc_absorb(&hashInstance, (uint8_t*)&block_0, (uint8_t*)&block_1, (uint8_t*)&block_2,
		                                (uint8_t*)&block_3, 8);

		OQS_SHA3_shake128_x4_inc_finalize(&hashInstance);
		size_t bytes = output_bytes - index;
		if (bytes > 4 * 168) {
			OQS_SHA3_shake128_x4_inc_squeeze(data, data + 168, data + 336, data + 504, 168, &hashInstance);
			data += 4 * 168;
		} else {
			// Last round
			alignas(32) uint8_t buf[4 * 168];
			OQS_SHA3_shake128_x4_inc_squeeze(buf, buf + 168, buf + 336, buf + 504, 168, &hashInstance);
			memcpy(data, buf, bytes);
		}
		index += 4 * 168;
		OQS_SHA3_shake128_x4_inc_ctx_release(&hashInstance);
	}
}

#else

/**
 * liboqs version of reference
 */

static inline void snova_pk_expander_init(snova_pk_expander_t* instance, const uint8_t* seed, size_t input_bytes) {
	(void)input_bytes;
	memcpy(instance, seed, 16);
}

static inline void snova_pk_expander(uint8_t* data, size_t num_bytes, snova_pk_expander_t* instance) {
	uint8_t *pt_seed_array = (uint8_t*)instance;
	size_t input_bytes = SEED_LENGTH_PUBLIC;
	size_t output_bytes = num_bytes;
	size_t index = 0;
	uint64_t block = 0;

	while (index < output_bytes) {
		OQS_SHA3_shake128_inc_ctx hashInstance;
		OQS_SHA3_shake128_inc_init(&hashInstance);
		OQS_SHA3_shake128_inc_absorb(&hashInstance, pt_seed_array, input_bytes);

		// Turn SHAKE128 into SHAKE128 CTR-XOF
		for (int idx = 0; idx < 8; idx++) {
			// Little endian
			uint8_t block_i = (block >> (8 * idx)) & 0xff;
			OQS_SHA3_shake128_inc_absorb(&hashInstance, &block_i, 1);
		}

		OQS_SHA3_shake128_inc_finalize(&hashInstance);
		size_t bytes = output_bytes - index;
		if (bytes > 168) {
			bytes = 168;
		}

		OQS_SHA3_shake128_inc_squeeze(data, bytes, &hashInstance);
		OQS_SHA3_shake128_inc_ctx_release(&hashInstance);

		block++;
		data += bytes;
		index += bytes;
	}
}

#endif

static void snova_pk_expand(uint8_t* out, size_t outlen, const uint8_t* in, size_t inlen) {
	snova_pk_expander_t instance;
	snova_pk_expander_init(&instance, in, inlen);
	snova_pk_expander(out, outlen, &instance);
}

#endif
