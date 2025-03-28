// SPDX-License-Identifier: MIT
/**
 * Vectorized SNOVA-SHAKE XOF
 *
 * The optimised implementation can be tested against the reference using the genkat function of SNOVA
 * for optimisation levels OPTIMISATION=0 versus OPTIMISATION=1 or OPTIMISATION=2.
 */

#include <assert.h>
#include <string.h>

#include "snova.h"

/**
 * Reference version
 */
void snova_shake(const uint8_t *pt_seed_array, size_t input_bytes, uint64_t *data64, size_t output_bytes) {
	assert(input_bytes <= 152);
	assert(output_bytes % 8 == 0);

	size_t index = 0;
	uint64_t block = 0;
	uint8_t *data = (uint8_t *)data64;

	while (index < output_bytes) {
		Keccak_HashInstance hashInstance;
		Keccak_HashInitialize_SHAKE128(&hashInstance);
		Keccak_HashUpdate(&hashInstance, pt_seed_array, 8 * input_bytes);

		// Turn SHAKE128 into SHAKE128 CTR-XOF
		for (int idx = 0; idx < 8; idx++) {
			// Little endian
			uint8_t block_i = (block >> (8 * idx)) & 0xff;
			Keccak_HashUpdate(&hashInstance, &block_i, 8);
		}

		Keccak_HashFinal(&hashInstance, NULL);
		size_t bytes = output_bytes - index;
		if (bytes > (hashInstance.sponge.rate / 8)) {
			bytes = (hashInstance.sponge.rate / 8);
		}

		Keccak_HashSqueeze(&hashInstance, data, 8 * bytes);

		block++;
		data += bytes;
		index += bytes;
	}
}
