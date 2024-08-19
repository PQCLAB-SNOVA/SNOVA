// SPDX-License-Identifier: MIT
/**
 * Vectorized SHAKE128 CTR-XOF based on XKCP.
 *
 * Copyright (c) 2024 SNOVA TEAM
 */

#include <stdlib.h>
#include <assert.h>

#include "snova_shake.h"

/**
 * Reference version
 */
void snova_shake_ref(const uint8_t *pt_seed_array, size_t input_bytes, uint64_t *data64,
                     size_t output_bytes)
{
    assert(input_bytes <= 152);
    assert(output_bytes % 8 == 0);

    size_t index = 0;
    uint64_t block = 0;
    uint8_t *data = (uint8_t *)data64;

    while (index < output_bytes)
    {
        Keccak_HashInstance hashInstance;
        Keccak_HashInitialize_SHAKE128(&hashInstance);
        Keccak_HashUpdate(&hashInstance, pt_seed_array, 8 * input_bytes);

        // Turn SHAKE128 into SHAKE128 CTR-XOF
        for (int idx = 0; idx < 8; idx++)
        {
            // Little endian
            uint8_t block_i = (block >> (8 * idx)) & 0xff;
            Keccak_HashUpdate(&hashInstance, &block_i, 8);
        }

        Keccak_HashFinal(&hashInstance, NULL);
        size_t bytes = output_bytes - index;
        if (bytes > (hashInstance.sponge.rate / 8))
            bytes = (hashInstance.sponge.rate / 8);

        Keccak_HashSqueeze(&hashInstance, data, 8 * bytes);

        block++;
        data += bytes;
        index += bytes;
    }
}
