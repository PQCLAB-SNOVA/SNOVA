// SPDX-License-Identifier: MIT
/**
 * Vectorized SHAKE128 CTR-XOF based on XKCP.
 *
 * Copyright (c) 2024 SNOVA TEAM
 */

#include <stdlib.h>

#include "snova_shake.h"

#ifndef PARALLELISM
#if __AVX512F__
#define PARALLELISM 8
#elif __AVX2__
#define PARALLELISM 4
#else
#define PARALLELISM 1
#endif
#endif

#if PARALLELISM == 8
#include "KeccakP-1600-times8-SnP.h"
#elif PARALLELISM == 4
#include "KeccakP-1600-times4-SnP.h"
#endif

#if PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN
#error "PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN"
#endif

/**
 * Squeeze bytes in parallel.
 */
int snova_shake_squeeze(Keccak_HashInstance *keccak_instance, uint64_t *data, size_t num_bytes)
{
    uint8_t prepared_state[200 * PARALLELISM];
    ALIGN(PARALLELISM * 8)
    uint8_t states[200 * PARALLELISM];
    uint64_t *states64 = (uint64_t *)states;
    uint64_t *data64 = data;

    KeccakWidth1600_SpongeInstance *sponge = &keccak_instance->sponge;
    uint32_t bytes_rate = sponge->rate / 8;

    uint64_t block = 0;
    uint64_t index = 0;
    size_t last_idx = num_bytes;

    for (int idx = 0; idx < PARALLELISM; idx++)
    {
        uint64_t *state64 = (uint64_t *)sponge->state;
        uint64_t *prep64 = (uint64_t *)prepared_state;

        for (int idx2 = 0; idx2 < 25; idx2++)
            prep64[idx + idx2 * PARALLELISM] = state64[idx2];
        // SHAKE padding. Use the (uint8_t *)prepared_state here
        prepared_state[idx * 8 + (sponge->byteIOIndex + 8) * PARALLELISM] ^= keccak_instance->delimitedSuffix;
        prepared_state[idx * 8 + (bytes_rate - 8) * PARALLELISM + 7] ^= 0x80;
    }

    while (index < last_idx)
    {
        uint32_t byteIOIndex = sponge->byteIOIndex;
        memcpy(states, prepared_state, PARALLELISM * 200);
        for (int idx = 0; idx < PARALLELISM; idx++)
        {
            states64[byteIOIndex * PARALLELISM / 8 + idx] ^= block;
            block++;
        }

#if PARALLELISM == 1
        KeccakP1600_Permute_24rounds(states);
#elif PARALLELISM == 4
        KeccakP1600times4_PermuteAll_24rounds(states);
#elif PARALLELISM == 8
        KeccakP1600times8_PermuteAll_24rounds(states);
#else
#error "PARALLELISM must be 1, 4 or 8"
#endif

        for (int idx = 0; idx < PARALLELISM; idx++)
        {
            size_t bytes = last_idx - index;
            if (bytes > bytes_rate)
                bytes = bytes_rate;

            uint64_t *source64 = states64 + idx;
            for (size_t idx2 = 0; idx2 < bytes / 8; idx2++)
            {
                *data64 = *source64;
                source64 += PARALLELISM;
                data64++;
            }
            index += bytes;

            if (index >= last_idx)
                break;
        }
    }

    return 0;
}

/**
 * Generate XOF data from a short seed
 */
void snova_shake_opt(const uint8_t *seed, size_t input_bytes, uint64_t *output, size_t output_bytes)
{
    Keccak_HashInstance keccak_instance;
    Keccak_HashInitialize_SHAKE128(&keccak_instance);
    Keccak_HashUpdate(&keccak_instance, seed, 8 * input_bytes);
    snova_shake_squeeze(&keccak_instance, output, output_bytes);
}
