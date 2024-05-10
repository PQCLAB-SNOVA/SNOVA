/**
 * Vectorized SHAKE128 CTR-XOF based on XKCP.
 */

#include <stdlib.h>
#include <assert.h>

#include "vexof.h"

// Sanity check
#if KeccakP1600_stateSizeInBytes != 200
#error "KeccakP1600_stateSizeInBytes must be 200"
#endif

#if PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN
#error "PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN"
#endif

#if PARALLELISM == 8
#include "../shake/KeccakP-1600-times8-SnP.h"
#elif PARALLELISM == 4
#include "../shake/KeccakP-1600-times4-SnP.h"
#endif

/**
 * Create VeXOF instance
 */
int VeXOF_HashInitialize(VeXOF_Instance *vexof_instance)
{
    vexof_instance->squeezing = 0;
    return Keccak_HashInitialize_SHAKE128(&vexof_instance->keccak_instance);
}

/**
 * Add bytes to instance
 */
int VeXOF_HashUpdate(VeXOF_Instance *vexof_instance, const uint8_t *data, size_t bytes)
{
    vexof_instance->squeezing = 0;
    return Keccak_HashUpdate(&vexof_instance->keccak_instance, data, 8 * bytes);
}

/**
 * Squeeze bytes in parallel.
 */
int VeXOF_Squeeze(VeXOF_Instance *vexof_instance, uint64_t *data, size_t num_bytes)
{
    assert(num_bytes % 8 == 0);

    KeccakWidth1600_SpongeInstance *sponge = &vexof_instance->keccak_instance.sponge;
    uint32_t bytes_rate = sponge->rate / 8;

    if (!vexof_instance->squeezing)
    {
        assert(sponge->byteIOIndex % 8 == 0);
        assert(sponge->byteIOIndex < (bytes_rate - 10));

        uint64_t *state64 = (uint64_t *)sponge->state;
        uint64_t *prep64 = (uint64_t *)vexof_instance->prepared_state;
        for (int idx = 0; idx < PARALLELISM; idx++)
        {
            for (int idx2 = 0; idx2 < 25; idx2++)
                prep64[idx + idx2 * PARALLELISM] = state64[idx2];
            // SHAKE padding. Use the (uint8_t *)prepared_state here
            vexof_instance->prepared_state[idx * 8 + (sponge->byteIOIndex + 8) * PARALLELISM] ^= vexof_instance->keccak_instance.delimitedSuffix;
            vexof_instance->prepared_state[idx * 8 + (bytes_rate - 8) * PARALLELISM + 7] ^= 0x80;
        }

        vexof_instance->block = 0;
        vexof_instance->index = 0;
        vexof_instance->squeezing = 1;
    }

    uint8_t *states = &vexof_instance->states_data[0];
    size_t last_idx = vexof_instance->index + num_bytes;
    uint64_t *states64 = (uint64_t *)states;
    uint64_t *data64 = data;

    // Squeeze bytes already created in a preceding invocation
    uint32_t mod_index = vexof_instance->index % (PARALLELISM * bytes_rate);
    if (mod_index)
    {
        size_t remaining = PARALLELISM * bytes_rate - mod_index;
        if (remaining > num_bytes)
            remaining = num_bytes;

        for (size_t idx = 0; idx < remaining / 8; idx++)
        {
            uint32_t idx1 = mod_index + 8 * idx;
            uint32_t idx2 = (idx1 % bytes_rate) / 8 * PARALLELISM + (idx1 / bytes_rate);
            data64[idx] = states64[idx2];
        }

        data64 += remaining / 8;
        vexof_instance->index += remaining;
    }

    while (vexof_instance->index < last_idx)
    {
        uint32_t byteIOIndex = sponge->byteIOIndex;
        memcpy(states, vexof_instance->prepared_state, PARALLELISM * 200);
        for (int idx = 0; idx < PARALLELISM; idx++)
        {
            *((uint64_t *)&states[byteIOIndex * PARALLELISM + 8 * idx]) ^= vexof_instance->block;
            vexof_instance->block++;
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
            size_t bytes = last_idx - vexof_instance->index;
            if (bytes > bytes_rate)
                bytes = bytes_rate;

            uint64_t *source64 = states64 + idx;
            for (size_t idx2 = 0; idx2 < bytes / 8; idx2++)
            {
                *data64 = *source64;
                source64 += PARALLELISM;
                data64++;
            }
            vexof_instance->index += bytes;

            if (vexof_instance->index >= last_idx)
                break;
        }
    }

    return 0;
}

/**
 * Generate XOF data from a short seed
 */
void vexof(const uint8_t *seed, size_t input_bytes, uint64_t *output, size_t output_bytes)
{
    assert(input_bytes <= 152);

    if (output_bytes > 168)
    {
        VeXOF_Instance vexofInstance;
        VeXOF_HashInitialize(&vexofInstance);
        VeXOF_HashUpdate(&vexofInstance, seed, input_bytes);
        VeXOF_Squeeze(&vexofInstance, output, output_bytes);
    }
    else
    {
        Keccak_HashInstance keccak_instance;
        Keccak_HashInitialize_SHAKE128(&keccak_instance);
        Keccak_HashUpdate(&keccak_instance, seed, 8 * input_bytes);
        // Add the block index (0)
        keccak_instance.sponge.byteIOIndex += 8;
        Keccak_HashFinal(&keccak_instance, NULL);
        Keccak_HashSqueeze(&keccak_instance, (uint8_t *)output, 8 * output_bytes);
    }
}
