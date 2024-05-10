/**
 * 2024 Jan Adriaan Leegwater
 *
 * Test optimized vexof implementation against reference
 */

#include <stdio.h>
#include "vexof.h"

#define NUM_XOF_BYTES 100000

/**
 * Reference version
 */
void vexof_ref(const uint8_t *pt_seed_array, size_t input_bytes, uint8_t *data,
               size_t output_bytes)
{
    size_t index = 0;
    uint64_t block = 0;

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

int main()
{
    uint8_t pt_public_key_seed[200];
    memset(pt_public_key_seed, 1, 200);

    uint8_t prng_output_reference[NUM_XOF_BYTES] = {0};
    uint64_t prng_output_vexof[NUM_XOF_BYTES / 8] = {0};
    uint8_t *prng_output_vexof8 = (uint8_t *)prng_output_vexof;

    vexof(pt_public_key_seed, 16, prng_output_vexof, NUM_XOF_BYTES);
    vexof_ref(pt_public_key_seed, 16, prng_output_reference, NUM_XOF_BYTES);

    int testok = 1;
    for (int idx = 0; idx < NUM_XOF_BYTES; idx++)
        if (prng_output_vexof8[idx] != prng_output_reference[idx])
        {
            printf("Test Failed @ %d: %02x %02x\n", idx, prng_output_vexof8[idx],
                   prng_output_reference[idx]);
            testok = 0;
            break;
        }

    if (testok)
    {
        printf("Test ok\n");
    }

    return 0;
}
