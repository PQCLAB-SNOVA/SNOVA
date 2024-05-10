/**
 * Vectorized SHAKE128 CTR-XOF based on XKCP.
 */

#ifndef VEXOF_H
#define VEXOF_H

#include "../shake/KeccakHash.h"

#ifndef PARALLELISM
/**
 * Setting PARALLELISM to 8 will not always improve the performance on AVX512 architectures.
 * Whether is does depends on specifics of the application. Testing is advised.
 */
#if __AVX512F__
#define PARALLELISM 8
#elif __AVX2__
#define PARALLELISM 4
#else
#define PARALLELISM 1
#endif
#endif

typedef struct
{
    Keccak_HashInstance keccak_instance;
    uint8_t prepared_state[200 * PARALLELISM];
    ALIGN(PARALLELISM * 8)
    uint8_t states_data[200 * PARALLELISM];
    int squeezing;
    uint64_t block;
    uint64_t index;
} VeXOF_Instance;

/**
 * Function to initialize the VeXOF instance.
 * @param  vexof_instance    Pointer to the VeXOF hash instance to be initialized.
 * @return KECCAK_SUCCESS if successful, KECCAK_FAIL otherwise.
 */
int VeXOF_HashInitialize(VeXOF_Instance *vexof_instance);

/**
 * Function to give input data to be absorbed. Can be called multiple times.
 * @param  vexof_instance    Pointer to the VeXOF instance.
 * @param  data              Pointer to the input data.
 * @param  num_bytes         The number of input bytes provided in the input data.
 * @return KECCAK_SUCCESS if successful, KECCAK_FAIL otherwise.
 */
int VeXOF_HashUpdate(VeXOF_Instance *vexof_instance, const uint8_t *data, size_t num_bytes);

/**
 * Function to squeeze output data. Can be called multiple times.
 * @param  vexof_instance    Pointer to the VeXOF instance.
 * @param  data              Pointer to the buffer where to store the output data.
 * @param  num_bytes         The number of output bytes desired.
 * @return KECCAK_SUCCESS if successful, KECCAK_FAIL otherwise.
 */
int VeXOF_Squeeze(VeXOF_Instance *vexof_instance, uint64_t *data, size_t num_bytes);

/**
 * Function to generate XOF data from a short seed.
 * @param  seed              Pointer to the seed data.
 * @param  input_bytes       The number of seed bytes.
 * @param  output            Pointer to the buffer where to store the output data.
 * @param  output_bytes      The number of output bytes desired.
 */
void vexof(const uint8_t *seed, size_t input_bytes, uint64_t *output, size_t output_bytes);

#endif
