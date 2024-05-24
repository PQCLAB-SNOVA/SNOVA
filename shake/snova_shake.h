/**
 * Vectorized SHAKE128 CTR-XOF based on XKCP.
 */

#ifndef VEXOF_H
#define VEXOF_H

#include "KeccakHash.h"

/**
 * Function to generate XOF data from a short seed.
 * @param  seed              Pointer to the seed data.
 * @param  input_bytes       The number of seed bytes.
 * @param  output            Pointer to the buffer where to store the output data.
 * @param  output_bytes      The number of output bytes desired.
 */
void snova_shake(const uint8_t *seed, size_t input_bytes, uint64_t *output, size_t output_bytes);

#endif
