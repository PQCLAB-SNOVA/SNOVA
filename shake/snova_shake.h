// SPDX-License-Identifier: MIT
/**
 * Vectorized SHAKE128 CTR-XOF based on XKCP.
 *
 * Copyright (c) 2024 SNOVA TEAM
 */

#ifndef SNOVA_SHAKE_H
#define SNOVA_SHAKE_H

#include "KeccakHash.h"

/**
 * Function to generate XOF data from a short seed.
 * @param  seed              Pointer to the seed data.
 * @param  input_bytes       The number of seed bytes.
 * @param  output            Pointer to the buffer where to store the output data.
 * @param  output_bytes      The number of output bytes desired.
 */
void snova_shake_ref(const uint8_t *seed, size_t input_bytes, uint64_t *output, size_t output_bytes);
void snova_shake_opt(const uint8_t *seed, size_t input_bytes, uint64_t *output, size_t output_bytes);

#endif
