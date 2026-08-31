// SPDX-License-Identifier: MIT
/*
 * Vectorized SHAKE128 CTR-XOF based on XKCP.
 * Copyright (c) 2024 SNOVA TEAM
 */

#ifndef SNOVA_SHAKE_H
#define SNOVA_SHAKE_H

#include "KeccakHash.h"

/*
 * Function to generate XOF data from a short seed.
 * @param  seed              Pointer to the seed data.
 * @param  input_bytes       The number of seed bytes.
 * @param  output            Pointer to the buffer where to store the output data.
 * @param  output_bytes      The number of output bytes desired.
 */
void snova_shake_ref(const uint8_t *seed, size_t input_bytes, uint64_t *output, size_t output_bytes);
void snova_shake_opt(const uint8_t *seed, size_t input_bytes, uint64_t *output, size_t output_bytes);

#if (defined(SNOVA_VERIFY_STREAM) && SNOVA_VERIFY_STREAM) || \
    (defined(SNOVA_SIGN_STREAM) && SNOVA_SIGN_STREAM) || \
    (defined(SNOVA_KEYGEN_STREAM) && SNOVA_KEYGEN_STREAM) || \
    (defined(SNOVA_PKX_PGEN) && SNOVA_PKX_PGEN)
#define SNOVA_XOF_STREAM_API 1
/*
 streaming verify: persistent-context refactor of snova_shake_squeeze's
 * per-batch loop. Each _next produces one permute batch (PARALLELISM * rate
 * bytes); batches concatenate in block-index order, i.e. byte-identical to
 * the one-shot squeeze by construction (gated by the chunker selftest). 
 */
typedef struct {
    uint64_t prepared64[25 * 8]; /* per-lane padded states (sized for <=8 lanes) */
    uint64_t block;              /* next CTR block index */
    uint32_t byteIOIndex;        /* seed length inside the rate (counter offset) */
    uint32_t bytes_rate;         /* SHAKE128 rate in bytes (168) */
} snova_xof_stream_t;

void snova_shake_stream_init(snova_xof_stream_t *ctx, const uint8_t *seed,
                             size_t input_bytes);
/*
 * Returns the number of bytes produced (PARALLELISM * rate). 
 */
size_t snova_shake_stream_next(snova_xof_stream_t *ctx, uint8_t *out);
#endif /* SNOVA_VERIFY_STREAM */

#endif
