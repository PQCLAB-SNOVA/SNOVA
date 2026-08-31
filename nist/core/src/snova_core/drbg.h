/**
 * @file drbg.h
 */
#ifndef SNOVA_DRBG_H
#define SNOVA_DRBG_H

#include <stddef.h>
#include <stdint.h>

/**
 * @param entropy_input 48-byte entropy.
 * @param personalization Personalization string or NULL.
 * @param security_strength Strength (256, kept for API compatibility).
 * @return none.
 */
void randombytes_init(const uint8_t *entropy_input,
                      const uint8_t *personalization, int security_strength);

/**
 * @param x Output buffer.
 * @param xlen Number of bytes.
 * @return 0 on success.
 */
int randombytes(uint8_t *x, size_t xlen);

#endif
