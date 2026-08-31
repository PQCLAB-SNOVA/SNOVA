/**
 * @file gf16.h
 */
#ifndef M0_GF16_H
#define M0_GF16_H

#include <stdint.h>

typedef uint8_t gf16_t;

void gf16_field_init(void);

static inline gf16_t gf16_add(gf16_t a, gf16_t b) { return (gf16_t)(a ^ b); }

gf16_t gf16_mul(gf16_t a, gf16_t b);

gf16_t gf16_inv(gf16_t a);

#endif
