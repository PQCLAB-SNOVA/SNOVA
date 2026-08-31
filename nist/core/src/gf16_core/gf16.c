/**
 * @file gf16.c
 */
#include "gf16.h"
#include "xgf16.h"

void gf16_field_init(void) {
}

gf16_t gf16_mul(gf16_t a, gf16_t b) {
    uint32_t p = xgf16_spread(a) * xgf16_spread(b);
    return xgf16_unspread(xgf16_reduce(p));
}

gf16_t gf16_inv(gf16_t a) {
    a &= 0x0F;
    uint32_t x  = xgf16_spread(a);
    uint32_t x2 = xgf16_reduce(x * x);
    uint32_t x4 = xgf16_reduce(x2 * x2);
    uint32_t x8 = xgf16_reduce(x4 * x4);
    uint32_t x12 = xgf16_reduce(x4 * x8);
    uint32_t x14 = xgf16_reduce(x2 * x12);
    return xgf16_unspread(x14);
}
