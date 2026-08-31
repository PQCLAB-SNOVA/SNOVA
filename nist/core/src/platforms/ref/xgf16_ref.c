/**
 * @file platforms/ref/xgf16_ref.c
 */
#include "gf16_simd_ref.h"
#include "../../gf16_core/gf16.h"

gf16_lane32_t vtl_mt[16];
gf16_lane32_t mtk2_16[256];
static int vtl_inited = 0;

void vtl_init(void) {
    if (vtl_inited) return;
    for (int k = 0; k < 16; ++k) {
        for (int i = 0; i < 16; ++i) {
            uint8_t prod = gf16_mul((gf16_t)k, (gf16_t)i);
            vtl_mt[k].b[i]      = prod;
            vtl_mt[k].b[i + 16] = prod;
        }
    }
    for (int i_scalar = 0; i_scalar < 16; ++i_scalar) {
        for (int j_scalar = 0; j_scalar < 16; ++j_scalar) {
            int pair = i_scalar | (j_scalar << 4);
            for (int idx = 0; idx < 16; ++idx) {
                uint8_t i_lo = (uint8_t)idx;
                uint8_t lo = gf16_mul((gf16_t)i_scalar, (gf16_t)i_lo);
                uint8_t hi = gf16_mul((gf16_t)j_scalar, (gf16_t)i_lo);
                mtk2_16[pair].b[idx]      = (uint8_t)(lo | (hi << 4));
                mtk2_16[pair].b[idx + 16] = mtk2_16[pair].b[idx];
            }
        }
    }
    vtl_inited = 1;
}

gf16_lane32_t gf16_32_mul_32(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    for (int i = 0; i < 32; ++i)
        r.b[i] = gf16_mul((gf16_t)(a.b[i] & 0x0F), (gf16_t)(b.b[i] & 0x0F));
    return r;
}

gf16_lane32_t gf16_32_mul_32_add(gf16_lane32_t a, gf16_lane32_t b, gf16_lane32_t acc) {
    gf16_lane32_t prod = gf16_32_mul_32(a, b);
    return gf16_lane32_xor(acc, prod);
}
