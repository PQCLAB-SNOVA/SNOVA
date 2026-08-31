/**
 * @file platforms/portable_opt/xgf16_opt.c
 */
#include "gf16_simd_opt.h"

gf16_lane32_t vtl_mt[16];
gf16_lane32_t mtk2_16[256];
static int vtl_inited = 0;

static uint8_t opt_gf16_mul_scalar(uint8_t a, uint8_t b) {
    return snova_opt_unspread(
        snova_opt_spread(a & 0x0F) * snova_opt_spread(b & 0x0F));
}

void vtl_init(void) {
    if (vtl_inited) return;
    for (int k = 0; k < 16; ++k) {
        for (int i = 0; i < 16; ++i) {
            uint8_t prod = opt_gf16_mul_scalar((uint8_t)k, (uint8_t)i);
            vtl_mt[k].b[i]      = prod;
            vtl_mt[k].b[i + 16] = prod;
        }
    }
    for (int i_scalar = 0; i_scalar < 16; ++i_scalar) {
        for (int j_scalar = 0; j_scalar < 16; ++j_scalar) {
            int pair = i_scalar | (j_scalar << 4);
            for (int idx = 0; idx < 16; ++idx) {
                uint8_t lo = opt_gf16_mul_scalar((uint8_t)i_scalar, (uint8_t)idx);
                uint8_t hi = opt_gf16_mul_scalar((uint8_t)j_scalar, (uint8_t)idx);
                mtk2_16[pair].b[idx]      = (uint8_t)(lo | (hi << 4));
                mtk2_16[pair].b[idx + 16] = mtk2_16[pair].b[idx];
            }
        }
    }
    vtl_inited = 1;
}
