/**
 * @file platforms/wasm_simd128/xgf16_wasm.c
 */
#include "gf16_simd_wasm.h"

gf16_lane32_t vtl_mt[16];
gf16_lane32_t mtk2_16[256];
static int vtl_inited = 0;

static uint8_t wasm_gf16_mul_scalar(uint8_t a, uint8_t b) {
    uint32_t as = (uint32_t)(a & 0x0Fu);
    as = (as | (as << 4));
    as = (as & 0x41u) | ((as << 2) & 0x208u);
    uint32_t bs = (uint32_t)(b & 0x0Fu);
    bs = (bs | (bs << 4));
    bs = (bs & 0x41u) | ((bs << 2) & 0x208u);
    uint32_t p = as * bs;
    uint32_t res = p & 0x49249249u;
    uint32_t up = p >> 12;
    res = res ^ up ^ (up << 3);
    up = res >> 12;
    res = res ^ up ^ (up << 3);
    up = res >> 12;
    res = res ^ up ^ (up << 3);
    res &= 0x249u;
    res = res | (res >> 4);
    return (uint8_t)((res & 0x5u) | ((res >> 2) & 0xau));
}

void vtl_init(void) {
    if (vtl_inited) return;
    for (int k = 0; k < 16; ++k) {
        uint8_t buf[32];
        for (int i = 0; i < 16; ++i) {
            uint8_t prod = wasm_gf16_mul_scalar((uint8_t)k, (uint8_t)i);
            buf[i]      = prod;
            buf[i + 16] = prod;
        }
        vtl_mt[k] = gf16_lane32_load(buf);
    }
    for (int i_scalar = 0; i_scalar < 16; ++i_scalar) {
        for (int j_scalar = 0; j_scalar < 16; ++j_scalar) {
            int pair = i_scalar | (j_scalar << 4);
            uint8_t buf[32];
            for (int idx = 0; idx < 16; ++idx) {
                uint8_t lo = wasm_gf16_mul_scalar((uint8_t)i_scalar, (uint8_t)idx);
                uint8_t hi = wasm_gf16_mul_scalar((uint8_t)j_scalar, (uint8_t)idx);
                uint8_t v = (uint8_t)(lo | (hi << 4));
                buf[idx]      = v;
                buf[idx + 16] = v;
            }
            mtk2_16[pair] = gf16_lane32_load(buf);
        }
    }
    vtl_inited = 1;
}
