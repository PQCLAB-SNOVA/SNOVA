/**
 * @file xgf16.h
 */
#ifndef M1_XGF16_H
#define M1_XGF16_H

#include <stdint.h>

#define XGF16_LANE_MASK 0x49249249u
#define XGF16_NIBBLE_MASK 0x249u

static inline uint32_t xgf16_spread(uint8_t a) {
    uint32_t m = (uint32_t)(a & 0x0F) | ((uint32_t)(a & 0x0F) << 4);
    return (m & 0x41u) | ((m << 2) & 0x208u);
}

static inline uint32_t xgf16_reduce(uint32_t v) {
    uint32_t r = v & XGF16_LANE_MASK;
    uint32_t u = r >> 12;
    r = r ^ u ^ (u << 3);
    u = r >> 12;  r = r ^ u ^ (u << 3);
    u = r >> 12;  r = r ^ u ^ (u << 3);
    return r & XGF16_NIBBLE_MASK;
}

static inline uint8_t xgf16_unspread(uint32_t v) {
    uint32_t r = v | (v >> 4);
    return (uint8_t)((r & 0x5u) | ((r >> 2) & 0xAu));
}

#endif
