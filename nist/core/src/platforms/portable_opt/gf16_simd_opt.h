/**
 * @file platforms/portable_opt/gf16_simd_opt.h
 */
#ifndef SNOVA_PLATFORMS_PORTABLE_OPT_GF16_SIMD_OPT_H
#define SNOVA_PLATFORMS_PORTABLE_OPT_GF16_SIMD_OPT_H

#include <stdint.h>
#include <string.h>

typedef struct { uint8_t b[32]; } gf16_lane32_t;
typedef struct { uint8_t b[16]; } gf16_lane16_t;

static inline gf16_lane32_t gf16_lane32_load(const void *p) {
    gf16_lane32_t r;
    memcpy(r.b, p, 32);
    return r;
}
#define gf16_lane32_loadu(p)        gf16_lane32_load(p)

static inline void gf16_lane32_store(void *p, gf16_lane32_t v) {
    memcpy(p, v.b, 32);
}
#define gf16_lane32_storeu(p, v)    gf16_lane32_store((p), (v))

static inline gf16_lane32_t gf16_lane32_zero(void) {
    gf16_lane32_t r;
    memset(r.b, 0, 32);
    return r;
}
static inline gf16_lane32_t gf16_lane32_set1(uint8_t v) {
    gf16_lane32_t r;
    memset(r.b, v, 32);
    return r;
}
static inline gf16_lane32_t gf16_lane32_set1_epi32(uint32_t u) {
    gf16_lane32_t r;
    for (int i = 0; i < 8; ++i) memcpy(&r.b[i * 4], &u, 4);
    return r;
}

static inline gf16_lane32_t gf16_lane32_xor(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    for (int i = 0; i < 32; ++i) r.b[i] = a.b[i] ^ b.b[i];
    return r;
}
static inline gf16_lane32_t gf16_lane32_and(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    for (int i = 0; i < 32; ++i) r.b[i] = a.b[i] & b.b[i];
    return r;
}
static inline gf16_lane32_t gf16_lane32_or(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    for (int i = 0; i < 32; ++i) r.b[i] = a.b[i] | b.b[i];
    return r;
}

static inline gf16_lane32_t gf16_lane32_slli16(gf16_lane32_t v, int n) {
    gf16_lane32_t r;
    for (int i = 0; i < 16; ++i) {
        uint16_t w;
        memcpy(&w, &v.b[i * 2], 2);
        w = (uint16_t)(w << n);
        memcpy(&r.b[i * 2], &w, 2);
    }
    return r;
}
static inline gf16_lane32_t gf16_lane32_srli16(gf16_lane32_t v, int n) {
    gf16_lane32_t r;
    for (int i = 0; i < 16; ++i) {
        uint16_t w;
        memcpy(&w, &v.b[i * 2], 2);
        w = (uint16_t)(w >> n);
        memcpy(&r.b[i * 2], &w, 2);
    }
    return r;
}

static inline gf16_lane32_t gf16_lane32_shuffle(gf16_lane32_t tab, gf16_lane32_t idx) {
    gf16_lane32_t r;
    for (int i = 0; i < 32; ++i) {
        uint8_t v = idx.b[i];
        int lane = i & 16;
        uint8_t lookup = tab.b[lane + (v & 0x0F)];
        uint8_t mask = (uint8_t)((v >> 7) - 1u);
        r.b[i] = lookup & mask;
    }
    return r;
}

static inline gf16_lane32_t gf16_lane32_cmpeq8(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    for (int i = 0; i < 32; ++i) r.b[i] = (a.b[i] == b.b[i]) ? 0xFF : 0x00;
    return r;
}

static inline uint8_t gf16_lane32_xor_reduce_nibble(gf16_lane32_t v) {
    uint8_t x = 0;
    for (int i = 0; i < 32; ++i) x ^= v.b[i];
    return (uint8_t)(x & 0x0F);
}

static inline uint32_t snova_opt_spread(uint8_t n) {
    uint32_t middle = (uint32_t)n | ((uint32_t)n << 4);
    return (middle & 0x41u) | ((middle << 2) & 0x208u);
}

static inline uint32_t snova_opt_reduce(uint32_t p) {
    uint32_t res = p & 0x49249249u;
    uint32_t upper = p >> 12;
    res = res ^ upper ^ (upper << 3);
    upper = res >> 12;
    res = res ^ upper ^ (upper << 3);
    upper = res >> 12;
    res = res ^ upper ^ (upper << 3);
    return res & 0x249u;
}

static inline uint8_t snova_opt_unspread(uint32_t p) {
    uint32_t r = snova_opt_reduce(p);
    r = r | (r >> 4);
    return (uint8_t)((r & 0x5u) | ((r >> 2) & 0xau));
}

static inline gf16_lane32_t gf16_32_mul_32(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    for (int i = 0; i < 32; ++i) {
        uint32_t as = snova_opt_spread(a.b[i] & 0x0F);
        uint32_t bs = snova_opt_spread(b.b[i] & 0x0F);
        r.b[i] = snova_opt_unspread(as * bs);
    }
    return r;
}

static inline gf16_lane32_t gf16_32_mul_32_add(gf16_lane32_t a, gf16_lane32_t b, gf16_lane32_t acc) {
    return gf16_lane32_xor(acc, gf16_32_mul_32(a, b));
}

extern gf16_lane32_t vtl_mt[16];
extern gf16_lane32_t mtk2_16[256];
void vtl_init(void);

static inline gf16_lane32_t vtl_ct_multtab(uint8_t k) {
    gf16_lane32_t r;
    uint32_t ks = snova_opt_spread(k & 0x0F);
    for (int i = 0; i < 16; ++i) {
        uint32_t is = snova_opt_spread((uint8_t)i);
        uint8_t prod = snova_opt_unspread(ks * is);
        r.b[i]      = prod;
        r.b[i + 16] = prod;
    }
    return r;
}

static inline gf16_lane32_t vtl_ct_multtab_pair(uint8_t k_lo, uint8_t k_hi) {
    gf16_lane32_t lo = vtl_ct_multtab(k_lo);
    gf16_lane32_t hi = vtl_ct_multtab(k_hi);
    gf16_lane32_t r;
    for (int i = 0; i < 32; ++i)
        r.b[i] = (uint8_t)((lo.b[i] & 0x0F) | ((hi.b[i] & 0x0F) << 4));
    return r;
}

static inline gf16_lane32_t gf16_32_mul_k(gf16_lane32_t a, uint8_t k) {
    gf16_lane32_t r;
    uint32_t ks = snova_opt_spread(k & 0x0F);
    for (int i = 0; i < 32; ++i) {
        uint32_t as = snova_opt_spread(a.b[i] & 0x0F);
        r.b[i] = snova_opt_unspread(as * ks);
    }
    return r;
}

static inline gf16_lane32_t gf16_32_mul_k_add(gf16_lane32_t a, uint8_t k, gf16_lane32_t acc) {
    return gf16_lane32_xor(acc, gf16_32_mul_k(a, k));
}

#endif
