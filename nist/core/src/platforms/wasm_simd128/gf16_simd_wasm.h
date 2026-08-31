/**
 * @file platforms/wasm_simd128/gf16_simd_wasm.h
 */
#ifndef SNOVA_PLATFORMS_WASM_SIMD128_GF16_SIMD_WASM_H
#define SNOVA_PLATFORMS_WASM_SIMD128_GF16_SIMD_WASM_H

#include <stdint.h>
#include <string.h>
#include <wasm_simd128.h>

typedef union {
    struct {
        v128_t v0;
        v128_t v1;
    };
    uint8_t b[32];
} gf16_lane32_t;

typedef struct {
    v128_t v;
} gf16_lane16_t;

static inline gf16_lane32_t gf16_lane32_load(const void *p) {
    gf16_lane32_t r;
    r.v0 = wasm_v128_load((const v128_t *)p);
    r.v1 = wasm_v128_load((const v128_t *)((const uint8_t *)p + 16));
    return r;
}
#define gf16_lane32_loadu(p)        gf16_lane32_load(p)

static inline void gf16_lane32_store(void *p, gf16_lane32_t v) {
    wasm_v128_store((v128_t *)p, v.v0);
    wasm_v128_store((v128_t *)((uint8_t *)p + 16), v.v1);
}
#define gf16_lane32_storeu(p, v)    gf16_lane32_store((p), (v))

static inline gf16_lane32_t gf16_lane32_zero(void) {
    gf16_lane32_t r;
    r.v0 = wasm_i64x2_const(0, 0);
    r.v1 = wasm_i64x2_const(0, 0);
    return r;
}
static inline gf16_lane32_t gf16_lane32_set1(uint8_t v) {
    gf16_lane32_t r;
    r.v0 = wasm_i8x16_splat((int8_t)v);
    r.v1 = r.v0;
    return r;
}
static inline gf16_lane32_t gf16_lane32_set1_epi32(uint32_t u) {
    gf16_lane32_t r;
    r.v0 = wasm_i32x4_splat((int32_t)u);
    r.v1 = r.v0;
    return r;
}

static inline gf16_lane32_t gf16_lane32_xor(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    r.v0 = wasm_v128_xor(a.v0, b.v0);
    r.v1 = wasm_v128_xor(a.v1, b.v1);
    return r;
}
static inline gf16_lane32_t gf16_lane32_and(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    r.v0 = wasm_v128_and(a.v0, b.v0);
    r.v1 = wasm_v128_and(a.v1, b.v1);
    return r;
}
static inline gf16_lane32_t gf16_lane32_or(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    r.v0 = wasm_v128_or(a.v0, b.v0);
    r.v1 = wasm_v128_or(a.v1, b.v1);
    return r;
}

static inline gf16_lane32_t gf16_lane32_slli16(gf16_lane32_t v, int n) {
    gf16_lane32_t r;
    r.v0 = wasm_i16x8_shl(v.v0, n);
    r.v1 = wasm_i16x8_shl(v.v1, n);
    return r;
}
static inline gf16_lane32_t gf16_lane32_srli16(gf16_lane32_t v, int n) {
    gf16_lane32_t r;
    r.v0 = wasm_u16x8_shr(v.v0, n);
    r.v1 = wasm_u16x8_shr(v.v1, n);
    return r;
}

static inline gf16_lane32_t gf16_lane32_shuffle(gf16_lane32_t tab, gf16_lane32_t idx) {
    gf16_lane32_t r;
    const v128_t nib = wasm_i8x16_splat((int8_t)0x0F);
    v128_t lo0 = wasm_v128_and(idx.v0, nib);
    v128_t lo1 = wasm_v128_and(idx.v1, nib);
    v128_t look0 = wasm_i8x16_swizzle(tab.v0, lo0);
    v128_t look1 = wasm_i8x16_swizzle(tab.v1, lo1);
    v128_t m0 = wasm_i8x16_shr(idx.v0, 7);
    v128_t m1 = wasm_i8x16_shr(idx.v1, 7);
    r.v0 = wasm_v128_andnot(look0, m0);
    r.v1 = wasm_v128_andnot(look1, m1);
    return r;
}

static inline gf16_lane32_t gf16_lane32_cmpeq8(gf16_lane32_t a, gf16_lane32_t b) {
    gf16_lane32_t r;
    r.v0 = wasm_i8x16_eq(a.v0, b.v0);
    r.v1 = wasm_i8x16_eq(a.v1, b.v1);
    return r;
}

static inline uint8_t gf16_lane32_xor_reduce_nibble(gf16_lane32_t v) {
    v128_t x = wasm_v128_xor(v.v0, v.v1);
    uint64_t lo = wasm_i64x2_extract_lane(x, 0);
    uint64_t hi = wasm_i64x2_extract_lane(x, 1);
    uint64_t s = lo ^ hi;
    s ^= s >> 32;
    s ^= s >> 16;
    s ^= s >> 8;
    return (uint8_t)(s & 0x0Fu);
}

static inline v128_t snova_wasm_spread_x4(v128_t n) {
    v128_t m = wasm_v128_or(n, wasm_i32x4_shl(n, 4));
    v128_t lo = wasm_v128_and(m, wasm_i32x4_splat(0x41));
    v128_t hi = wasm_v128_and(wasm_i32x4_shl(m, 2), wasm_i32x4_splat(0x208));
    return wasm_v128_or(lo, hi);
}

static inline v128_t snova_wasm_reduce_x4(v128_t p) {
    v128_t mask_l = wasm_i32x4_splat(0x49249249);
    v128_t res = wasm_v128_and(p, mask_l);
    v128_t upper = wasm_u32x4_shr(p, 12);
    res = wasm_v128_xor(wasm_v128_xor(res, upper), wasm_i32x4_shl(upper, 3));
    upper = wasm_u32x4_shr(res, 12);
    res = wasm_v128_xor(wasm_v128_xor(res, upper), wasm_i32x4_shl(upper, 3));
    upper = wasm_u32x4_shr(res, 12);
    res = wasm_v128_xor(wasm_v128_xor(res, upper), wasm_i32x4_shl(upper, 3));
    return wasm_v128_and(res, wasm_i32x4_splat(0x249));
}

static inline v128_t snova_wasm_unspread_x4(v128_t p) {
    v128_t r = snova_wasm_reduce_x4(p);
    r = wasm_v128_or(r, wasm_u32x4_shr(r, 4));
    v128_t mask_5 = wasm_i32x4_splat(0x5);
    v128_t mask_a = wasm_i32x4_splat(0xa);
    return wasm_v128_or(wasm_v128_and(r, mask_5),
                        wasm_v128_and(wasm_u32x4_shr(r, 2), mask_a));
}

static inline void snova_wasm_spread_16(v128_t bytes16, v128_t out[4]) {
    v128_t lo16 = wasm_u16x8_extend_low_u8x16(bytes16);
    v128_t hi16 = wasm_u16x8_extend_high_u8x16(bytes16);
    v128_t b0_3 = wasm_u32x4_extend_low_u16x8(lo16);
    v128_t b4_7 = wasm_u32x4_extend_high_u16x8(lo16);
    v128_t b8_11 = wasm_u32x4_extend_low_u16x8(hi16);
    v128_t b12_15 = wasm_u32x4_extend_high_u16x8(hi16);
    v128_t mask = wasm_i32x4_splat(0x0F);
    out[0] = snova_wasm_spread_x4(wasm_v128_and(b0_3, mask));
    out[1] = snova_wasm_spread_x4(wasm_v128_and(b4_7, mask));
    out[2] = snova_wasm_spread_x4(wasm_v128_and(b8_11, mask));
    out[3] = snova_wasm_spread_x4(wasm_v128_and(b12_15, mask));
}

static inline v128_t snova_wasm_pack_16(v128_t in[4]) {
    v128_t u0 = snova_wasm_unspread_x4(in[0]);
    v128_t u1 = snova_wasm_unspread_x4(in[1]);
    v128_t u2 = snova_wasm_unspread_x4(in[2]);
    v128_t u3 = snova_wasm_unspread_x4(in[3]);
    v128_t p01 = wasm_u16x8_narrow_i32x4(u0, u1);
    v128_t p23 = wasm_u16x8_narrow_i32x4(u2, u3);
    return wasm_u8x16_narrow_i16x8(p01, p23);
}

static inline gf16_lane32_t gf16_32_mul_32(gf16_lane32_t a, gf16_lane32_t b) {
    v128_t as[4], bs[4];
    snova_wasm_spread_16(a.v0, as);
    snova_wasm_spread_16(b.v0, bs);
    v128_t lo[4];
    for (int k = 0; k < 4; ++k) lo[k] = wasm_i32x4_mul(as[k], bs[k]);
    snova_wasm_spread_16(a.v1, as);
    snova_wasm_spread_16(b.v1, bs);
    v128_t hi[4];
    for (int k = 0; k < 4; ++k) hi[k] = wasm_i32x4_mul(as[k], bs[k]);
    gf16_lane32_t r;
    r.v0 = snova_wasm_pack_16(lo);
    r.v1 = snova_wasm_pack_16(hi);
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
    uint8_t buf[32];
    for (int i = 0; i < 16; ++i) {
        uint32_t ks = (uint32_t)(k & 0x0F);
        ks = (ks | (ks << 4));
        ks = (ks & 0x41u) | ((ks << 2) & 0x208u);
        uint32_t is = (uint32_t)i;
        is = (is | (is << 4));
        is = (is & 0x41u) | ((is << 2) & 0x208u);
        uint32_t p = ks * is;
        uint32_t res = p & 0x49249249u;
        uint32_t up = p >> 12;
        res = res ^ up ^ (up << 3);
        up = res >> 12;
        res = res ^ up ^ (up << 3);
        up = res >> 12;
        res = res ^ up ^ (up << 3);
        res = (res & 0x249u);
        res = res | (res >> 4);
        uint8_t prod = (uint8_t)((res & 0x5u) | ((res >> 2) & 0xau));
        buf[i] = prod;
        buf[i + 16] = prod;
    }
    r = gf16_lane32_load(buf);
    return r;
}

static inline gf16_lane32_t vtl_ct_multtab_pair(uint8_t k_lo, uint8_t k_hi) {
    gf16_lane32_t lo = vtl_ct_multtab(k_lo);
    gf16_lane32_t hi = vtl_ct_multtab(k_hi);
    v128_t mask = wasm_i8x16_splat(0x0F);
    gf16_lane32_t r;
    r.v0 = wasm_v128_or(wasm_v128_and(lo.v0, mask),
                         wasm_i16x8_shl(wasm_v128_and(hi.v0, mask), 4));
    r.v1 = wasm_v128_or(wasm_v128_and(lo.v1, mask),
                         wasm_i16x8_shl(wasm_v128_and(hi.v1, mask), 4));
    return r;
}

static inline gf16_lane32_t gf16_32_mul_k(gf16_lane32_t a, uint8_t k) {
    return gf16_lane32_shuffle(vtl_mt[k & 0x0F], a);
}

static inline gf16_lane32_t gf16_32_mul_k_add(gf16_lane32_t a, uint8_t k, gf16_lane32_t acc) {
    return gf16_lane32_xor(acc, gf16_32_mul_k(a, k));
}

#endif
