/**
 * @file snova_core/pgen.h
 */
#ifndef SNOVA_PGEN_H
#define SNOVA_PGEN_H

#include <stdint.h>
#include <string.h>

#if SNOVA_PK_EXPAND_SHAKE
#include "../third_party/shake_xkcp/snova_shake.h"
#else
#include "../platforms/x86_avx2/aes_ni.h"
#endif

typedef struct {
    uint8_t     block;
    uint16_t    mi, ni;
    uint16_t    ncols;
    const gf_t *cells;
} snova_prow_t;

#define SNOVA_PGEN_MAXCOLS   (SNOVA_v > SNOVA_o ? SNOVA_v : SNOVA_o)
#if SNOVA_q == 16
#define SNOVA_PGEN_ROWBYTES  ((SNOVA_PGEN_MAXCOLS * SNOVA_l2 + 1) / 2 + 1)
#else
#define SNOVA_PGEN_ROWBYTES  (SNOVA_PGEN_MAXCOLS * SNOVA_l2)
#endif
#define SNOVA_PGEN_CHUNK     1344

#define SNOVA_PGEN_ROWS_P11  (SNOVA_m1 * SNOVA_v)
#define SNOVA_PGEN_ROWS_P12  (SNOVA_m1 * SNOVA_v)
#define SNOVA_PGEN_ROWS_P21  (SNOVA_m1 * SNOVA_o)
#define SNOVA_PGEN_ROWS_ALL  (SNOVA_PGEN_ROWS_P11 + SNOVA_PGEN_ROWS_P12 + SNOVA_PGEN_ROWS_P21)

typedef struct {
#if SNOVA_PK_EXPAND_SHAKE
    snova_xof_stream_t sh;
#else
    uint8_t  key[16];
    uint64_t next_block;
#endif
    _Alignas(32) uint8_t buf[SNOVA_PGEN_CHUNK];
    size_t   pos, avail;
    int      seq;
    _Alignas(32) uint8_t packed[SNOVA_PGEN_ROWBYTES];
    _Alignas(32) gf_t    cells[SNOVA_PGEN_MAXCOLS * SNOVA_l2 + 16];
#if SNOVA_q == 16
    uint8_t  carry_nib;
    uint8_t  has_carry;
#endif
} snova_pgen_t;

static inline void snova_pgen_read(snova_pgen_t *g, uint8_t *dst, size_t n) {
    while (n) {
        if (g->pos == g->avail) {
#if SNOVA_PK_EXPAND_SHAKE
            g->avail = snova_shake_stream_next(&g->sh, g->buf);
#else
            aes128_ctr_zero_at(g->buf, SNOVA_PGEN_CHUNK, g->key, g->next_block);
            g->next_block += SNOVA_PGEN_CHUNK / 16;
            g->avail = SNOVA_PGEN_CHUNK;
#endif
            g->pos = 0;
        }
        size_t take = g->avail - g->pos;
        if (take > n) take = n;
        memcpy(dst, g->buf + g->pos, take);
        g->pos += take;
        dst += take;
        n -= take;
    }
}

static inline void snova_pgen_init(snova_pgen_t *g, const uint8_t pkseed[16]) {
#if SNOVA_PK_EXPAND_SHAKE
    snova_shake_stream_init(&g->sh, pkseed, 16);
#else
    memcpy(g->key, pkseed, 16);
    g->next_block = 0;
#endif
    g->pos = 0;
    g->avail = 0;
    g->seq = 0;
    memset(g->cells + SNOVA_PGEN_MAXCOLS * SNOVA_l2, 0, 16 * sizeof(gf_t));
#if SNOVA_q == 16
    g->has_carry = 0;
    g->carry_nib = 0;
#endif
}

static inline int snova_pgen_peek(const snova_pgen_t *g, int *block, int *mi, int *ni, int *ncols) {
    int s = g->seq;
    if (s >= SNOVA_PGEN_ROWS_ALL) return 0;
    if (s < SNOVA_PGEN_ROWS_P11) {
        *block = 0; *mi = s / SNOVA_v; *ni = s % SNOVA_v; *ncols = SNOVA_v;
    } else if (s < SNOVA_PGEN_ROWS_P11 + SNOVA_PGEN_ROWS_P12) {
        int t = s - SNOVA_PGEN_ROWS_P11;
        *block = 1; *mi = t / SNOVA_v; *ni = t % SNOVA_v; *ncols = SNOVA_o;
    } else {
        int t = s - SNOVA_PGEN_ROWS_P11 - SNOVA_PGEN_ROWS_P12;
        *block = 2; *mi = t / SNOVA_o; *ni = t % SNOVA_o; *ncols = SNOVA_v;
    }
    return 1;
}

static inline int snova_pgen_next_row_into(snova_pgen_t *g, snova_prow_t *out, gf_t *dst) {
    int block, mi, ni, ncols;
    if (!snova_pgen_peek(g, &block, &mi, &ni, &ncols)) return 0;
    size_t ngf = (size_t)ncols * SNOVA_l2;
#if SNOVA_q == 16
    size_t i = 0;
    if (g->has_carry) { dst[i++] = g->carry_nib; g->has_carry = 0; }
    size_t rem = ngf - i;
    size_t nbytes = rem / 2 + (rem & 1);
    snova_pgen_read(g, g->packed, nbytes);
    size_t pairs = rem / 2;
    for (size_t b = 0; b < pairs; ++b) {
        dst[i++] = (gf_t)(g->packed[b] & 0x0f);
        dst[i++] = (gf_t)(g->packed[b] >> 4);
    }
    if (rem & 1) {
        dst[i++] = (gf_t)(g->packed[pairs] & 0x0f);
        g->carry_nib = (uint8_t)(g->packed[pairs] >> 4);
        g->has_carry = 1;
    }
#else
    snova_pgen_read(g, g->packed, ngf);
    convert_bytes_to_GF(dst, g->packed, ngf);
#endif
    out->block = (uint8_t)block;
    out->mi = (uint16_t)mi;
    out->ni = (uint16_t)ni;
    out->ncols = (uint16_t)ncols;
    out->cells = dst;
    g->seq += 1;
    return 1;
}

static inline int snova_pgen_next_row(snova_pgen_t *g, snova_prow_t *out) {
    return snova_pgen_next_row_into(g, out, g->cells);
}

static inline void snova_pgen_fill_pblocks(const uint8_t pkseed[16], gf_t *P_matrix) {
    snova_pgen_t g;
    snova_pgen_init(&g, pkseed);
    const size_t p12_base = (size_t)SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
    const size_t p21_base = (size_t)SNOVA_m1 * SNOVA_v * (SNOVA_v + SNOVA_o) * SNOVA_l2;
    snova_prow_t r;
    while (snova_pgen_next_row(&g, &r)) {
        size_t off;
        if (r.block == 0)      off = ((size_t)(r.mi * SNOVA_v + r.ni) * SNOVA_v) * SNOVA_l2;
        else if (r.block == 1) off = p12_base + ((size_t)(r.mi * SNOVA_v + r.ni) * SNOVA_o) * SNOVA_l2;
        else                   off = p21_base + ((size_t)(r.mi * SNOVA_o + r.ni) * SNOVA_v) * SNOVA_l2;
        memcpy(P_matrix + off, r.cells, (size_t)r.ncols * SNOVA_l2 * sizeof(gf_t));
    }
}

static inline void snova_pgen_fill_p11(const uint8_t pkseed[16], gf_t *P11) {
    snova_pgen_t g;
    snova_pgen_init(&g, pkseed);
    snova_prow_t r;
    while (g.seq < SNOVA_PGEN_ROWS_P11 && snova_pgen_next_row(&g, &r)) {
        size_t off = ((size_t)(r.mi * SNOVA_v + r.ni) * SNOVA_v) * SNOVA_l2;
        memcpy(P11 + off, r.cells, (size_t)r.ncols * SNOVA_l2 * sizeof(gf_t));
    }
}

#endif
