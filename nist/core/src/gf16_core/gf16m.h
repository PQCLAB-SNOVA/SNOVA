/**
 * @file gf16m.h
 */
#ifndef M0_GF16M_H
#define M0_GF16M_H

#include "../snova_params.h"
#include "gf16.h"

typedef gf16_t gf16m_t[SNOVA_SQ_RANK];

#define GF16M_AT(m, i, j) ((m)[(i) * SNOVA_RANK + (j)])

void gf16m_ring_init(void);

void gf16m_zero(gf16m_t a);
void gf16m_copy(gf16m_t dst, const gf16m_t src);
int  gf16m_eq(const gf16m_t a, const gf16m_t b);
void gf16m_identity(gf16m_t a);
void gf16m_add(const gf16m_t a, const gf16m_t b, gf16m_t c);
void gf16m_mul(const gf16m_t a, const gf16m_t b, gf16m_t c);
void gf16m_scale(const gf16m_t a, gf16_t k, gf16m_t c);
void gf16m_transpose(const gf16m_t a, gf16m_t at);
gf16_t gf16m_det(const gf16m_t a);

const gf16_t *gf16m_S(void);
const gf16_t *gf16m_Spow(int k);

void gf16m_make_invertible(gf16m_t a);

#endif
