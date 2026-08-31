/**
 * @file gf16m_mul.c
 */
#include "gf16m.h"
#include "xgf16.h"
#include <string.h>

void gf16m_mul(const gf16m_t a, const gf16m_t b, gf16m_t c) {
    uint32_t xa[SNOVA_SQ_RANK], xb[SNOVA_SQ_RANK];
    for (int k = 0; k < SNOVA_SQ_RANK; ++k) {
        xa[k] = xgf16_spread(a[k]);
        xb[k] = xgf16_spread(b[k]);
    }

    gf16m_t t;
    for (int i = 0; i < SNOVA_RANK; ++i)
        for (int j = 0; j < SNOVA_RANK; ++j) {
            uint32_t acc = 0;
            for (int k = 0; k < SNOVA_RANK; ++k)
                acc ^= xa[i * SNOVA_RANK + k] * xb[k * SNOVA_RANK + j];
            GF16M_AT(t, i, j) = xgf16_unspread(xgf16_reduce(acc));
        }
    memcpy(c, t, SNOVA_SQ_RANK);
}
