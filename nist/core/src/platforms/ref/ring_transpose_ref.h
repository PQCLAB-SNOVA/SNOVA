/**
 * @file platforms/ref/ring_transpose_ref.h
 */
#ifndef SNOVA_PLATFORMS_REF_RING_TRANSPOSE_REF_H
#define SNOVA_PLATFORMS_REF_RING_TRANSPOSE_REF_H

#include <stdint.h>

static inline void ring_transpose_L4_chunk(
    gf16_lane32_t v0, gf16_lane32_t v1, gf16_lane32_t v2, gf16_lane32_t v3,
    gf16_lane32_t *o0, gf16_lane32_t *o1, gf16_lane32_t *o2, gf16_lane32_t *o3)
{
    gf16_lane32_t r0, r1, r2, r3;
    const gf16_lane32_t *vs[4] = { &v0, &v1, &v2, &v3 };
    gf16_lane32_t *rs[4] = { &r0, &r1, &r2, &r3 };
    for (int cell = 0; cell < 8; ++cell)
        for (int ei = 0; ei < 4; ++ei)
            for (int ej = 0; ej < 4; ++ej)
                rs[ei]->b[cell * 4 + ej] = vs[ej]->b[cell * 4 + ei];
    *o0 = r0;
    *o1 = r1;
    *o2 = r2;
    *o3 = r3;
}

#endif
