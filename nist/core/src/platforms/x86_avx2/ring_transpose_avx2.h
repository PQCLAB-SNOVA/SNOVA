/**
 * @file ring_transpose.h
 */
#ifndef M7_RING_TRANSPOSE_H
#define M7_RING_TRANSPOSE_H

#include <immintrin.h>

static inline void ring_transpose_L4_chunk(
    __m256i v0, __m256i v1, __m256i v2, __m256i v3,
    __m256i *o0, __m256i *o1, __m256i *o2, __m256i *o3)
{
    __m256i t01l = _mm256_unpacklo_epi8(v0, v1);
    __m256i t01h = _mm256_unpackhi_epi8(v0, v1);
    __m256i t23l = _mm256_unpacklo_epi8(v2, v3);
    __m256i t23h = _mm256_unpackhi_epi8(v2, v3);

    __m256i z0 = _mm256_unpacklo_epi16(t01l, t23l);
    __m256i z1 = _mm256_unpackhi_epi16(t01l, t23l);
    __m256i z2 = _mm256_unpacklo_epi16(t01h, t23h);
    __m256i z3 = _mm256_unpackhi_epi16(t01h, t23h);

    __m256i a01l = _mm256_unpacklo_epi32(z0, z1);
    __m256i a01h = _mm256_unpackhi_epi32(z0, z1);
    __m256i a23l = _mm256_unpacklo_epi32(z2, z3);
    __m256i a23h = _mm256_unpackhi_epi32(z2, z3);

    *o0 = _mm256_unpacklo_epi64(a01l, a23l);
    *o1 = _mm256_unpackhi_epi64(a01l, a23l);
    *o2 = _mm256_unpacklo_epi64(a01h, a23h);
    *o3 = _mm256_unpackhi_epi64(a01h, a23h);
}

#endif
