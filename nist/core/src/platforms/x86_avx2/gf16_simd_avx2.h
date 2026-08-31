/**
 * @file platforms/x86_avx2/gf16_simd_avx2.h
 */
#ifndef SNOVA_PLATFORMS_X86_AVX2_GF16_SIMD_AVX2_H
#define SNOVA_PLATFORMS_X86_AVX2_GF16_SIMD_AVX2_H

#include <stdint.h>
#include <immintrin.h>

typedef __m256i gf16_lane32_t;
typedef __m128i gf16_lane16_t;

#define gf16_lane32_load(p)        _mm256_load_si256((const __m256i *)(p))
#define gf16_lane32_loadu(p)       _mm256_loadu_si256((const __m256i *)(p))
#define gf16_lane32_store(p, v)    _mm256_store_si256((__m256i *)(p), (v))
#define gf16_lane32_storeu(p, v)   _mm256_storeu_si256((__m256i *)(p), (v))

#define gf16_lane32_zero()         _mm256_setzero_si256()
#define gf16_lane32_set1(b)        _mm256_set1_epi8((char)(b))
#define gf16_lane32_set1_epi32(u)  _mm256_set1_epi32((int)(u))

#define gf16_lane32_xor(a, b)      _mm256_xor_si256((a), (b))
#define gf16_lane32_and(a, b)      _mm256_and_si256((a), (b))
#define gf16_lane32_or(a, b)       _mm256_or_si256((a), (b))

#define gf16_lane32_slli16(v, n)   _mm256_slli_epi16((v), (n))
#define gf16_lane32_srli16(v, n)   _mm256_srli_epi16((v), (n))

#define gf16_lane32_shuffle(t, i)  _mm256_shuffle_epi8((t), (i))

#define gf16_lane32_cmpeq8(a, b)   _mm256_cmpeq_epi8((a), (b))

static inline uint8_t gf16_lane32_xor_reduce_nibble(gf16_lane32_t v) {
    __m128i lo = _mm256_castsi256_si128(v);
    __m128i hi = _mm256_extracti128_si256(v, 1);
    __m128i x = _mm_xor_si128(lo, hi);
    x = _mm_xor_si128(x, _mm_srli_si128(x, 8));
    x = _mm_xor_si128(x, _mm_srli_si128(x, 4));
    x = _mm_xor_si128(x, _mm_srli_si128(x, 2));
    x = _mm_xor_si128(x, _mm_srli_si128(x, 1));
    return (uint8_t)(_mm_extract_epi8(x, 0) & 0x0F);
}

#if defined(__GFNI__)
#define SNOVA_GF16_HAS_GFNI 1
#define gf16_lane32_gf2p8mul(a, b) _mm256_gf2p8mul_epi8((a), (b))
static inline gf16_lane32_t gf16_lane32_gfni_cleanup(gf16_lane32_t v) {
    const __m256i m0f = _mm256_set1_epi8(0x0f);
    __m256i vhi = _mm256_and_si256(v, _mm256_set1_epi8((char)0xf0));
    __m256i a = _mm256_srli_epi16(vhi, 3);
    __m256i b = _mm256_srli_epi16(v, 4);
    return _mm256_and_si256(_mm256_xor_si256(_mm256_xor_si256(v, a), b), m0f);
}
#endif

#if defined(__GFNI__)
#define SNOVA_GF16_HAS_CELLMUL 1
#define SNOVA_CELL_MUL256(a, b) _mm256_gf2p8mul_epi8((a), (b))
#define SNOVA_CELL_MUL128(a, b) _mm_gf2p8mul_epi8((a), (b))
#elif defined(SNOVA_QRP16_CELLMUL)
#include "../../gf16_core/gf16_qrp16.h"
#define SNOVA_GF16_HAS_CELLMUL 1
#define SNOVA_CELL_MUL256(a, b) gf16_qrp16_256_byte_mul((a), (b))
#define SNOVA_CELL_MUL128(a, b) gf16_qrp16_128_byte_mul((a), (b))
static inline gf16_lane32_t gf16_lane32_gfni_cleanup(gf16_lane32_t v) {
    const __m256i m0f = _mm256_set1_epi8(0x0f);
    __m256i vhi = _mm256_and_si256(v, _mm256_set1_epi8((char)0xf0));
    __m256i a = _mm256_srli_epi16(vhi, 3);
    __m256i b = _mm256_srli_epi16(v, 4);
    return _mm256_and_si256(_mm256_xor_si256(_mm256_xor_si256(v, a), b), m0f);
}
#endif

#endif
