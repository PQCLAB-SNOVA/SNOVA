/**
 * @file gf16_core/gf16_mullo16.h
 */
#ifndef SNOVA_GF16_CORE_GF16_MULLO16_H
#define SNOVA_GF16_CORE_GF16_MULLO16_H

#include <stdint.h>
#include <immintrin.h>

static inline __m256i cl_gf16_expand_u16x16(__m256i x) {
    __m256i v = _mm256_or_si256(_mm256_or_si256(x, _mm256_slli_epi16(x, 3)),
                                _mm256_or_si256(_mm256_slli_epi16(x, 6), _mm256_slli_epi16(x, 9)));
    return _mm256_and_si256(v, _mm256_set1_epi16(0x1111));
}
static inline __m256i cl_gf16_compress_u16x16(__m256i a) {
    const __m256i m0f = _mm256_set1_epi16(0x000f);
    __m256i val = _mm256_xor_si256(
        _mm256_xor_si256(_mm256_and_si256(a, m0f),
                         _mm256_srli_epi16(_mm256_and_si256(a, _mm256_set1_epi16(0x00f0)), 3)),
        _mm256_xor_si256(_mm256_srli_epi16(_mm256_and_si256(a, _mm256_set1_epi16(0x0f00)), 6),
                         _mm256_srli_epi16(_mm256_and_si256(a, _mm256_set1_epi16((short)0xf000)), 9)));
    val = _mm256_xor_si256(_mm256_xor_si256(val,
              _mm256_srli_epi16(_mm256_and_si256(val, _mm256_set1_epi16(0x00f0)), 3)),
              _mm256_srli_epi16(val, 4));
    return _mm256_and_si256(val, m0f);
}
static inline __m128i cl_gf16_pack_u16_to_bytes(__m256i c) {
    return _mm_packus_epi16(_mm256_castsi256_si128(c), _mm256_extracti128_si256(c, 1));
}

static inline uint16_t cl_expand_scalar16(uint8_t x) {
    return (uint16_t)((x & 1) | ((x & 2) << 3) | ((x & 4) << 6) | ((x & 8) << 9));
}

#endif
