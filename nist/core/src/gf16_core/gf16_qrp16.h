/**
 * @file gf16_core/gf16_qrp16.h
 */
#ifndef SNOVA_GF16_QRP16_H
#define SNOVA_GF16_QRP16_H

#include <stdint.h>

static inline uint64_t gf16_u64_mul(uint64_t a, uint64_t b) {
    uint64_t t = 0, mask;
    for (int i = 0; i < 4; i++) {
        mask = (b & 0x1111111111111111ull) * 0xf;
        t ^= (a & mask);
        mask = ((a >> 3) & 0x1111111111111111ull);
        a = ((a ^ (mask * 0x9)) << 1) ^ mask;
        b >>= 1;
    }
    return t;
}

#if defined(__AVX2__)
#include <immintrin.h>

static inline __m256i gf16_qrp16_256_mul_(__m256i a, __m256i b, __m256i unit) {
    __m256i t = _mm256_setzero_si256();
    for (int i = 0; i < 4; i++) {
        __m256i lb = _mm256_and_si256(b, unit);
        __m256i mask = _mm256_or_si256(
            _mm256_or_si256(lb, _mm256_slli_epi64(lb, 1)),
            _mm256_or_si256(_mm256_slli_epi64(lb, 2), _mm256_slli_epi64(lb, 3)));
        t = _mm256_xor_si256(t, _mm256_and_si256(a, mask));
        __m256i m3 = _mm256_and_si256(_mm256_srli_epi64(a, 3), unit);
        a = _mm256_xor_si256(
            _mm256_slli_epi64(
                _mm256_xor_si256(a, _mm256_or_si256(m3, _mm256_slli_epi64(m3, 3))), 1),
            m3);
        b = _mm256_srli_epi64(b, 1);
    }
    return t;
}

static inline __m256i gf16_qrp16_256_nib_mul(__m256i a, __m256i b) {
    return gf16_qrp16_256_mul_(a, b, _mm256_set1_epi64x(0x1111111111111111ll));
}

static inline __m128i gf16_qrp16_256_nib_outer4_pack(
    __m128i a0, __m128i a1, __m128i a2, __m128i a3,
    __m128i b0, __m128i b1, __m128i b2, __m128i b3) {
    const __m128i A_lo = _mm_or_si128(a0, _mm_slli_epi16(a1, 4));
    const __m128i A_hi = _mm_or_si128(a2, _mm_slli_epi16(a3, 4));
    const __m128i B_lo = _mm_or_si128(b0, _mm_slli_epi16(b1, 4));
    const __m128i B_hi = _mm_or_si128(b2, _mm_slli_epi16(b3, 4));
    const __m256i A = _mm256_set_m128i(A_hi, A_lo);
    const __m256i B = _mm256_set_m128i(B_hi, B_lo);
    const __m256i P = gf16_qrp16_256_nib_mul(A, B);
    return _mm_xor_si128(_mm256_castsi256_si128(P),
                         _mm256_extracti128_si256(P, 1));
}

static inline __m128i gf16_nibpack_fold128(__m128i t) {
    const __m128i m = _mm_set1_epi8(0x0f);
    return _mm_xor_si128(_mm_and_si128(t, m),
                         _mm_and_si128(_mm_srli_epi16(t, 4), m));
}

static inline __m256i gf16_qrp16_256_byte_mul(__m256i a, __m256i b) {
    const __m256i unit = _mm256_set1_epi8(0x01);
    const __m256i zero = _mm256_setzero_si256();
    __m256i t = zero;
    for (int i = 0; i < 4; i++) {
        __m256i lb = _mm256_and_si256(b, unit);
        __m256i mask = _mm256_sub_epi8(zero, lb);
        t = _mm256_xor_si256(t, _mm256_and_si256(a, mask));
        __m256i m3 = _mm256_and_si256(_mm256_srli_epi64(a, 3), unit);
        a = _mm256_xor_si256(
            _mm256_slli_epi64(
                _mm256_xor_si256(a, _mm256_or_si256(m3, _mm256_slli_epi64(m3, 3))), 1),
            m3);
        b = _mm256_srli_epi64(b, 1);
    }
    return t;
}

static inline __m128i gf16_qrp16_128_mul_(__m128i a, __m128i b, __m128i unit) {
    __m128i t = _mm_setzero_si128();
    for (int i = 0; i < 4; i++) {
        __m128i lb = _mm_and_si128(b, unit);
        __m128i mask = _mm_or_si128(
            _mm_or_si128(lb, _mm_slli_epi64(lb, 1)),
            _mm_or_si128(_mm_slli_epi64(lb, 2), _mm_slli_epi64(lb, 3)));
        t = _mm_xor_si128(t, _mm_and_si128(a, mask));
        __m128i m3 = _mm_and_si128(_mm_srli_epi64(a, 3), unit);
        a = _mm_xor_si128(
            _mm_slli_epi64(
                _mm_xor_si128(a, _mm_or_si128(m3, _mm_slli_epi64(m3, 3))), 1),
            m3);
        b = _mm_srli_epi64(b, 1);
    }
    return t;
}

static inline __m128i gf16_qrp16_128_byte_mul(__m128i a, __m128i b) {
    const __m128i unit = _mm_set1_epi8(0x01);
    const __m128i zero = _mm_setzero_si128();
    __m128i t = zero;
    for (int i = 0; i < 4; i++) {
        __m128i lb = _mm_and_si128(b, unit);
        __m128i mask = _mm_sub_epi8(zero, lb);
        t = _mm_xor_si128(t, _mm_and_si128(a, mask));
        __m128i m3 = _mm_and_si128(_mm_srli_epi64(a, 3), unit);
        a = _mm_xor_si128(
            _mm_slli_epi64(
                _mm_xor_si128(a, _mm_or_si128(m3, _mm_slli_epi64(m3, 3))), 1),
            m3);
        b = _mm_srli_epi64(b, 1);
    }
    return t;
}

#endif

#endif
