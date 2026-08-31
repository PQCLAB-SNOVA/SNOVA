/**
 * @file gf16_core/gf16_qrp16_ct.c
 */
#include <stdint.h>
#include <stdio.h>
#include "gf16_qrp16.h"
#include <valgrind/memcheck.h>

int main(void) {
    volatile uint64_t sink64 = 0;
    uint64_t a = 0, b = 0;
    VALGRIND_MAKE_MEM_UNDEFINED(&a, sizeof a);
    VALGRIND_MAKE_MEM_UNDEFINED(&b, sizeof b);
    sink64 ^= gf16_u64_mul(a, b);

#if defined(__AVX2__)
    uint8_t abuf[32], bbuf[32];
    VALGRIND_MAKE_MEM_UNDEFINED(abuf, sizeof abuf);
    VALGRIND_MAKE_MEM_UNDEFINED(bbuf, sizeof bbuf);
    __m256i av = _mm256_loadu_si256((const __m256i *)abuf);
    __m256i bv = _mm256_loadu_si256((const __m256i *)bbuf);
    uint8_t out[32];
    _mm256_storeu_si256((__m256i *)out, gf16_qrp16_256_nib_mul(av, bv));
    for (int i = 0; i < 32; i++) sink64 ^= out[i];
    _mm256_storeu_si256((__m256i *)out, gf16_qrp16_256_byte_mul(av, bv));
    for (int i = 0; i < 32; i++) sink64 ^= out[i];
#endif

    VALGRIND_MAKE_MEM_DEFINED((void *)&sink64, sizeof sink64);
    printf("gf16_qrp16_ct: ran 3 variants on tainted operands "
           "(CT verdict = memcheck error count; sink=%llu)\n",
           (unsigned long long)sink64);
    return 0;
}
