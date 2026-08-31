/**
 * @file xgf16_avx.c
 */
#include "xgf16_avx.h"
#include "../../gf16_core/xgf16.h"
#include "../../gf16_core/gf16.h"

#include <string.h>

__m256i vtl_mt[16];

__m256i mtk2_16[256];

__m256i vtl_basis_t1, vtl_basis_t2, vtl_basis_t4, vtl_basis_t8;
__m256i vtl_basis_m1, vtl_basis_m2, vtl_basis_m4, vtl_basis_m8;
__m256i vtl_basis_zero;
static int vtl_ready = 0;
#define vtl_t1 vtl_basis_t1
#define vtl_t2 vtl_basis_t2
#define vtl_t4 vtl_basis_t4
#define vtl_t8 vtl_basis_t8
#define vtl_m1_mask vtl_basis_m1
#define vtl_m2_mask vtl_basis_m2
#define vtl_m4_mask vtl_basis_m4
#define vtl_m8_mask vtl_basis_m8
#define vtl_zero vtl_basis_zero

static __m256i build_table_for(uint8_t k) {
    uint8_t buf[32];
    for (int i = 0; i < 16; ++i)
        buf[i] = (uint8_t)(gf16_mul((gf16_t)k, (gf16_t)i) & 0x0F);
    memcpy(buf + 16, buf, 16);
    return _mm256_loadu_si256((const __m256i *)buf);
}

void vtl_init(void) {
    if (vtl_ready) return;
    gf16_field_init();
    for (int k = 0; k < 16; ++k) vtl_mt[k] = build_table_for((uint8_t)k);
    vtl_t1 = vtl_mt[1];
    vtl_t2 = vtl_mt[2];
    vtl_t4 = vtl_mt[4];
    vtl_t8 = vtl_mt[8];
    vtl_m1_mask = _mm256_set1_epi32(1);
    vtl_m2_mask = _mm256_set1_epi32(2);
    vtl_m4_mask = _mm256_set1_epi32(4);
    vtl_m8_mask = _mm256_set1_epi32(8);
    vtl_zero = _mm256_setzero_si256();

    for (int i = 0; i < 16; ++i)
        for (int j = 0; j < 16; ++j) {
            uint8_t buf[32];
            for (int k = 0; k < 16; ++k) {
                uint8_t lo = (uint8_t)(gf16_mul((gf16_t)j, (gf16_t)k) & 0x0F);
                uint8_t hi = (uint8_t)(gf16_mul((gf16_t)i, (gf16_t)k) & 0x0F);
                buf[k] = (hi << 4) | lo;
            }
            memcpy(buf + 16, buf, 16);
            mtk2_16[i * 16 + j] = _mm256_loadu_si256((const __m256i *)buf);
        }

    vtl_ready = 1;
}

static __m256i gf16_32_mul_32_core(__m256i a, __m256i b, __m256i acc) {
    __m256i ax[4];
    ax[0] = a;
    ax[1] = _mm256_shuffle_epi8(vtl_mt[2], a);
    ax[2] = _mm256_shuffle_epi8(vtl_mt[4], a);
    ax[3] = _mm256_shuffle_epi8(vtl_mt[8], a);
    static const uint8_t mb[4] = {0x01, 0x02, 0x04, 0x08};
    for (int k = 0; k < 4; ++k) {
        __m256i mask = _mm256_set1_epi8((char)mb[k]);
        __m256i sel  = _mm256_cmpeq_epi8(_mm256_and_si256(b, mask), mask);
        acc = _mm256_xor_si256(acc, _mm256_and_si256(ax[k], sel));
    }
    return acc;
}

__m256i gf16_32_mul_32(__m256i a, __m256i b) {
    return gf16_32_mul_32_core(a, b, _mm256_setzero_si256());
}

__m256i gf16_32_mul_32_add(__m256i a, __m256i b, __m256i acc) {
    return gf16_32_mul_32_core(a, b, acc);
}
