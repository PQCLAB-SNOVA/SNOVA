/**
 * @file aes_ni.c
 */
#include "aes_ni.h"
#include <immintrin.h>
#include <wmmintrin.h>
#include <string.h>

static __m128i aes128_assist(__m128i tmp, __m128i ak) {
    ak = _mm_shuffle_epi32(ak, 0xff);
    __m128i t3 = _mm_slli_si128(tmp, 4);
    tmp = _mm_xor_si128(tmp, t3);
    t3 = _mm_slli_si128(t3, 4);
    tmp = _mm_xor_si128(tmp, t3);
    t3 = _mm_slli_si128(t3, 4);
    tmp = _mm_xor_si128(tmp, t3);
    return _mm_xor_si128(tmp, ak);
}

void aes128_keyexp(const uint8_t key[16], uint8_t rk_out[176]) {
    __m128i rk[11];
    rk[0] = _mm_loadu_si128((const __m128i *)key);
    rk[1]  = aes128_assist(rk[0],  _mm_aeskeygenassist_si128(rk[0],  0x01));
    rk[2]  = aes128_assist(rk[1],  _mm_aeskeygenassist_si128(rk[1],  0x02));
    rk[3]  = aes128_assist(rk[2],  _mm_aeskeygenassist_si128(rk[2],  0x04));
    rk[4]  = aes128_assist(rk[3],  _mm_aeskeygenassist_si128(rk[3],  0x08));
    rk[5]  = aes128_assist(rk[4],  _mm_aeskeygenassist_si128(rk[4],  0x10));
    rk[6]  = aes128_assist(rk[5],  _mm_aeskeygenassist_si128(rk[5],  0x20));
    rk[7]  = aes128_assist(rk[6],  _mm_aeskeygenassist_si128(rk[6],  0x40));
    rk[8]  = aes128_assist(rk[7],  _mm_aeskeygenassist_si128(rk[7],  0x80));
    rk[9]  = aes128_assist(rk[8],  _mm_aeskeygenassist_si128(rk[8],  0x1b));
    rk[10] = aes128_assist(rk[9],  _mm_aeskeygenassist_si128(rk[9],  0x36));
    memcpy(rk_out, rk, 176);
}

void aes128_encrypt(const uint8_t rk_in[176], const uint8_t in[16], uint8_t out[16]) {
    __m128i rk[11];
    memcpy(rk, rk_in, 176);
    __m128i s = _mm_loadu_si128((const __m128i *)in);
    s = _mm_xor_si128(s, rk[0]);
    s = _mm_aesenc_si128(s, rk[1]);
    s = _mm_aesenc_si128(s, rk[2]);
    s = _mm_aesenc_si128(s, rk[3]);
    s = _mm_aesenc_si128(s, rk[4]);
    s = _mm_aesenc_si128(s, rk[5]);
    s = _mm_aesenc_si128(s, rk[6]);
    s = _mm_aesenc_si128(s, rk[7]);
    s = _mm_aesenc_si128(s, rk[8]);
    s = _mm_aesenc_si128(s, rk[9]);
    s = _mm_aesenclast_si128(s, rk[10]);
    _mm_storeu_si128((__m128i *)out, s);
}

static void aes256_a1(__m128i *t1, __m128i *t2) {
    *t2 = _mm_shuffle_epi32(*t2, 0xff);
    __m128i t3 = _mm_slli_si128(*t1, 4);
    *t1 = _mm_xor_si128(*t1, t3);
    t3 = _mm_slli_si128(t3, 4);
    *t1 = _mm_xor_si128(*t1, t3);
    t3 = _mm_slli_si128(t3, 4);
    *t1 = _mm_xor_si128(*t1, t3);
    *t1 = _mm_xor_si128(*t1, *t2);
}

static void aes256_a2(__m128i *t1, __m128i *t3) {
    __m128i t2 = _mm_aeskeygenassist_si128(*t1, 0x00);
    __m128i t4 = _mm_shuffle_epi32(t2, 0xaa);
    __m128i t5 = _mm_slli_si128(*t3, 4);
    *t3 = _mm_xor_si128(*t3, t5);
    t5 = _mm_slli_si128(t5, 4);
    *t3 = _mm_xor_si128(*t3, t5);
    t5 = _mm_slli_si128(t5, 4);
    *t3 = _mm_xor_si128(*t3, t5);
    *t3 = _mm_xor_si128(*t3, t4);
}

void aes256_keyexp(const uint8_t key[32], uint8_t rk_out[240]) {
    __m128i rk[15];
    rk[0] = _mm_loadu_si128((const __m128i *)key);
    rk[1] = _mm_loadu_si128((const __m128i *)(key + 16));
    __m128i t1 = rk[0], t3 = rk[1], t2;
#define R256(i, rc)                                                     \
    do {                                                                \
        t2 = _mm_aeskeygenassist_si128(t3, (rc));                       \
        aes256_a1(&t1, &t2); rk[(i)] = t1;                              \
        aes256_a2(&t1, &t3); rk[(i) + 1] = t3;                          \
    } while (0)
    R256(2,  0x01); R256(4,  0x02); R256(6,  0x04);
    R256(8,  0x08); R256(10, 0x10); R256(12, 0x20);
    t2 = _mm_aeskeygenassist_si128(t3, 0x40);
    aes256_a1(&t1, &t2); rk[14] = t1;
#undef R256
    memcpy(rk_out, rk, 240);
}

void aes256_encrypt(const uint8_t rk_in[240], const uint8_t in[16], uint8_t out[16]) {
    __m128i rk[15];
    memcpy(rk, rk_in, 240);
    __m128i s = _mm_loadu_si128((const __m128i *)in);
    s = _mm_xor_si128(s, rk[0]);
    for (int i = 1; i < 14; ++i) s = _mm_aesenc_si128(s, rk[i]);
    s = _mm_aesenclast_si128(s, rk[14]);
    _mm_storeu_si128((__m128i *)out, s);
}

void aes128_ctr_zero(uint8_t *out, size_t outlen, const uint8_t key[16]) {
    uint8_t rk_bytes[176];
    aes128_keyexp(key, rk_bytes);
    __m128i rk[11];
    memcpy(rk, rk_bytes, 176);

    const __m128i bswap_mask = _mm_set_epi8(
        8, 9, 10, 11, 12, 13, 14, 15,
        7, 6, 5, 4, 3, 2, 1, 0);

    __m128i ctr_vec = _mm_setzero_si128();

    while (outlen >= 64) {
        __m128i ctr_le = _mm_shuffle_epi8(ctr_vec, bswap_mask);
        __m128i n0 = ctr_vec;
        __m128i n1 = _mm_shuffle_epi8(
            _mm_add_epi64(ctr_le, _mm_set_epi64x(1, 0)), bswap_mask);
        __m128i n2 = _mm_shuffle_epi8(
            _mm_add_epi64(ctr_le, _mm_set_epi64x(2, 0)), bswap_mask);
        __m128i n3 = _mm_shuffle_epi8(
            _mm_add_epi64(ctr_le, _mm_set_epi64x(3, 0)), bswap_mask);

        __m128i t0 = _mm_xor_si128(n0, rk[0]);
        __m128i t1 = _mm_xor_si128(n1, rk[0]);
        __m128i t2 = _mm_xor_si128(n2, rk[0]);
        __m128i t3 = _mm_xor_si128(n3, rk[0]);
#define AES_NI_ENC4(IDX) \
        t0 = _mm_aesenc_si128(t0, rk[IDX]); \
        t1 = _mm_aesenc_si128(t1, rk[IDX]); \
        t2 = _mm_aesenc_si128(t2, rk[IDX]); \
        t3 = _mm_aesenc_si128(t3, rk[IDX])
        AES_NI_ENC4(1); AES_NI_ENC4(2); AES_NI_ENC4(3);
        AES_NI_ENC4(4); AES_NI_ENC4(5); AES_NI_ENC4(6);
        AES_NI_ENC4(7); AES_NI_ENC4(8); AES_NI_ENC4(9);
#undef AES_NI_ENC4
        t0 = _mm_aesenclast_si128(t0, rk[10]);
        t1 = _mm_aesenclast_si128(t1, rk[10]);
        t2 = _mm_aesenclast_si128(t2, rk[10]);
        t3 = _mm_aesenclast_si128(t3, rk[10]);

        _mm_storeu_si128((__m128i *)(out +  0), t0);
        _mm_storeu_si128((__m128i *)(out + 16), t1);
        _mm_storeu_si128((__m128i *)(out + 32), t2);
        _mm_storeu_si128((__m128i *)(out + 48), t3);

        ctr_vec = _mm_shuffle_epi8(
            _mm_add_epi64(ctr_le, _mm_set_epi64x(4, 0)), bswap_mask);
        out += 64;
        outlen -= 64;
    }

    uint8_t ctr[16];
    _mm_storeu_si128((__m128i *)ctr, ctr_vec);
    while (outlen > 0) {
        __m128i s = _mm_loadu_si128((const __m128i *)ctr);
        s = _mm_xor_si128(s, rk[0]);
        s = _mm_aesenc_si128(s, rk[1]);
        s = _mm_aesenc_si128(s, rk[2]);
        s = _mm_aesenc_si128(s, rk[3]);
        s = _mm_aesenc_si128(s, rk[4]);
        s = _mm_aesenc_si128(s, rk[5]);
        s = _mm_aesenc_si128(s, rk[6]);
        s = _mm_aesenc_si128(s, rk[7]);
        s = _mm_aesenc_si128(s, rk[8]);
        s = _mm_aesenc_si128(s, rk[9]);
        s = _mm_aesenclast_si128(s, rk[10]);
        if (outlen >= 16) {
            _mm_storeu_si128((__m128i *)out, s);
            out += 16; outlen -= 16;
        } else {
            uint8_t blk[16];
            _mm_storeu_si128((__m128i *)blk, s);
            memcpy(out, blk, outlen);
            outlen = 0;
        }
        for (int i = 15; i >= 0; --i) if (++ctr[i]) break;
    }
}

#if (defined(SNOVA_VERIFY_STREAM) && SNOVA_VERIFY_STREAM) || \
    (defined(SNOVA_SIGN_STREAM) && SNOVA_SIGN_STREAM) || \
    (defined(SNOVA_KEYGEN_STREAM) && SNOVA_KEYGEN_STREAM) || \
    (defined(SNOVA_PKX_PGEN) && SNOVA_PKX_PGEN)
void aes128_ctr_zero_at(uint8_t *out, size_t outlen, const uint8_t key[16],
                        uint64_t block_offset) {
    uint8_t rk_bytes[176];
    aes128_keyexp(key, rk_bytes);
    __m128i rk[11];
    memcpy(rk, rk_bytes, 176);

    const __m128i bswap_mask = _mm_set_epi8(
        8, 9, 10, 11, 12, 13, 14, 15,
        7, 6, 5, 4, 3, 2, 1, 0);

    uint8_t ctr0[16] = {0};
    for (int i = 0; i < 8; ++i) ctr0[15 - i] = (uint8_t)(block_offset >> (8 * i));
    __m128i ctr_vec = _mm_loadu_si128((const __m128i *)ctr0);

    while (outlen >= 64) {
        __m128i ctr_le = _mm_shuffle_epi8(ctr_vec, bswap_mask);
        __m128i n0 = ctr_vec;
        __m128i n1 = _mm_shuffle_epi8(
            _mm_add_epi64(ctr_le, _mm_set_epi64x(1, 0)), bswap_mask);
        __m128i n2 = _mm_shuffle_epi8(
            _mm_add_epi64(ctr_le, _mm_set_epi64x(2, 0)), bswap_mask);
        __m128i n3 = _mm_shuffle_epi8(
            _mm_add_epi64(ctr_le, _mm_set_epi64x(3, 0)), bswap_mask);

        __m128i t0 = _mm_xor_si128(n0, rk[0]);
        __m128i t1 = _mm_xor_si128(n1, rk[0]);
        __m128i t2 = _mm_xor_si128(n2, rk[0]);
        __m128i t3 = _mm_xor_si128(n3, rk[0]);
#define AES_NI_ENC4_AT(IDX) \
        t0 = _mm_aesenc_si128(t0, rk[IDX]); \
        t1 = _mm_aesenc_si128(t1, rk[IDX]); \
        t2 = _mm_aesenc_si128(t2, rk[IDX]); \
        t3 = _mm_aesenc_si128(t3, rk[IDX])
        AES_NI_ENC4_AT(1); AES_NI_ENC4_AT(2); AES_NI_ENC4_AT(3);
        AES_NI_ENC4_AT(4); AES_NI_ENC4_AT(5); AES_NI_ENC4_AT(6);
        AES_NI_ENC4_AT(7); AES_NI_ENC4_AT(8); AES_NI_ENC4_AT(9);
#undef AES_NI_ENC4_AT
        t0 = _mm_aesenclast_si128(t0, rk[10]);
        t1 = _mm_aesenclast_si128(t1, rk[10]);
        t2 = _mm_aesenclast_si128(t2, rk[10]);
        t3 = _mm_aesenclast_si128(t3, rk[10]);

        _mm_storeu_si128((__m128i *)(out +  0), t0);
        _mm_storeu_si128((__m128i *)(out + 16), t1);
        _mm_storeu_si128((__m128i *)(out + 32), t2);
        _mm_storeu_si128((__m128i *)(out + 48), t3);

        ctr_vec = _mm_shuffle_epi8(
            _mm_add_epi64(ctr_le, _mm_set_epi64x(4, 0)), bswap_mask);
        out += 64;
        outlen -= 64;
    }

    uint8_t ctr[16];
    _mm_storeu_si128((__m128i *)ctr, ctr_vec);
    while (outlen > 0) {
        __m128i s = _mm_loadu_si128((const __m128i *)ctr);
        s = _mm_xor_si128(s, rk[0]);
        s = _mm_aesenc_si128(s, rk[1]);
        s = _mm_aesenc_si128(s, rk[2]);
        s = _mm_aesenc_si128(s, rk[3]);
        s = _mm_aesenc_si128(s, rk[4]);
        s = _mm_aesenc_si128(s, rk[5]);
        s = _mm_aesenc_si128(s, rk[6]);
        s = _mm_aesenc_si128(s, rk[7]);
        s = _mm_aesenc_si128(s, rk[8]);
        s = _mm_aesenc_si128(s, rk[9]);
        s = _mm_aesenclast_si128(s, rk[10]);
        if (outlen >= 16) {
            _mm_storeu_si128((__m128i *)out, s);
            out += 16; outlen -= 16;
        } else {
            uint8_t blk[16];
            _mm_storeu_si128((__m128i *)blk, s);
            memcpy(out, blk, outlen);
            outlen = 0;
        }
        for (int i = 15; i >= 0; --i) if (++ctr[i]) break;
    }
}
#endif
