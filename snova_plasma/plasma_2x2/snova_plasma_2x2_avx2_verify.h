#ifndef PLASMA_2x2_AVX2_VERIFY_H
#define PLASMA_2x2_AVX2_VERIFY_H

#include <immintrin.h>

// Align Public and Intermediate to 256 bits
#define n_SNOVA16 ((n_SNOVA + 15) / 16)
#define n_SNOVA2 (n_SNOVA16 * 16)
#define m_SNOVA32 (((m_SNOVA + 31) / 32) * 32)

/**
 * Parallel Operation Transposed Matrix
 */
static inline void gf16m_transpose_u16_16way(__m256i* ap, const __m256i a) {
    __m256i t0 = _mm256_slli_epi64(a, 4) & _mm256_set1_epi16(0x0f00);
    __m256i t1 = _mm256_srli_epi64(a, 4) & _mm256_set1_epi16(0x00f0);
    uint16_t mask = 0xf00f;
    *ap = (a & _mm256_set1_epi16(mask)) ^ t0 ^ t1;
}

void evaluation_2x2_avx2_vtl(map_group1_u16* restrict map1_u16, uint16_t* restrict sig_var_u16, uint16_t* restrict P22,
                               uint16_t* restrict hash_in_GF16Matrix_u16) {
    // ------- u16 new -------

    uint16_t Left_NL[n_SNOVA_mult4][lsq_SNOVA] __attribute__((aligned(32))) = {0};
    uint16_t Right_NL[n_SNOVA_mult4][lsq_SNOVA] __attribute__((aligned(32))) = {0};
    uint16_t Public[n_SNOVA * n_SNOVA2] __attribute__((aligned(32)));

    uint8_t Intermediate8[lsq_SNOVA * m_SNOVA_mult2 * n_SNOVA2 * l_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t Inter_shuffled[lsq_SNOVA * n_SNOVA * l_SNOVA * m_SNOVA32] __attribute__((aligned(32))) = {0};
    uint8_t hash_in_GF16Matrix0[2 * m_SNOVA32] __attribute__((aligned(32))) = {0};
    uint8_t hash_in_GF16Matrix1[2 * m_SNOVA32] __attribute__((aligned(32))) = {0};

    __m256i *Intermediate256 = (__m256i *)Intermediate8;
    __m256i *Public256 = (__m256i *)Public;

    __m256i l_mask = _mm256_set1_epi64x(0x0f0f0f0f0f0f0f0full);
    __m256i mask_0 = _mm256_set1_epi64x(0xf00ff00ff00ff00full);
    __m256i mask_1 = _mm256_set1_epi64x(0x0f000f000f000f00ull);

    // Convert

    __m256i Aalpha = _mm256_setr_epi16(
        map1_u16->Aalpha[0], map1_u16->Aalpha[1], map1_u16->Aalpha[2], map1_u16->Aalpha[3],
        map1_u16->Aalpha[0], map1_u16->Aalpha[1], map1_u16->Aalpha[2], map1_u16->Aalpha[3],
        map1_u16->Aalpha[0], map1_u16->Aalpha[1], map1_u16->Aalpha[2], map1_u16->Aalpha[3],
        map1_u16->Aalpha[0], map1_u16->Aalpha[1], map1_u16->Aalpha[2], map1_u16->Aalpha[3]);
    __m256i Balpha = _mm256_setr_epi16(
        map1_u16->Balpha[0], map1_u16->Balpha[1], map1_u16->Balpha[2], map1_u16->Balpha[3],
        map1_u16->Balpha[0], map1_u16->Balpha[1], map1_u16->Balpha[2], map1_u16->Balpha[3],
        map1_u16->Balpha[0], map1_u16->Balpha[1], map1_u16->Balpha[2], map1_u16->Balpha[3],
        map1_u16->Balpha[0], map1_u16->Balpha[1], map1_u16->Balpha[2], map1_u16->Balpha[3]);
    __m256i Qalpha1 = _mm256_setr_epi16(
        map1_u16->Qalpha1[0], map1_u16->Qalpha1[1], map1_u16->Qalpha1[2], map1_u16->Qalpha1[3],
        map1_u16->Qalpha1[0], map1_u16->Qalpha1[1], map1_u16->Qalpha1[2], map1_u16->Qalpha1[3],
        map1_u16->Qalpha1[0], map1_u16->Qalpha1[1], map1_u16->Qalpha1[2], map1_u16->Qalpha1[3],
        map1_u16->Qalpha1[0], map1_u16->Qalpha1[1], map1_u16->Qalpha1[2], map1_u16->Qalpha1[3]);
    __m256i Qalpha2 = _mm256_setr_epi16(
        map1_u16->Qalpha2[0], map1_u16->Qalpha2[1], map1_u16->Qalpha2[2], map1_u16->Qalpha2[3],
        map1_u16->Qalpha2[0], map1_u16->Qalpha2[1], map1_u16->Qalpha2[2], map1_u16->Qalpha2[3],
        map1_u16->Qalpha2[0], map1_u16->Qalpha2[1], map1_u16->Qalpha2[2], map1_u16->Qalpha2[3],
        map1_u16->Qalpha2[0], map1_u16->Qalpha2[1], map1_u16->Qalpha2[2], map1_u16->Qalpha2[3]);

    // Parallel Operation Transposed Matrix
    uint16_t sig_var_tr_u16[n_SNOVA_mult16] __attribute__((aligned(32))) = {0};
    __m256i* sig_var_256 = (__m256i*)sig_var_u16;
    __m256i* sig_var_tr_256 = (__m256i*)sig_var_tr_u16;

    for (int index = 0; index < n_SNOVA16; ++index) {
        gf16m_transpose_u16_16way(sig_var_tr_256 + index, sig_var_256[index]);
    }
   
    for (int index = 0; index < n_SNOVA; index += 4)
    {
        __m256i t256;
        gf16m_u64_mul_2x2_16way_4_16((uint64_t *)(sig_var_tr_u16 + index), &Qalpha1, &t256);
        gf16m_u64_mul_2x2_16way_16_16(&Aalpha, &t256, (__m256i *)Left_NL[index]);

        gf16m_u64_mul_2x2_16way_4_16((uint64_t *)(sig_var_u16 + index), &Balpha, &t256);
        gf16m_u64_mul_2x2_16way_16_16(&Qalpha2, &t256, (__m256i *)Right_NL[index]);
    }

    uint16_t *left16 = (uint16_t *)Left_NL;
    for (int idx = 0; idx < lsq_SNOVA * n_SNOVA_mult4; idx++)
    {
        uint16_t val = left16[idx];
        left16[idx] = (val & 0xf00f) ^ ((val >> 4) & 0x00f0) ^ ((val << 4) & 0x0f00);
    }

    // Prepare Public
    for (int mi = 0; mi < m_SNOVA; ++mi) {
        for (int dj = 0; dj < v_SNOVA; dj++)
            for (int dk = 0; dk < v_SNOVA; dk += 16)
                _mm256_storeu_si256((__m256i *)&Public[dj * n_SNOVA2 + dk],
                    _mm256_loadu_si256((__m256i *)&map1_u16->P11[mi][dj][dk]));

        for (int dj = 0; dj < v_SNOVA; dj++)
            for (int dk = v_SNOVA; dk < n_SNOVA; dk++)
                Public[dj * n_SNOVA2 + dk] = map1_u16->P12[mi][dj][dk - v_SNOVA];

        for (int dj = v_SNOVA; dj < n_SNOVA; dj++)
            for (int dk = 0; dk < v_SNOVA; dk += 16)
                _mm256_storeu_si256((__m256i *)&Public[dj * n_SNOVA2 + dk],
                    _mm256_loadu_si256((__m256i *)&map1_u16->P21[mi][dj - v_SNOVA][dk]));

        for (int dj = v_SNOVA; dj < n_SNOVA; dj++)
            for (int dk = v_SNOVA; dk < n_SNOVA; dk++)
                Public[dj * n_SNOVA2 + dk] = P22[((mi * o_SNOVA + dj - v_SNOVA) * o_SNOVA + dk - v_SNOVA)];

        // Transpose at 2x2 level
        for (int idx = 0; idx < n_SNOVA * n_SNOVA16; idx++)
        {
            __m256i pub_0 = Public256[idx];
            Public256[idx] = (pub_0 & mask_0) ^ (_mm256_slli_epi64(pub_0, 4) & mask_1) ^ _mm256_srli_epi64(pub_0 & mask_1, 4);
        }

        // Main loop
        for (int dj = 0; dj < n_SNOVA; ++dj)
            for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                for (int dk_j1 = 0; dk_j1 < n_SNOVA16; dk_j1 ++)
                {
                    __m256i k_lh0 = mtk2_16[Left_NL[dj][alpha] & 0xff];
                    __m256i k_lh1 = mtk2_16[(Left_NL[dj][alpha] >> 8) & 0xff];
                    __m256i pubval = *(__m256i *)&Public[dj * n_SNOVA2 + 16 * dk_j1];

                    Intermediate256[mi * lsq_SNOVA * n_SNOVA16 + alpha * n_SNOVA16 + dk_j1] ^=
                        _mm256_shuffle_epi8(k_lh0, pubval & l_mask) ^
                        _mm256_shuffle_epi8(k_lh1, _mm256_srli_epi16(pubval, 4) & l_mask);
                }
    }

    for (int dk = 0; dk < n_SNOVA; ++dk)
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int mi = 0; mi < m_SNOVA; mi += 2)
            {
                uint16_t val = *(uint16_t *)&Intermediate8[mi * lsq_SNOVA  * n_SNOVA2 * l_SNOVA + alpha * n_SNOVA2 * l_SNOVA + dk * l_SNOVA];
                uint16_t val2 = *(uint16_t *)&Intermediate8[(mi + 1) * lsq_SNOVA  * n_SNOVA2 * l_SNOVA + alpha * n_SNOVA2 * l_SNOVA + dk * l_SNOVA];
                *(uint16_t *)&Inter_shuffled[(alpha * n_SNOVA + dk)* l_SNOVA * m_SNOVA32 + 0 * m_SNOVA32 + mi] = (val & 0xff) ^ ((val2 << 8) & 0xff00);
                *(uint16_t *)&Inter_shuffled[(alpha * n_SNOVA + dk)* l_SNOVA * m_SNOVA32 + 1 * m_SNOVA32 + mi] = ((val >> 8) & 0xff) ^ (val2 & 0xff00);
            }

    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int dk = 0; dk < n_SNOVA; ++dk)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                for (int mi = 0; mi < m_SNOVA; mi += 32)
                {
                    __m256i k_lh0 = mtk2_16[(Right_NL[dk][alpha] >> 8 * j1) & 0xff];
                    __m256i intval0 = *(__m256i *)&Inter_shuffled[(alpha * n_SNOVA + dk)* l_SNOVA * m_SNOVA32 + j1 * m_SNOVA32 + mi];

                    *(__m256i *)&hash_in_GF16Matrix0[mi] ^= _mm256_shuffle_epi8(k_lh0, intval0 & l_mask);
                    *(__m256i *)&hash_in_GF16Matrix1[mi] ^= _mm256_shuffle_epi8(k_lh0, _mm256_srli_epi64(intval0, 4) & l_mask);
                }

    for (int mi = 0; mi < m_SNOVA; ++mi)
        hash_in_GF16Matrix_u16[mi] = hash_in_GF16Matrix0[mi] ^ (hash_in_GF16Matrix1[mi] << 8);
}

int verify_signture_2x2(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature, const uint8_t* pk) {
     uint8_t signed_hash[bytes_hash];
    const uint8_t* pt_salt = pt_signature + bytes_signature;
    public_key* pk_stru = (public_key*)pk;

    // ------- u16 new -------
    map_group1_u16 map1_u16;
    uint16_t hash_in_GF16Matrix_u16[m_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t* hash_in_bytes = (uint8_t*)hash_in_GF16Matrix_u16;
    uint16_t sig_var_u16[n_SNOVA_mult16] __attribute__((aligned(32)));

    createSignedHash(pt_digest, bytes_digest, pk_stru->pt_public_key_seed, pt_salt, signed_hash);

    memcpy(sig_var_u16, pt_signature, bytes_signature); // without salt

#if (o_SNOVA * l_SNOVA) & 0x1 == 1
    signed_hash[bytes_hash - 1] &= 0x0f;
#endif

    // generate PRNG part of public key
    gen_A_B_Q_P_2x2(&map1_u16, pk_stru->pt_public_key_seed);

    // evaluate signature GF16Matrix array
    evaluation_2x2_avx2_vtl(&map1_u16, sig_var_u16, (uint16_t*)(pk_stru->P22), hash_in_GF16Matrix_u16);

    int result = 0;
    for (int i = 0; i < bytes_hash; ++i) {
        if (hash_in_bytes[i] != signed_hash[i]) {
            result = -1;
            break;
        }
    }
    return result;
}

#endif
