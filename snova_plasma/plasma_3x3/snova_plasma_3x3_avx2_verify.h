/**
 * AVX2 optimized version using Vectorized Table Lookup for all multiplications.
 */

#ifndef PLASMA_3x3_AVX2_VERIFY_H
#define PLASMA_3x3_AVX2_VERIFY_H

#define nl_SNOVA32 ((n_SNOVA * l_SNOVA + 31) / 32)
#define nl_SNOVA (nl_SNOVA32 * 32)

#define ml_SNOVA32 ((m_SNOVA * l_SNOVA + 31) / 32)
#define ml_SNOVA (ml_SNOVA32 * 32)

#define mnl_SNOVA32 (m_SNOVA * nl_SNOVA32)
#define mnl_SNOVA (mnl_SNOVA32 * 32)

int verify_signture_vtl(const uint8_t *pt_digest, uint64_t bytes_digest, const uint8_t *pt_signature, const uint8_t *pk)
{
    uint8_t hash_in_bytes[bytes_hash];
    uint8_t signed_hash[bytes_hash];
    const uint8_t *pt_salt = pt_signature + bytes_signature;

    gf16m_t signature_in_GF16Matrix[n_SNOVA + 24] = {0};
    gf16m_t P22[m_SNOVA][o_SNOVA][o_SNOVA];
    map_group1 map1;
    public_key *pk_stru = (public_key *)pk;

    Keccak_HashInstance hashInstance;
    Keccak_HashInitialize_SHAKE256(&hashInstance);
    Keccak_HashUpdate(&hashInstance, pk_stru->pt_public_key_seed, 8 * seed_length_public);
    Keccak_HashUpdate(&hashInstance, pt_digest, 8 * bytes_digest);
    Keccak_HashUpdate(&hashInstance, pt_salt, 8 * bytes_salt);
    Keccak_HashFinal(&hashInstance, NULL);
    Keccak_HashSqueeze(&hashInstance, signed_hash, 8 * bytes_hash);

#if (o_SNOVA * l_SNOVA) & 0x1 == 1
    signed_hash[bytes_hash - 1] &= 0x0f;
#endif

    // generate PRNG part of public key
    gen_A_B_Q_P(&map1, pk_stru->pt_public_key_seed);

    // read  P22
    input_P22((uint8_t *)P22, pk_stru->P22);

    convert_bytes_to_GF16s(pt_signature, (gf16_t *)signature_in_GF16Matrix, GF16s_signature);

    // Align  256 bits
    uint8_t left8[(lsq_SNOVA * l_SNOVA + 1) * nl_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t right8[(lsq_SNOVA * l_SNOVA + 1) * nl_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t signature8[l_SNOVA * nl_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t inter_a8[lsq_SNOVA * l_SNOVA * nl_SNOVA] __attribute__((aligned(32)));
    uint8_t inter_b8[lsq_SNOVA * l_SNOVA * nl_SNOVA] __attribute__((aligned(32))) = {0};

    __m256i *signature256 = (__m256i *)signature8;
    __m256i *inter_a256 = (__m256i *)inter_a8;
    __m256i *inter_b256 = (__m256i *)inter_b8;
    __m256i *left256 = (__m256i *)left8;
    __m256i *right256 = (__m256i *)right8;

    uint8_t public8[l_SNOVA * nl_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t inter8[(lsq_SNOVA * l_SNOVA + 1) / 2 * mnl_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t r_inter8[(l_SNOVA + 1) * lsq_SNOVA * ml_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t res8[l_SNOVA * ml_SNOVA] __attribute__((aligned(32))) = {0};

    __m256i *public256 = (__m256i *)public8;
    __m256i *inter256 = (__m256i *)inter8;
    __m256i *r_inter256 = (__m256i *)r_inter8;
    __m256i *res256 = (__m256i *)res8;

    uint8_t finalres[m_SNOVA * lsq_SNOVA];

    // Prepare signature
    for (int index = 0; index < n_SNOVA; ++index)
        for (int k1 = 0; k1 < rank; k1++)
            for (int i1 = 0; i1 < rank; i1++)
                signature8[k1 * nl_SNOVA + index * l_SNOVA + i1] = signature_in_GF16Matrix[index][k1 * l_SNOVA + i1];

    // First left mult
    memset(inter_a8, 0, sizeof(inter_a8));
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int k1 = 0; k1 < rank; k1++)
            for (int j1 = 0; j1 < rank; j1++)
                for (int index_i1 = 0; index_i1 < nl_SNOVA32; index_i1++)
                {
                    __m256i k_lh0 = mtk2_16[map1.Qalpha1[alpha][k1 * l_SNOVA + j1]];
                    inter_a256[(alpha * l_SNOVA + j1) * nl_SNOVA32 + index_i1] ^=
                        _mm256_shuffle_epi8(k_lh0, signature256[k1 * nl_SNOVA32 + index_i1]);
                }

    // Shuffle intermediate
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int index = 0; index < n_SNOVA; ++index)
            for (int j1 = 0; j1 < rank; j1++)
                for (int i1 = 0; i1 < rank; i1++)
                    inter_b8[(alpha * l_SNOVA + j1) * nl_SNOVA + index * l_SNOVA + i1] = inter_a8[(alpha * l_SNOVA + i1) * nl_SNOVA + index * l_SNOVA + j1];

    // Second left mult
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int k1 = 0; k1 < rank; k1++)
            for (int j1 = 0; j1 < rank; j1++)
                for (int index_i1 = 0; index_i1 < nl_SNOVA32; index_i1++)
                {
                    __m256i k_lh0 = mtk2_16[map1.Aalpha[alpha][k1 * l_SNOVA + j1]];
                    left256[(alpha * l_SNOVA + k1) * nl_SNOVA32 + index_i1] ^=
                        _mm256_shuffle_epi8(k_lh0, inter_b256[(alpha * l_SNOVA + j1) * nl_SNOVA32 + index_i1]);
                }

    for (int alpha_i1 = 0; alpha_i1 < (lsq_SNOVA * l_SNOVA + 1) / 2; alpha_i1++)
        for (int dj_k1 = 0; dj_k1 < nl_SNOVA; dj_k1++)
            left8[alpha_i1 * n_SNOVA * l_SNOVA + dj_k1] = left8[2 * alpha_i1 * nl_SNOVA + dj_k1] ^
                                                                    (left8[(2 * alpha_i1 + 1) * nl_SNOVA + dj_k1] << 4);

    // First right mult
    memset(inter_a8, 0, sizeof(inter_a8));
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int k1 = 0; k1 < rank; k1++)
            for (int j1 = 0; j1 < rank; j1++)
                for (int index_i1 = 0; index_i1 < nl_SNOVA32; index_i1++)
                {
                    __m256i k_lh0 = mtk2_16[map1.Qalpha2[alpha][j1 * l_SNOVA + k1]];
                    inter_a256[(alpha * l_SNOVA + j1) * nl_SNOVA32 + index_i1] ^=
                        _mm256_shuffle_epi8(k_lh0, signature256[k1 * nl_SNOVA32 + index_i1]);
              }

    // Shuffle intermediate
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int index = 0; index < n_SNOVA; ++index)
            for (int j1 = 0; j1 < rank; j1++)
                for (int i1 = 0; i1 < rank; i1++)
                    inter_b8[(alpha * l_SNOVA + j1) * nl_SNOVA + index * l_SNOVA + i1] = inter_a8[(alpha * l_SNOVA + i1) * nl_SNOVA + index * l_SNOVA + j1];

    // Second right mult
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int k1 = 0; k1 < rank; k1++)
            for (int j1 = 0; j1 < rank; j1++)
                for (int index_i1 = 0; index_i1 < nl_SNOVA32; index_i1++)
                {
                    __m256i k_lh0 = mtk2_16[map1.Balpha[alpha][j1 * l_SNOVA + k1]];
                    right256[(alpha * l_SNOVA + k1) * nl_SNOVA32 + index_i1] ^=
                        _mm256_shuffle_epi8(k_lh0, inter_b256[(alpha * l_SNOVA + j1) * nl_SNOVA32 + index_i1]);
                }

    // Main loop
    for (int mi = 0; mi < m_SNOVA; ++mi)
    {
        for (int dj = 0; dj < v_SNOVA; dj++) {
            for (int dk = 0; dk < v_SNOVA; dk++)
                for (int i1 = 0; i1 < rank; i1++)
                    for (int j1 = 0; j1 < rank; j1++)
                        public8[i1 * nl_SNOVA + dk * l_SNOVA + j1] = map1.P11[mi][dj][dk][i1 * l_SNOVA + j1];

            for (int dk = v_SNOVA; dk < n_SNOVA; dk++)
                for (int i1 = 0; i1 < rank; i1++)
                    for (int j1 = 0; j1 < rank; j1++)
                        public8[i1 * nl_SNOVA + dk * l_SNOVA + j1] = map1.P12[mi][dj][dk - v_SNOVA][i1 * l_SNOVA + j1];

            for (int alpha_i1 = 0; alpha_i1 < (lsq_SNOVA * l_SNOVA + 1) / 2; alpha_i1++)
                for (int dk_j1_32 = 0; dk_j1_32 < nl_SNOVA32; dk_j1_32++)
                    for (int k1 = 0; k1 < l_SNOVA; k1++)
                    {
                        __m256i k_lh0 = mtk2_16[left8[alpha_i1 * n_SNOVA * l_SNOVA + dj * l_SNOVA + k1]];
                        __m256i pubval256 = public256[k1 * nl_SNOVA32 + dk_j1_32];
                        inter256[alpha_i1 * mnl_SNOVA32 + mi * nl_SNOVA32 + dk_j1_32] ^= _mm256_shuffle_epi8(k_lh0, pubval256);
                    }
        }

        for (int dj = v_SNOVA; dj < n_SNOVA; dj++) {
            for (int dk = 0; dk < v_SNOVA; dk++)
                for (int i1 = 0; i1 < rank; i1++)
                    for (int j1 = 0; j1 < rank; j1++)
                        public8[i1 * nl_SNOVA + dk * l_SNOVA + j1] = map1.P21[mi][dj - v_SNOVA][dk][i1 * l_SNOVA + j1];

            for (int dk = v_SNOVA; dk < n_SNOVA; dk++)
                for (int i1 = 0; i1 < rank; i1++)
                    for (int j1 = 0; j1 < rank; j1++)
                        public8[i1 * nl_SNOVA + dk * l_SNOVA + j1] = P22[mi][dj - v_SNOVA][dk - v_SNOVA][i1 * l_SNOVA + j1];

            for (int alpha_i1 = 0; alpha_i1 < (lsq_SNOVA * l_SNOVA + 1) / 2; alpha_i1++)
                for (int dk_j1_32 = 0; dk_j1_32 < nl_SNOVA32; dk_j1_32++)
                    for (int k1 = 0; k1 < l_SNOVA; k1++)
                    {
                        __m256i k_lh0 = mtk2_16[left8[alpha_i1 * n_SNOVA * l_SNOVA + dj * l_SNOVA + k1]];
                        __m256i pubval256 = public256[k1 * nl_SNOVA32 + dk_j1_32];
                        inter256[(alpha_i1 * m_SNOVA + mi) * nl_SNOVA32 + dk_j1_32] ^= _mm256_shuffle_epi8(k_lh0, pubval256);
                    }
        }
    }

    for (int dj = 0; dj < n_SNOVA; dj++) {
        // Shuffle intermediate
        for (int alpha_i1 = 0; alpha_i1 < (lsq_SNOVA * l_SNOVA + 1) / 2; alpha_i1++)
        {
            int i1_a = (2 * alpha_i1) % l_SNOVA;
            int alpha_a = (2 * alpha_i1) / l_SNOVA;
            int i1_b = (2 * alpha_i1 + 1) % l_SNOVA;
            int alpha_b = (2 * alpha_i1 + 1) / l_SNOVA;

            for (int j1 = 0; j1 < l_SNOVA; j1++)
                for (int mi = 0; mi < m_SNOVA; ++mi)
                    {
                        uint8_t ival = inter8[alpha_i1 * mnl_SNOVA + mi * nl_SNOVA + dj * l_SNOVA + j1];

                        r_inter8[alpha_a * l_SNOVA * ml_SNOVA + j1 * ml_SNOVA + i1_a * m_SNOVA + mi] = ival & 0xf;
                        r_inter8[alpha_b * l_SNOVA * ml_SNOVA + j1 * ml_SNOVA + i1_b * m_SNOVA + mi] = (ival >> 4) & 0xf;
                    }
        }

        // Second multiplication loop
        for (int alpha = 0; alpha < lsq_SNOVA; alpha++)
                for (int j1 = 0; j1 < l_SNOVA; j1++)
                    for (int k1 = 0; k1 < l_SNOVA; ++k1)
                        for (int mi_i1 = 0; mi_i1 < ml_SNOVA32; mi_i1++)
                        {
                            __m256i k_lh0 = mtk2_16[right8[alpha * l_SNOVA * nl_SNOVA + k1 * nl_SNOVA + dj * l_SNOVA + j1] & 0x0f];
                            __m256i ival256 = r_inter256[alpha * l_SNOVA * ml_SNOVA32 + j1 * ml_SNOVA32 + mi_i1];
                            res256[k1 * ml_SNOVA32 + mi_i1] ^= _mm256_shuffle_epi8(k_lh0, ival256);
                        }
    }

    // Finish up
    for (int mi = 0; mi < m_SNOVA; ++mi)
        for (int i1 = 0; i1 < l_SNOVA; i1++)
            for (int k1 = 0; k1 < l_SNOVA; ++k1)
                finalres[mi * lsq_SNOVA + i1 * l_SNOVA + k1] = res8[k1 * ml_SNOVA + i1 * m_SNOVA + mi];

    convert_GF16s_to_bytes(hash_in_bytes, finalres, m_SNOVA * lsq_SNOVA);

    for (int i = 0; i < bytes_hash; ++i)
        if (hash_in_bytes[i] != signed_hash[i])
            return -1;

    return 0;
}

#endif
