/**
 * Optimized version in plain C. The compiler will use AVX instructions if available.
 */

#ifndef SNOVA_OPT_H
#define SNOVA_OPT_H

/**
 * Generate private key (F part)
 */
void gen_F_opt(map_group2 *map2, map_group1 *map1, T12_t T12)
{
    uint32_t xF11a[m_SNOVA * v_SNOVA * v_SNOVA * lsq_SNOVA] = {0};
    uint32_t xT12a[v_SNOVA * o_SNOVA * lsq_SNOVA] = {0};
    uint32_t xF11b[m_SNOVA * v_SNOVA * v_SNOVA * lsq_SNOVA] = {0};
    uint32_t xT12b[v_SNOVA * o_SNOVA * lsq_SNOVA] = {0};

    uint32_t xtemp1[m_SNOVA * v_SNOVA * o_SNOVA * l_SNOVA * l_SNOVA] = {0};
    uint32_t xtemp2[m_SNOVA * v_SNOVA * o_SNOVA * l_SNOVA * l_SNOVA] = {0};

    // Prepare

    for (int dj = 0; dj < v_SNOVA; ++dj)
        for (int dk = 0; dk < o_SNOVA; ++dk)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                {
                    uint32_t val = gf16_from_nibble(T12[dj][dk][i1 * l_SNOVA + j1]);
                    xT12a[((dj * l_SNOVA + i1) * o_SNOVA + dk) * l_SNOVA + j1] = val;
                    xT12b[((dj * l_SNOVA + j1) * o_SNOVA + dk) * l_SNOVA + i1] = val;
                }

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int di = 0; di < v_SNOVA; di++)
            for (int dk = 0; dk < v_SNOVA; dk++)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    {
                        uint32_t val = gf16_from_nibble(map1->P11[mi][di][dk][i1 * l_SNOVA + j1]);
                        xF11a[(dk * l_SNOVA + j1) * m_SNOVA * v_SNOVA * l_SNOVA +
                              mi * v_SNOVA * l_SNOVA + di * l_SNOVA + i1] = val;
                        xF11b[(di * l_SNOVA + i1) * m_SNOVA * v_SNOVA * l_SNOVA +
                              mi * v_SNOVA * l_SNOVA + dk * l_SNOVA + j1] = val;
                    }

    // Actual calculations

    memcpy(map2->F11, map1->P11, m_SNOVA * v_SNOVA * v_SNOVA * lsq_SNOVA);
    memcpy(map2->F12, map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * lsq_SNOVA);
    memcpy(map2->F21, map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * lsq_SNOVA);

    for (int dj_j1 = 0; dj_j1 < o_SNOVA * l_SNOVA; ++dj_j1)
        for (int mi_di_i1 = 0; mi_di_i1 < m_SNOVA * v_SNOVA * l_SNOVA; ++mi_di_i1)
            for (int dk_k1 = 0; dk_k1 < v_SNOVA * l_SNOVA; ++dk_k1)
                xtemp1[dj_j1 * m_SNOVA * v_SNOVA * l_SNOVA + mi_di_i1] ^=
                    xF11a[dk_k1 * m_SNOVA * v_SNOVA * l_SNOVA + mi_di_i1] *
                    xT12a[dk_k1 * o_SNOVA * l_SNOVA + dj_j1];

    for (int dj_j1 = 0; dj_j1 < o_SNOVA * l_SNOVA; ++dj_j1)
        for (int mi_di_i1 = 0; mi_di_i1 < m_SNOVA * v_SNOVA * l_SNOVA; ++mi_di_i1)
            for (int dk_k1 = 0; dk_k1 < v_SNOVA * l_SNOVA; ++dk_k1)
                xtemp2[dj_j1 * m_SNOVA * v_SNOVA * l_SNOVA + mi_di_i1] ^=
                    xF11b[dk_k1 * m_SNOVA * v_SNOVA * l_SNOVA + mi_di_i1] *
                    xT12b[dk_k1 * o_SNOVA * l_SNOVA + dj_j1];

    // Convert back

    for (int mi = 0; mi < m_SNOVA; ++mi)
        for (int di = 0; di < v_SNOVA; ++di)
            for (int dj = 0; dj < o_SNOVA; ++dj)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        map2->F12[mi][di][dj][i1 * l_SNOVA + j1] ^=
                            gf16_to_nibble(xtemp1[(dj * l_SNOVA + j1) * m_SNOVA * v_SNOVA * l_SNOVA +
                                                  mi * v_SNOVA * l_SNOVA + di * l_SNOVA + i1]);

    for (int mi = 0; mi < m_SNOVA; ++mi)
        for (int dj = 0; dj < o_SNOVA; ++dj)
            for (int di = 0; di < v_SNOVA; ++di)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    {
                        map2->F21[mi][dj][di][i1 * l_SNOVA + j1] ^=
                            gf16_to_nibble(xtemp2[(dj * l_SNOVA + i1) * m_SNOVA * v_SNOVA * l_SNOVA +
                                                  mi * v_SNOVA * l_SNOVA + di * l_SNOVA + j1]);
                    }
}

/**
 * Computes signature
 */
int sign_digest_core_opt(uint8_t *pt_signature, const uint8_t *digest,
                         uint64_t bytes_digest, uint8_t *array_salt,
                         Aalpha_t Aalpha, Balpha_t Balpha, Qalpha1_t Qalpha1,
                         Qalpha2_t Qalpha2, T12_t T12, F11_t F11, F12_t F12,
                         F21_t F21, const uint8_t pt_public_key_seed[seed_length_public],
                         const uint8_t pt_private_key_seed[seed_length_private])
{
    uint8_t vinegar_gf16[n_SNOVA][lsq_SNOVA] = {0};
    uint32_t xVinegar_gf16[n_SNOVA][lsq_SNOVA] = {0};
    uint32_t temp_xgf16 = 0;

    uint32_t xSolution[m_SNOVA * lsq_SNOVA] = {0};

    uint8_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
    uint8_t signature_in_GF16Matrix[n_SNOVA][lsq_SNOVA];
    uint8_t signed_hash[bytes_hash];

    int flag_redo = 1;
    uint8_t num_sign = 0;

    memset(pt_signature, 0, (bytes_signature + bytes_salt));

    createSignedHash(digest, bytes_digest, pt_public_key_seed, array_salt, signed_hash);
    convert_bytes_to_GF16s(signed_hash, hash_in_GF16, GF16s_hash);

    // Prepare

    uint32_t xT12[v_SNOVA][o_SNOVA][lsq_SNOVA] = {0};
    uint32_t xGauss[m_SNOVA * lsq_SNOVA][m_SNOVA * lsq_SNOVA + 1] = {0};

    uint32_t xF11[m_SNOVA * v_SNOVA * v_SNOVA * l_SNOVA * l_SNOVA] = {0};
    uint32_t xF12[m_SNOVA * v_SNOVA * o_SNOVA * l_SNOVA * l_SNOVA] = {0};
    uint32_t xF21[m_SNOVA * o_SNOVA * v_SNOVA * l_SNOVA * l_SNOVA] = {0};

    uint32_t xAalpha[lsq_SNOVA * lsq_SNOVA] = {0};
    uint32_t xBalpha[lsq_SNOVA * lsq_SNOVA] = {0};
    uint32_t xQalpha1[lsq_SNOVA * lsq_SNOVA] = {0};
    uint32_t xQalpha2[lsq_SNOVA * lsq_SNOVA] = {0};

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int jdx = 0; jdx < v_SNOVA; jdx++)
            for (int kdx = 0; kdx < v_SNOVA; kdx++)
                for (int k1 = 0; k1 < l_SNOVA; ++k1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        xF11[jdx * v_SNOVA * m_SNOVA * l_SNOVA * l_SNOVA + k1 * m_SNOVA * v_SNOVA * l_SNOVA + mi * v_SNOVA * l_SNOVA + kdx * l_SNOVA + j1] =
                            gf16_from_nibble(F11[mi][jdx][kdx][k1 * l_SNOVA + j1]);

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int jdx = 0; jdx < v_SNOVA; jdx++)
            for (int kdx = 0; kdx < o_SNOVA; kdx++)
                for (int k1 = 0; k1 < l_SNOVA; ++k1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        xF12[jdx * l_SNOVA * m_SNOVA * o_SNOVA * l_SNOVA + k1 * m_SNOVA * o_SNOVA * l_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + j1] =
                            gf16_from_nibble(F12[mi][jdx][kdx][k1 * l_SNOVA + j1]);

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int jdx = 0; jdx < v_SNOVA; jdx++)
            for (int kdx = 0; kdx < o_SNOVA; kdx++)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int k1 = 0; k1 < l_SNOVA; ++k1)
                        xF21[(jdx * l_SNOVA + k1) * m_SNOVA * o_SNOVA * l_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + i1] =
                            gf16_from_nibble(F21[mi][kdx][jdx][i1 * l_SNOVA + k1]);

    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int idx = 0; idx < lsq_SNOVA; ++idx)
        {
            xQalpha1[alpha * lsq_SNOVA + idx] = gf16_from_nibble(Qalpha1[alpha][idx]);
            xQalpha2[alpha * lsq_SNOVA + idx] = gf16_from_nibble(Qalpha2[alpha][idx]);
            xAalpha[alpha * lsq_SNOVA + idx] = gf16_from_nibble(Aalpha[alpha][idx]);
            xBalpha[alpha * lsq_SNOVA + idx] = gf16_from_nibble(Balpha[alpha][idx]);
        }

    for (int dj = 0; dj < v_SNOVA; ++dj)
        for (int dk = 0; dk < o_SNOVA; ++dk)
            for (int idx = 0; idx < lsq_SNOVA; ++idx)
                xT12[dj][dk][idx] = gf16_from_nibble(T12[dj][dk][idx]);

    // Try to find a solution

    do
    {
        uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1] = {0};

        uint32_t xLeft[lsq_SNOVA * v_SNOVA * l_SNOVA * l_SNOVA] = {0};
        uint32_t xRight[lsq_SNOVA * v_SNOVA * l_SNOVA * l_SNOVA] = {0};
        uint32_t xFvv_in_GF16Matrix[m_SNOVA][l_SNOVA][l_SNOVA] = {0};

        uint32_t xtemp3[m_SNOVA * lsq_SNOVA * v_SNOVA * l_SNOVA * l_SNOVA] = {0};
        uint32_t xtemp_l[m_SNOVA * o_SNOVA * lsq_SNOVA * l_SNOVA * l_SNOVA] = {0};
        uint32_t xtemp_r[m_SNOVA * o_SNOVA * lsq_SNOVA * l_SNOVA * l_SNOVA] = {0};
        uint32_t xTemp[m_SNOVA][o_SNOVA][l_SNOVA][l_SNOVA][l_SNOVA][l_SNOVA] = {0};
        uint32_t xTemp_left[m_SNOVA][o_SNOVA][lsq_SNOVA][lsq_SNOVA] = {0};
        uint32_t xTemp_right[m_SNOVA][o_SNOVA][lsq_SNOVA][lsq_SNOVA] = {0};

        num_sign++;
        flag_redo = 0;

        // generate the vinegar value
        Keccak_HashInstance hashInstance;
        Keccak_HashInitialize_SHAKE256(&hashInstance);
        Keccak_HashUpdate(&hashInstance, pt_private_key_seed, 8 * seed_length_private);
        Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
        Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
        Keccak_HashUpdate(&hashInstance, &num_sign, 8);
        Keccak_HashFinal(&hashInstance, NULL);
        Keccak_HashSqueeze(&hashInstance, vinegar_in_byte, 8 * ((v_SNOVA * lsq_SNOVA + 1) >> 1));

        convert_bytes_to_GF16s(vinegar_in_byte, (uint8_t *)vinegar_gf16, v_SNOVA * lsq_SNOVA);
        for (int jdx = 0; jdx < v_SNOVA; ++jdx)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    xVinegar_gf16[jdx][i1 * l_SNOVA + j1] = gf16_from_nibble(vinegar_gf16[jdx][i1 * l_SNOVA + j1]);

        // evaluate the vinegar part of central map
        // 4 * V * L^5

        uint32_t xTemp_Q1[lsq_SNOVA][v_SNOVA][lsq_SNOVA] = {0};
        uint32_t xTemp_Q2[lsq_SNOVA][v_SNOVA][lsq_SNOVA] = {0};

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int jdx = 0; jdx < v_SNOVA; ++jdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        for (int k1 = 0; k1 < l_SNOVA; ++k1)
                            xTemp_Q1[alpha][jdx][i1 * l_SNOVA + j1] ^=
                                xVinegar_gf16[jdx][k1 * l_SNOVA + i1] * xQalpha1[alpha * lsq_SNOVA + k1 * l_SNOVA + j1];

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int jdx = 0; jdx < v_SNOVA; ++jdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        xTemp_Q1[alpha][jdx][i1 * l_SNOVA + j1] &= 0x49249249;

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int jdx = 0; jdx < v_SNOVA; ++jdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        for (int k1 = 0; k1 < l_SNOVA; ++k1)
                            xLeft[alpha * l_SNOVA * v_SNOVA * l_SNOVA + i1 * v_SNOVA * l_SNOVA + jdx * l_SNOVA + j1] ^=
                                xAalpha[alpha * lsq_SNOVA + i1 * l_SNOVA + k1] * xTemp_Q1[alpha][jdx][k1 * l_SNOVA + j1];

        for (int idx = 0; idx < v_SNOVA * lsq_SNOVA * lsq_SNOVA; ++idx)
            xLeft[idx] = gf16_reduce(xLeft[idx]);

        // Same for right

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int jdx = 0; jdx < v_SNOVA; ++jdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        for (int k1 = 0; k1 < l_SNOVA; ++k1)
                            xTemp_Q2[alpha][jdx][i1 * l_SNOVA + j1] ^=
                                xQalpha2[alpha * lsq_SNOVA + i1 * l_SNOVA + k1] * xVinegar_gf16[jdx][k1 * l_SNOVA + j1];

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int jdx = 0; jdx < v_SNOVA; ++jdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        xTemp_Q2[alpha][jdx][i1 * l_SNOVA + j1] &= 0x49249249;

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int jdx = 0; jdx < v_SNOVA; ++jdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        for (int k1 = 0; k1 < l_SNOVA; ++k1)
                            xRight[alpha * l_SNOVA * v_SNOVA * l_SNOVA + j1 * v_SNOVA * l_SNOVA + jdx * l_SNOVA + i1] ^=
                                xTemp_Q2[alpha][jdx][i1 * l_SNOVA + k1] * xBalpha[alpha * lsq_SNOVA + k1 * l_SNOVA + j1];

        for (int idx = 0; idx < v_SNOVA * lsq_SNOVA * lsq_SNOVA; ++idx)
            xRight[idx] = gf16_reduce(xRight[idx]);

        // Main multiplication
        // V^2 * O * L^5
        for (int mi_kdx_j1 = 0; mi_kdx_j1 < m_SNOVA * v_SNOVA * l_SNOVA; ++mi_kdx_j1)
            for (int alpha_i1 = 0; alpha_i1 < lsq_SNOVA * l_SNOVA; ++alpha_i1)
                for (int jdx_k1 = 0; jdx_k1 < v_SNOVA * l_SNOVA; ++jdx_k1)
                    xtemp3[alpha_i1 * m_SNOVA * v_SNOVA * l_SNOVA + mi_kdx_j1] ^=
                        xLeft[alpha_i1 * v_SNOVA * l_SNOVA + jdx_k1] * xF11[jdx_k1 * v_SNOVA * m_SNOVA * l_SNOVA + mi_kdx_j1];

        for (int idx = 0; idx < lsq_SNOVA * l_SNOVA * m_SNOVA * v_SNOVA * l_SNOVA; ++idx)
            xtemp3[idx] &= 0x49249249;

        // V * O * L^5
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                        for (int jdx_k1 = 0; jdx_k1 < v_SNOVA * l_SNOVA; ++jdx_k1)
                            xFvv_in_GF16Matrix[mi][i1][j1] ^=
                                xtemp3[alpha * m_SNOVA * v_SNOVA * l_SNOVA * l_SNOVA + i1 * m_SNOVA * v_SNOVA * l_SNOVA + mi * v_SNOVA * l_SNOVA + jdx_k1] *
                                xRight[(alpha * l_SNOVA + j1) * v_SNOVA * l_SNOVA + jdx_k1];

        // compute the coefficients of Xo and put into Gauss matrix and compute
        // the coefficients of Xo^t and add into Gauss matrix
        //
        // 2 * V * O^2 * L^5
        for (int alpha_i1 = 0; alpha_i1 < lsq_SNOVA * l_SNOVA; ++alpha_i1)
            for (int mi_kdx_j1 = 0; mi_kdx_j1 < m_SNOVA * o_SNOVA * l_SNOVA; ++mi_kdx_j1)
                for (int jdk_k1 = 0; jdk_k1 < v_SNOVA * l_SNOVA; ++jdk_k1)
                    xtemp_l[alpha_i1 * m_SNOVA * o_SNOVA * l_SNOVA + mi_kdx_j1] ^=
                        xLeft[alpha_i1 * v_SNOVA * l_SNOVA + jdk_k1] * xF12[jdk_k1 * m_SNOVA * o_SNOVA * l_SNOVA + mi_kdx_j1];

        for (int idx = 0; idx < lsq_SNOVA * l_SNOVA * m_SNOVA * o_SNOVA * l_SNOVA; ++idx)
            xtemp_l[idx] = gf16_reduce(xtemp_l[idx]);

        for (int alpha_j1 = 0; alpha_j1 < lsq_SNOVA * l_SNOVA; ++alpha_j1)
            for (int mi_kdx_i1 = 0; mi_kdx_i1 < m_SNOVA * o_SNOVA * l_SNOVA; ++mi_kdx_i1)
                for (int jdk_k1 = 0; jdk_k1 < v_SNOVA * l_SNOVA; ++jdk_k1)
                    xtemp_r[alpha_j1 * m_SNOVA * o_SNOVA * l_SNOVA + mi_kdx_i1] ^=
                        xRight[alpha_j1 * v_SNOVA * l_SNOVA + jdk_k1] * xF21[jdk_k1 * m_SNOVA * o_SNOVA * l_SNOVA + mi_kdx_i1];

        for (int idx = 0; idx < lsq_SNOVA * l_SNOVA * m_SNOVA * o_SNOVA * l_SNOVA; ++idx)
            xtemp_r[idx] = gf16_reduce(xtemp_r[idx]);

        // Calculate Temp -> Gauss matrix
        // O^2 * L^5
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            for (int k1 = 0; k1 < l_SNOVA; ++k1)
                                xTemp_left[mi][kdx][alpha][i1 * l_SNOVA + j1] ^=
                                    xtemp_l[alpha * l_SNOVA * m_SNOVA * o_SNOVA * l_SNOVA + i1 * m_SNOVA * o_SNOVA * l_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + k1] *
                                    xQalpha2[alpha * lsq_SNOVA + k1 * l_SNOVA + j1];

        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            xTemp_left[mi][kdx][alpha][i1 * l_SNOVA + j1] &= 0x49249249;

        // Outer product
        // O^2 * L^6
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            for (int i2 = 0; i2 < l_SNOVA; ++i2)
                                for (int j2 = 0; j2 < l_SNOVA; ++j2)
                                    xTemp[mi][kdx][i1][j2][i2][j1] ^=
                                        xTemp_left[mi][kdx][alpha][i1 * l_SNOVA + i2] * xBalpha[alpha * lsq_SNOVA + j2 * l_SNOVA + j1];

        // Same for Right
        // O^2 * L^5
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            for (int k1 = 0; k1 < l_SNOVA; ++k1)
                                xTemp_right[mi][kdx][alpha][i1 * l_SNOVA + j1] ^=
                                    xtemp_r[alpha * l_SNOVA * m_SNOVA * o_SNOVA * l_SNOVA + j1 * m_SNOVA * o_SNOVA * l_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + k1] *
                                    xQalpha1[alpha * lsq_SNOVA + i1 * l_SNOVA + k1];

        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            xTemp_right[mi][kdx][alpha][i1 * l_SNOVA + j1] &= 0x49249249;

        // Outer product
        // O^2 * L^6
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            for (int i2 = 0; i2 < l_SNOVA; ++i2)
                                for (int j2 = 0; j2 < l_SNOVA; ++j2)
                                    xTemp[mi][kdx][i1][j2][i2][j1] ^=
                                        xAalpha[alpha * lsq_SNOVA + i1 * l_SNOVA + j2] * xTemp_right[mi][kdx][alpha][i2 * l_SNOVA + j1];

        // Compose Gauss matrix
        // put hash value in the last column of Gauss matrix
        for (int index = 0; index < (m_SNOVA * lsq_SNOVA); index++)
            xGauss[index][m_SNOVA * lsq_SNOVA] = gf16_from_nibble(hash_in_GF16[index]);

        // Reorder xTemp
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        for (int i2 = 0; i2 < l_SNOVA; ++i2)
                            for (int j2 = 0; j2 < l_SNOVA; ++j2)
                                xGauss[mi * lsq_SNOVA + i1 * l_SNOVA + j1][kdx * lsq_SNOVA + i2 * l_SNOVA + j2] = gf16_reduce(xTemp[mi][kdx][i1][j2][i2][j1]);

        // last column of Gauss matrix
        for (int mi = 0; mi < m_SNOVA; mi++)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    xGauss[mi * lsq_SNOVA + i1 * l_SNOVA + j1][m_SNOVA * lsq_SNOVA] ^= gf16_reduce(xFvv_in_GF16Matrix[mi][i1][j1]);

        // Gauss elimination in constant time
        for (int mi2 = 0; mi2 < m_SNOVA * lsq_SNOVA; ++mi2)
        {
            // Find index to swap in constant time
            int swapidx = -1;
            for (int j2 = mi2; j2 < m_SNOVA * lsq_SNOVA; ++j2)
                swapidx += ((swapidx >> 31) & ct_xgf16_is_not_zero(xGauss[j2][mi2])) * (1 + j2);

            flag_redo |= swapidx >> 31;

            // Always swap
            swapidx += ((swapidx >> 31) & 1) * l_SNOVA;
            for (int k2 = mi2; k2 < m_SNOVA * lsq_SNOVA + 1; ++k2)
            {
                temp_xgf16 = xGauss[mi2][k2];
                xGauss[mi2][k2] = xGauss[swapidx][k2];
                xGauss[swapidx][k2] = temp_xgf16;
            }

            temp_xgf16 = gf16_inv(xGauss[mi2][mi2]);
            for (int k2 = mi2; k2 < m_SNOVA * lsq_SNOVA + 1; ++k2)
                xGauss[mi2][k2] = gf16_reduce(xGauss[mi2][k2] * temp_xgf16);

            for (int j2 = mi2 + 1; j2 < m_SNOVA * lsq_SNOVA; ++j2)
            {
                // Constant time version
                temp_xgf16 = ct_xgf16_is_not_zero(xGauss[j2][mi2]) * xGauss[j2][mi2];
                for (int k2 = mi2; k2 < m_SNOVA * lsq_SNOVA + 1; ++k2)
                    xGauss[j2][k2] = gf16_reduce(xGauss[j2][k2] ^ (xGauss[mi2][k2] * temp_xgf16));
            }
        }
    } while (flag_redo);

    for (int mi2 = m_SNOVA * lsq_SNOVA - 1; mi2 >= 0; --mi2)
    {
        temp_xgf16 = 0;
        for (int k2 = mi2 + 1; k2 < m_SNOVA * lsq_SNOVA; ++k2)
            temp_xgf16 ^= xGauss[mi2][k2] * xSolution[k2];

        xSolution[mi2] = xGauss[mi2][m_SNOVA * lsq_SNOVA] ^ gf16_reduce(temp_xgf16);
    }

    for (int index = 0; index < o_SNOVA; ++index)
        for (int i1 = 0; i1 < l_SNOVA; ++i1)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                vinegar_gf16[index + v_SNOVA][i1 * l_SNOVA + j1] = gf16_to_nibble(xSolution[index * lsq_SNOVA + i1 * l_SNOVA + j1]);

    // Establish Signature

    for (int dj = 0; dj < v_SNOVA; ++dj)
    {
        uint32_t xSig[lsq_SNOVA] = {0};

        for (int dk = 0; dk < o_SNOVA; ++dk)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int k1 = 0; k1 < l_SNOVA; ++k1)
                        xSig[i1 * l_SNOVA + j1] ^= xT12[dj][dk][i1 * l_SNOVA + k1] * xSolution[dk * lsq_SNOVA + k1 * l_SNOVA + j1];

        for (int idx = 0; idx < lsq_SNOVA; ++idx)
            signature_in_GF16Matrix[dj][idx] = vinegar_gf16[dj][idx] ^ gf16_to_nibble(xSig[idx]);
    }

    for (int index = 0; index < o_SNOVA; ++index)
        for (int idx = 0; idx < lsq_SNOVA; ++idx)
            signature_in_GF16Matrix[v_SNOVA + index][idx] = vinegar_gf16[v_SNOVA + index][idx];

    // output signature
    convert_GF16s_to_bytes(pt_signature, (gf16_t *)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);
    for (int i1 = 0; i1 < bytes_salt; ++i1)
        pt_signature[bytes_signature + i1] = array_salt[i1];

    return 0;
}

/**
 * Verifies signature
 */
int verify_signture_opt(const uint8_t *pt_digest, uint64_t bytes_digest, const uint8_t *pt_signature, const uint8_t *pk)
{
    uint8_t hash_in_bytes[bytes_hash];
    uint8_t signed_hash[bytes_hash];
    const uint8_t *pt_salt = pt_signature + bytes_signature;

    gf16m_t Left[lsq_SNOVA][n_SNOVA], Right[lsq_SNOVA][n_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA];
    gf16m_t P22[m_SNOVA][o_SNOVA][o_SNOVA];

    map_group1 map1;
    gf16m_t temp1, temp2;

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

    convert_bytes_to_GF16s(pt_signature, (gf16_t *)signature_in_GF16Matrix, GF16s_signature);
    // generate PRNG part of public key
    gen_A_B_Q_P(&map1, pk_stru->pt_public_key_seed);
    // read  P22
    input_P22((uint8_t *)P22, pk_stru->P22);

    // evaluate signature GF16Matrix array
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int index = 0; index < n_SNOVA; ++index)
        {
            gf16m_transpose(signature_in_GF16Matrix[index], temp1);
            gf16m_mul(temp1, map1.Qalpha1[alpha], temp2);
            gf16m_mul(map1.Aalpha[alpha], temp2, Left[alpha][index]);
            gf16m_mul(map1.Qalpha2[alpha], signature_in_GF16Matrix[index], temp2);
            gf16m_mul(temp2, map1.Balpha[alpha], Right[alpha][index]);
        }

    // Align Public and Intermediate to 256 bits
    uint8_t Public[m_SNOVA][n_SNOVA][n_SNOVA][lsq_SNOVA];
    uint32_t xPublic[m_SNOVA][n_SNOVA][lsq_SNOVA][((n_SNOVA + 7) / 8) * 8];
    uint32_t Intermediate[m_SNOVA][lsq_SNOVA][lsq_SNOVA][((n_SNOVA + 7) / 8) * 8] = {0};
    uint32_t xLeft[lsq_SNOVA][n_SNOVA][lsq_SNOVA];
    uint32_t xRight[lsq_SNOVA][lsq_SNOVA][n_SNOVA];
    uint32_t res[m_SNOVA][lsq_SNOVA] = {0};

    // Prepare Left and Right
    for (int idx1 = 0; idx1 < lsq_SNOVA; idx1++)
        for (int idx2 = 0; idx2 < n_SNOVA; idx2++)
            for (size_t idx = 0; idx < lsq_SNOVA; idx++)
            {
                xLeft[idx1][idx2][idx] = gf16_from_nibble(Left[idx1][idx2][idx]);
                xRight[idx1][idx][idx2] = gf16_from_nibble(Right[idx1][idx2][idx]);
            }

    // Prepare Public
    for (int i = 0; i < m_SNOVA; ++i)
    {
        // Convert remaining map1 matrices and combine into one
        for (int idx2 = 0; idx2 < v_SNOVA; idx2++)
            for (int idx3 = 0; idx3 < v_SNOVA; idx3++)
                for (size_t idx = 0; idx < lsq_SNOVA; idx++)
                    Public[i][idx2][idx3][idx] = (map1.P11[i][idx2][idx3][idx]);

        for (int idx2 = 0; idx2 < v_SNOVA; idx2++)
            for (int idx3 = v_SNOVA; idx3 < n_SNOVA; idx3++)
                for (size_t idx = 0; idx < lsq_SNOVA; idx++)
                    Public[i][idx2][idx3][idx] = (map1.P12[i][idx2][idx3 - v_SNOVA][idx]);

        for (int idx2 = v_SNOVA; idx2 < n_SNOVA; idx2++)
            for (int idx3 = 0; idx3 < v_SNOVA; idx3++)
                for (size_t idx = 0; idx < lsq_SNOVA; idx++)
                    Public[i][idx2][idx3][idx] = (map1.P21[i][idx2 - v_SNOVA][idx3][idx]);

        for (int idx2 = v_SNOVA; idx2 < n_SNOVA; idx2++)
            for (int idx3 = v_SNOVA; idx3 < n_SNOVA; idx3++)
                for (size_t idx = 0; idx < lsq_SNOVA; idx++)
                    Public[i][idx2][idx3][idx] = (P22[i][idx2 - v_SNOVA][idx3 - v_SNOVA][idx]);
    }

    for (int i = 0; i < m_SNOVA; ++i)
        for (int idx2 = 0; idx2 < n_SNOVA; idx2++)
            for (size_t idx = 0; idx < lsq_SNOVA; idx++)
                for (int idx3 = 0; idx3 < n_SNOVA; idx3++)
                    xPublic[i][idx2][idx][idx3] = gf16_from_nibble(Public[i][idx2][idx3][idx]);

    // Main loop
    for (int i = 0; i < m_SNOVA; ++i)
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int dj = 0; dj < n_SNOVA; ++dj)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int dk = 0; dk < n_SNOVA; ++dk)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            for (int k1 = 0; k1 < l_SNOVA; ++k1)
                                Intermediate[i][alpha][i1 * l_SNOVA + j1][dk] ^= xLeft[alpha][dj][i1 * l_SNOVA + k1] *
                                                                                 xPublic[i][dj][k1 * l_SNOVA + j1][dk];

    // Reduce for next multiplication
    for (int i = 0; i < m_SNOVA; ++i)
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int dk = 0; dk < n_SNOVA; ++dk)
                        Intermediate[i][alpha][i1 * l_SNOVA + j1][dk] &= 0x49249249;

    // Second loop
    for (int i = 0; i < m_SNOVA; ++i)
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int k1 = 0; k1 < l_SNOVA; ++k1)
                        for (int dk = 0; dk < n_SNOVA; ++dk)
                            res[i][i1 * l_SNOVA + j1] ^= Intermediate[i][alpha][i1 * l_SNOVA + k1][dk] *
                                                         xRight[alpha][k1 * l_SNOVA + j1][dk];

    // Finish up
    for (int mi = 0; mi < m_SNOVA; ++mi)
        for (int i1 = 0; i1 < l_SNOVA; ++i1)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                ((gf16_t *)signature_in_GF16Matrix)[mi * lsq_SNOVA + i1 * l_SNOVA + j1] = gf16_to_nibble(res[mi][i1 * l_SNOVA + j1]);
    convert_GF16s_to_bytes(hash_in_bytes, (gf16_t *)signature_in_GF16Matrix, m_SNOVA * lsq_SNOVA);

    int result = 0;
    for (int i = 0; i < bytes_hash; ++i)
        if (hash_in_bytes[i] != signed_hash[i])
        {
            result = -1;
            break;
        }

    return result;
}

#endif
