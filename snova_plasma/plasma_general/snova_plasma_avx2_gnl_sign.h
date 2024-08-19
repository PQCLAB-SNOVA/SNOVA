/**
 * VTL version for any l_SNOVA.
 */

#ifndef PLASMA_GNL_AVX2_SIGN_H
#define PLASMA_GNL_AVX2_SIGN_H

#include <immintrin.h>
#include <stdint.h>

#define mvl_SNOVA32 ((m_SNOVA * v_SNOVA * l_SNOVA + 31) / 32)
#define mvl_SNOVA (mvl_SNOVA32 * 32)

#define ml_SNOVA32 ((m_SNOVA * l_SNOVA + 31) / 32)
#define ml_SNOVA (ml_SNOVA32 * 32)

#define vl_SNOVA (v_SNOVA * l_SNOVA)
#define lcube_SNOVA (((lsq_SNOVA * l_SNOVA + 31) / 32) * 32)

#define vl4_SNOVA32 ((vl_SNOVA * lcube_SNOVA + 31) / 32)

#define GAUSS_ROW (m_SNOVA * lsq_SNOVA + 1)
#define GAUSS_ROW32 ((GAUSS_ROW + 31) / 32)
#define GAUSS_ROW_mult32 (GAUSS_ROW32 * 32)
#define GAUSS_COL (m_SNOVA * lsq_SNOVA)
#define GAUSS_COL32 ((GAUSS_COL + 31) / 32 * 32)

#ifndef v_SNOVA_mult4
#define v_SNOVA_mult4 ((v_SNOVA + 3) / 4 * 4)
#endif


/**
 * Computes signature
 */
int sign_digest_core_gnl_vtl(uint8_t *pt_signature, const uint8_t *digest,
                         uint64_t bytes_digest, uint8_t *array_salt,
                         Aalpha_t Aalpha, Balpha_t Balpha, Qalpha1_t Qalpha1,
                         Qalpha2_t Qalpha2, T12_t T12, F11_t F11, F12_t F12,
                         F21_t F21, const uint8_t pt_public_key_seed[seed_length_public],
                         const uint8_t pt_private_key_seed[seed_length_private])
{
    uint8_t vinegar_gf16[n_SNOVA][lsq_SNOVA] = {0};
    uint32_t xVinegar_gf16[n_SNOVA][lsq_SNOVA] = {0};
    uint8_t solution[GAUSS_COL32] = {0};
    uint8_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
    uint8_t signed_hash[bytes_hash];

    int flag_redo = 1;
    uint8_t num_sign = 0;

    memset(pt_signature, 0, (bytes_signature + bytes_salt));
    createSignedHash(digest, bytes_digest, pt_public_key_seed, array_salt, signed_hash);
    convert_bytes_to_GF16s(signed_hash, hash_in_GF16, GF16s_hash);

    // Prepare

    uint8_t Gauss[m_SNOVA * lsq_SNOVA][GAUSS_ROW_mult32] __attribute__((aligned(32)));

    __m256i f11_256[mvl_SNOVA32 * vl_SNOVA] = {0};
    __m256i f12_256[mol_SNOVA32 * vl_SNOVA] = {0};
    __m256i f21_256[mol_SNOVA32 * vl_SNOVA] = {0};

    uint8_t *f11_8 = (uint8_t *)f11_256;
    uint8_t *f12_8 = (uint8_t *)f12_256;
    uint8_t *f21_8 = (uint8_t *)f21_256;

    uint32_t xAalpha[lsq_SNOVA * lsq_SNOVA];
    uint32_t xBalpha[lsq_SNOVA * lsq_SNOVA];
    uint32_t xQalpha1[lsq_SNOVA * lsq_SNOVA];
    uint32_t xQalpha2[lsq_SNOVA * lsq_SNOVA];

    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        for (int idx = 0; idx < lsq_SNOVA; ++idx)
        {
            xQalpha1[alpha * lsq_SNOVA + idx] = gf16_from_nibble(Qalpha1[alpha][idx]);
            xQalpha2[alpha * lsq_SNOVA + idx] = gf16_from_nibble(Qalpha2[alpha][idx]);
            xAalpha[alpha * lsq_SNOVA + idx] = gf16_from_nibble(Aalpha[alpha][idx]);
            xBalpha[alpha * lsq_SNOVA + idx] = gf16_from_nibble(Balpha[alpha][idx]);
        }

    for (int jdx = 0; jdx < v_SNOVA; jdx++)
        for (int k1 = 0; k1 < l_SNOVA; ++k1)
            for (int mi = 0; mi < m_SNOVA; mi++)
                for (int kdx = 0; kdx < v_SNOVA; kdx++)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        f11_8[jdx * l_SNOVA * mvl_SNOVA + k1 * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + kdx * l_SNOVA + j1] =
                            F11[mi][jdx][kdx][k1 * l_SNOVA + j1];

    for (int jdx = 0; jdx < v_SNOVA; jdx++)
        for (int k1 = 0; k1 < l_SNOVA; ++k1)
            for (int mi = 0; mi < m_SNOVA; mi++)
                for (int kdx = 0; kdx < o_SNOVA; kdx++)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                        f12_8[jdx * l_SNOVA * mol_SNOVA + k1 * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + j1] =
                            F12[mi][jdx][kdx][k1 * l_SNOVA + j1];

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int jdx = 0; jdx < v_SNOVA; jdx++)
            for (int kdx = 0; kdx < o_SNOVA; kdx++)
                for (int k1 = 0; k1 < l_SNOVA; ++k1)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        f21_8[jdx * l_SNOVA * mol_SNOVA + k1 * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + i1] =
                            F21[mi][kdx][jdx][i1 * l_SNOVA + k1];

    // Try to find a solution
    do
    {
        uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1] = {0};

        uint32_t xLeft[lsq_SNOVA * v_SNOVA * l_SNOVA * l_SNOVA] = {0};
        uint32_t xRight[lsq_SNOVA * v_SNOVA * l_SNOVA * l_SNOVA] = {0};

        uint32_t xTemp_Q1[lsq_SNOVA][v_SNOVA][lsq_SNOVA] = {0};
        uint32_t xTemp_Q2[lsq_SNOVA][v_SNOVA][lsq_SNOVA] = {0};

        // TODO: Reuse some of these variables to save memory usage
        __m256i left256[vl4_SNOVA32] = {0};
        __m256i right256[vl4_SNOVA32] = {0};

        __m256i temp3_256[mvl_SNOVA32 * lcube_SNOVA] = {0};
        __m256i temp_l256[mol_SNOVA32 * lcube_SNOVA] = {0};
        __m256i temp_r256[mol_SNOVA32 * lcube_SNOVA] = {0};

        __m256i temp_left256[mol_SNOVA32 * lcube_SNOVA] = {0};
        __m256i temp_right256[mol_SNOVA32 * lcube_SNOVA] = {0};
        __m256i res_left256[mol_SNOVA32 * lcube_SNOVA] = {0};
        __m256i res_right256[mol_SNOVA32 * lcube_SNOVA] = {0};
        __m256i gtemp256[mol_SNOVA32 * lcube_SNOVA] = {0};

        __m256i fvv_res256[ml_SNOVA32 * l_SNOVA] = {0};
        __m256i fvv_temp256[ml_SNOVA32 * v_SNOVA * lcube_SNOVA] = {0};

        uint8_t *left8 = (uint8_t *)left256;
        uint8_t *right8 = (uint8_t *)right256;
        uint8_t *temp3_8 = (uint8_t *)temp3_256;
        uint8_t *temp_l8 = (uint8_t *)temp_l256;
        uint8_t *temp_r8 = (uint8_t *)temp_r256;
        uint8_t *temp_left8 = (uint8_t *)temp_left256;
        uint8_t *temp_right8 = (uint8_t *)temp_right256;
        uint8_t *gtemp8 = (uint8_t *)gtemp256;
        uint8_t *fvv_temp8 = (uint8_t *)fvv_temp256;
        uint8_t *fvv_res8 = (uint8_t *)fvv_res256;

        num_sign++;
        flag_redo = 0;

        memset(Gauss, 0, sizeof(Gauss));

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

        for (int jdx = 0; jdx < v_SNOVA; ++jdx)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
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
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int k1 = 0; k1 < l_SNOVA; ++k1)
                            xLeft[alpha * vl_SNOVA * l_SNOVA + i1 * v_SNOVA * l_SNOVA + jdx * l_SNOVA + j1] ^=
                                xAalpha[alpha * lsq_SNOVA + i1 * l_SNOVA + k1] * xTemp_Q1[alpha][jdx][k1 * l_SNOVA + j1];

        for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; ++alpha_i1)
            for (int jdx_k1 = 0; jdx_k1 < vl_SNOVA; ++jdx_k1)
                xLeft[alpha_i1 * vl_SNOVA + jdx_k1] = gf16_reduce(xLeft[alpha_i1 * vl_SNOVA + jdx_k1]);

        for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; ++alpha_i1)
            for (int jdx_k1 = 0; jdx_k1 < vl_SNOVA; ++jdx_k1)
                left8[alpha_i1 * vl_SNOVA + jdx_k1] = xgf16_to_nibble(xLeft[alpha_i1 * vl_SNOVA + jdx_k1]);

        // Same for right
        for (int jdx = 0; jdx < v_SNOVA; ++jdx)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
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
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int k1 = 0; k1 < l_SNOVA; ++k1)
                            xRight[alpha * vl_SNOVA * l_SNOVA + j1 * v_SNOVA * l_SNOVA + jdx * l_SNOVA + i1] ^=
                                xTemp_Q2[alpha][jdx][i1 * l_SNOVA + k1] * xBalpha[alpha * lsq_SNOVA + k1 * l_SNOVA + j1];

        for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; ++alpha_i1)
            for (int jdx_k1 = 0; jdx_k1 < vl_SNOVA; ++jdx_k1)
                xRight[alpha_i1 * vl_SNOVA + jdx_k1] = gf16_reduce(xRight[alpha_i1 * vl_SNOVA + jdx_k1]);

        for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; ++alpha_i1)
            for (int jdx_k1 = 0; jdx_k1 < vl_SNOVA; ++jdx_k1)
                right8[alpha_i1 * vl_SNOVA + jdx_k1] = xgf16_to_nibble(xRight[alpha_i1 * vl_SNOVA + jdx_k1]);

        // Prepare main

        for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; alpha_i1 += 2)
            for (int jdx_k1 = 0; jdx_k1 < vl_SNOVA; jdx_k1++)
            {
                __m256i k_lh = vtl_multtab(left8[alpha_i1 * vl_SNOVA + jdx_k1]) ^
                               _mm256_slli_epi16(vtl_multtab(left8[alpha_i1 * vl_SNOVA + jdx_k1 + vl_SNOVA]), 4);

                for (int mi_kdx_j1 = 0; mi_kdx_j1 < mvl_SNOVA32; mi_kdx_j1++)
                    temp3_256[alpha_i1 * mvl_SNOVA32 + mi_kdx_j1] ^= _mm256_shuffle_epi8(k_lh, f11_256[jdx_k1 * mvl_SNOVA32 + mi_kdx_j1]);
            }

        // Separate twin result
        __m256i l_mask = _mm256_set1_epi64x(0x0f0f0f0f0f0f0f0full);
        for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; alpha_i1 += 2)
            for (int mi_kdx_j1 = 0; mi_kdx_j1 < mvl_SNOVA32; mi_kdx_j1++)
            {
                __m256i res = temp3_256[alpha_i1 * mvl_SNOVA32 + mi_kdx_j1];
                temp3_256[alpha_i1 * mvl_SNOVA32 + mi_kdx_j1] = res & l_mask;
                temp3_256[(alpha_i1 + 1) * mvl_SNOVA32 + mi_kdx_j1] = _mm256_srli_epi16(res, 4) & l_mask;
            }

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int kdx = 0; kdx < v_SNOVA; kdx++)
                for (int mi = 0; mi < m_SNOVA; mi++)
                    for (int i1 = 0; i1 < l_SNOVA; i1++)
                        for (int j1 = 0; j1 < l_SNOVA; j1++)
                        {
                            fvv_temp8[alpha * l_SNOVA * v_SNOVA * ml_SNOVA + j1 * v_SNOVA * ml_SNOVA + kdx * ml_SNOVA + mi * l_SNOVA + i1] =
                                temp3_8[alpha * l_SNOVA * mvl_SNOVA + i1 * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + kdx * l_SNOVA + j1];
                        }

        // V * O * L^5
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int jdx = 0; jdx < v_SNOVA; ++jdx)
                for (int k1 = 0; k1 < l_SNOVA; ++k1)
                    for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    {
                        __m256i k_lh = vtl_multtab(right8[(alpha * l_SNOVA + j1) * vl_SNOVA + jdx * l_SNOVA + k1]);
                        for (int mi_i1 = 0; mi_i1 < ml_SNOVA32; mi_i1++)
                            fvv_res256[j1 * ml_SNOVA32 + mi_i1] ^=
                                _mm256_shuffle_epi8(k_lh, fvv_temp256[((alpha * l_SNOVA + k1) * v_SNOVA + jdx) * ml_SNOVA32 + mi_i1]);
                    }

        // Compose Gauss matrix
        // last column of Gauss matrix
        for (int mi = 0; mi < m_SNOVA; mi++)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    Gauss[mi * lsq_SNOVA + i1 * l_SNOVA + j1][m_SNOVA * lsq_SNOVA] =
                        hash_in_GF16[mi * lsq_SNOVA + i1 * l_SNOVA + j1] ^ fvv_res8[j1 * ml_SNOVA + mi * l_SNOVA + i1];

        // compute the coefficients of Xo and put into Gauss matrix and compute
        // the coefficients of Xo^t and add into Gauss matrix
        //
        // 2 * V * O^2 * L^5

        for (int jdx_k1 = 0; jdx_k1 < vl_SNOVA; jdx_k1++)
            for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; alpha_i1++)
                {
                    __m256i k_lh = vtl_multtab(left8[alpha_i1 * vl_SNOVA + jdx_k1]);
                    for (int mi_kdx_j1 = 0; mi_kdx_j1 < mol_SNOVA32; mi_kdx_j1++)
                        temp_l256[alpha_i1 * mol_SNOVA32 + mi_kdx_j1] ^= _mm256_shuffle_epi8(k_lh, f12_256[jdx_k1 * mol_SNOVA32 + mi_kdx_j1]);
                }

        for (int jdx_k1 = 0; jdx_k1 < vl_SNOVA; jdx_k1++)
            for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; alpha_i1++)
                {
                    __m256i k_lh = vtl_multtab(right8[alpha_i1 * vl_SNOVA + jdx_k1]);
                    for (int mi_kdx_j1 = 0; mi_kdx_j1 < mol_SNOVA32; mi_kdx_j1++)
                        temp_r256[alpha_i1 * mol_SNOVA32 + mi_kdx_j1] ^= _mm256_shuffle_epi8(k_lh, f21_256[jdx_k1 * mol_SNOVA32 + mi_kdx_j1]);
                }

        // Shuffle Left
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int mi = 0; mi < m_SNOVA; ++mi)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            temp_left8[alpha * l_SNOVA * mol_SNOVA + j1 * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + i1] =
                                temp_l8[(alpha * l_SNOVA + i1) * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + j1];

        // O^2 * L^5
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int k1 = 0; k1 < l_SNOVA; ++k1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                {
                    __m256i k_lh = mtk2_16[Qalpha2[alpha][k1 * l_SNOVA + j1]];
                    for (int mi_kdx_i1 = 0; mi_kdx_i1 < mol_SNOVA32; mi_kdx_i1++)
                        res_left256[(alpha * l_SNOVA + j1) * mol_SNOVA32 + mi_kdx_i1] ^=
                            _mm256_shuffle_epi8(k_lh, temp_left256[(alpha * l_SNOVA + k1) * mol_SNOVA32 + mi_kdx_i1]);
                }

        // Outer product, O^2 * L^6
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                for (int j2 = 0; j2 < l_SNOVA; ++j2)
                {
                    __m256i k_lh = mtk2_16[Balpha[alpha][j2 * l_SNOVA + j1]];
                    for (int i2_mi_kdx_i1 = 0; i2_mi_kdx_i1 < l_SNOVA * mol_SNOVA32; i2_mi_kdx_i1++)
                        gtemp256[(j2 * l_SNOVA + j1) * l_SNOVA * mol_SNOVA32 + i2_mi_kdx_i1] ^=
                            _mm256_shuffle_epi8(k_lh, res_left256[alpha * l_SNOVA * mol_SNOVA32 + i2_mi_kdx_i1]);
                }

        // Shuffle to Gauss matrix
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int i2 = 0; i2 < l_SNOVA; ++i2)
                        for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                            for (int j2 = 0; j2 < l_SNOVA; ++j2)
                                Gauss[mi * lsq_SNOVA + i1 * l_SNOVA + j1][kdx * lsq_SNOVA + i2 * l_SNOVA + j2] ^=
                                    gtemp8[(j2 * l_SNOVA + j1) * l_SNOVA * mol_SNOVA + i2 * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + i1];

        // Same for Right
        memset(gtemp256, 0, sizeof(gtemp256));

        // Shuffle
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int mi = 0; mi < m_SNOVA; ++mi)
                    for (int i1 = 0; i1 < l_SNOVA; ++i1)
                        for (int j1 = 0; j1 < l_SNOVA; ++j1)
                            temp_right8[alpha * l_SNOVA * mol_SNOVA + j1 * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + i1] =
                                temp_r8[(alpha * l_SNOVA + i1) * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + j1];

        // O^2 * L^5
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int k1 = 0; k1 < l_SNOVA; ++k1)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                {
                    __m256i k_lh = mtk2_16[Qalpha1[alpha][i1 * l_SNOVA + k1]];
                    for (int mi_kdx_j1 = 0; mi_kdx_j1 < mol_SNOVA32; mi_kdx_j1++)
                        res_right256[(alpha * l_SNOVA + i1) * mol_SNOVA32 + mi_kdx_j1] ^=
                            _mm256_shuffle_epi8(k_lh, temp_right256[(alpha * l_SNOVA + k1) * mol_SNOVA32 + mi_kdx_j1]);
                }

        // Outer product, O^2 * L^6
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j2 = 0; j2 < l_SNOVA; ++j2)
                {
                    __m256i k_lh = mtk2_16[Aalpha[alpha][i1 * l_SNOVA + j2]];
                    for (int i2_mi_kdx_j1 = 0; i2_mi_kdx_j1 < l_SNOVA * mol_SNOVA32; i2_mi_kdx_j1++)
                        gtemp256[(j2 * l_SNOVA + i1) * l_SNOVA * mol_SNOVA32 + i2_mi_kdx_j1] ^=
                            _mm256_shuffle_epi8(k_lh, res_right256[alpha * l_SNOVA * mol_SNOVA32 + i2_mi_kdx_j1]);
                }

        // Store in Gauss
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int j2 = 0; j2 < l_SNOVA; ++j2)
                        for (int i2 = 0; i2 < l_SNOVA; ++i2)
                            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                                Gauss[mi * lsq_SNOVA + i1 * l_SNOVA + j1][kdx * lsq_SNOVA + i2 * l_SNOVA + j2] ^=
                                    gtemp8[(j2 * l_SNOVA + i1) * l_SNOVA * mol_SNOVA + i2 * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA + j1];

        /**********************************************/

        // Gauss elimination in constant time
        for (int mi2 = 0; mi2 < m_SNOVA * lsq_SNOVA; ++mi2)
        {
            int swap = ct_gf16_is_not_zero(Gauss[mi2][mi2]) - 1;
            for (int j2 = mi2 + 1; j2 < m_SNOVA * lsq_SNOVA; ++j2) {
                __m256i swap256 = _mm256_set1_epi32(swap);
                __m256i *gdest = (__m256i *)&Gauss[mi2][0];
                __m256i *gsource = (__m256i *)&Gauss[j2][0];
                for (int k2 = 0; k2 < GAUSS_ROW32; ++k2)
                    gdest[k2] ^= gsource[k2] & swap256;

                swap = ct_gf16_is_not_zero(Gauss[mi2][mi2]) - 1;
            }
            flag_redo |= swap;

            // Constant time GF16 inverse
            __m256i res256 = _mm256_shuffle_epi8(avx_inv_table, _mm256_set1_epi8(Gauss[mi2][mi2]));
            uint8_t t_GF16 = _mm256_extract_epi8(res256, 0);

            int kstart = (mi2 / 32) * 32;
            for (int k = kstart; k < GAUSS_ROW_mult32; k += 32)
                gf16_32_mul_k(Gauss[mi2] + k, t_GF16, Gauss[mi2] + k);

            for (int j2 = mi2 + 1; j2 < m_SNOVA * lsq_SNOVA; ++j2)
            {
                t_GF16 = Gauss[j2][mi2];
                for (int k2 = kstart; k2 < GAUSS_ROW_mult32; k2 += 32)
                    gf16_32_mul_k_add(Gauss[mi2] + k2, t_GF16, (Gauss[j2] + k2));
            }
        }

        // Cleanup
        if (!flag_redo) {
            SNOVA_CLEAR(vinegar_in_byte)
            SNOVA_CLEAR(xLeft)
            SNOVA_CLEAR(xRight)
            SNOVA_CLEAR(left256)
            SNOVA_CLEAR(right256)
            SNOVA_CLEAR(temp3_256)
            SNOVA_CLEAR(temp_l256)
            SNOVA_CLEAR(temp_r256)
            SNOVA_CLEAR(temp_left256)
            SNOVA_CLEAR(temp_right256)
            SNOVA_CLEAR(res_left256)
            SNOVA_CLEAR(res_right256)
            SNOVA_CLEAR(gtemp256)
            SNOVA_CLEAR(fvv_res256)
            SNOVA_CLEAR(fvv_temp256)
        }
    } while (flag_redo);

    // printf("times of Gauss elimination : %d\n", num_sign);

    for (int mil2 = m_SNOVA * lsq_SNOVA - 1; mil2 >= 0; --mil2)
    {
        uint8_t t_GF16 = 0;
        uint8_t Gauss_last_col = Gauss[mil2][m_SNOVA * lsq_SNOVA];
        __m256i t_GF16_256 = _mm256_setzero_si256();

        Gauss[mil2][m_SNOVA * lsq_SNOVA] = 0;

        int kstart = ((mil2 + 1) / 32) * 32;
        for (int k2 = kstart; k2 < GAUSS_COL32; k2 += 32)
            gf16_32_mul_32_add(Gauss[mil2] + k2, solution + k2, (uint8_t *)&t_GF16_256);

        uint8_t *t_GF16_256_8_ptr = (uint8_t *)(&t_GF16_256);
        for (int k2 = 0; k2 < 32; k2++)
            t_GF16 ^= t_GF16_256_8_ptr[k2];

        solution[mil2] = Gauss_last_col ^ t_GF16;
    }

    // Establish Signature

    for (int index = 0; index < o_SNOVA; ++index)
        for (int i1 = 0; i1 < l_SNOVA; ++i1)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                vinegar_gf16[index + v_SNOVA][i1 * l_SNOVA + j1] = solution[index * lsq_SNOVA + i1 * l_SNOVA + j1];

    gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
    gf16m_t signature_in_GF16Matrix[n_SNOVA] = {0};
    gf16m_t gf16m_secret_temp0;

    memcpy((uint8_t *)X_in_GF16Matrix, (uint8_t *)vinegar_gf16, n_SNOVA * lsq_SNOVA);

    for (int index = 0; index < v_SNOVA; ++index)
    {
        gf16m_clone(signature_in_GF16Matrix[index], X_in_GF16Matrix[index]);
        for (int i = 0; i < o_SNOVA; ++i)
        {
            gf16m_mul(T12[index][i], X_in_GF16Matrix[v_SNOVA + i], gf16m_secret_temp0);
            gf16m_add(signature_in_GF16Matrix[index], gf16m_secret_temp0, signature_in_GF16Matrix[index]);
        }
    }

    for (int index = 0; index < o_SNOVA; ++index)
        gf16m_clone(signature_in_GF16Matrix[v_SNOVA + index], X_in_GF16Matrix[v_SNOVA + index]);

    convert_GF16s_to_bytes(pt_signature, (gf16_t *)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);

    // output signature
    memcpy(pt_signature + bytes_signature, array_salt, bytes_salt);

    // Cleanup
    SNOVA_CLEAR(vinegar_gf16)
    SNOVA_CLEAR(xVinegar_gf16)
    SNOVA_CLEAR(solution)
    SNOVA_CLEAR(hash_in_GF16)
    SNOVA_CLEAR(signed_hash)
    SNOVA_CLEAR(Gauss)
    SNOVA_CLEAR(f11_256)
    SNOVA_CLEAR(f12_256)
    SNOVA_CLEAR(f21_256)
    SNOVA_CLEAR(gf16m_secret_temp0)

    return 0;
}

#endif
