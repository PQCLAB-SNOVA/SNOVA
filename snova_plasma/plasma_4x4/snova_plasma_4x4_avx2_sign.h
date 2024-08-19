/*
Use uint64_t to store a GF16 4x4 matrix (with each half-byte storing one element of GF16)."
Optimizing GF16 4x4 matrix multiplication using AVX2.
*/

#ifndef PLASMA_4x4_AVX2_SIGN_H
#define PLASMA_4x4_AVX2_SIGN_H

#include <immintrin.h>
#include <stdint.h>

#include "../../deriv_params.h"
#include "../../gf16.h"
#include "../../snova_kernel.h"
#include "../snova_plasma_avx2.h"
// uint64_t to store a GF16 4x4 matrix
#include "snova_plasma_4x4_avx2.h"
#include "snova_plasma_4x4_data.h"

/**
 * Computes signature
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param Aalpha -
 * @param Balpha -
 * @param Qalpha1 -
 * @param Qalpha2 -
 * @param T12 -
 * @param F11 -
 * @param F12 -
 * @param F21 -
 * @param pt_public_key_seed - pointer to output public key seed.
 * @param pt_private_key_seed - pointer to output private key seed.
 */
int sign_digest_core_4x4_avx2_vtl(uint8_t *pt_signature, const uint8_t *digest, uint64_t bytes_digest, uint8_t *array_salt, //
                                  Aalpha_t Aalpha, Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2,                   //
                                  T12_t T12, F11_t F11, F12_t F12, F21_t F21,                                               //
                                  const uint8_t pt_public_key_seed[seed_length_public],                                     //
                                  const uint8_t pt_private_key_seed[seed_length_private])
{
    const __m256i l_mask = _mm256_set1_epi64x(0x0f0f0f0f0f0f0f0full);
    const __m128i transpose4x4 = _mm_set_epi64x(0x0f0b07030e0a0602ll, 0x0d0905010c080400ll);

    uint8_t vinegar_gf16[n_SNOVA][lsq_SNOVA] = {0};
    uint8_t solution[GAUSS_ROW_mult32] = {0};
    uint8_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
    uint8_t signed_hash[bytes_hash];

    int flag_redo = 1;
    uint8_t num_sign = 0;

    memset(pt_signature, 0, (bytes_signature + bytes_salt));
    createSignedHash(digest, bytes_digest, pt_public_key_seed, array_salt, signed_hash);
    convert_bytes_to_GF16s(signed_hash, hash_in_GF16, GF16s_hash);

    // Prepare

    uint8_t Gauss[m_SNOVA * lsq_SNOVA][GAUSS_ROW_mult32] __attribute__((aligned(32)));

    uint8_t f11_8[mvl_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t f12_8[mol_SNOVA * vl_SNOVA] __attribute__((aligned(32))) = {0};
    uint8_t f21_8[mol_SNOVA * vl_SNOVA] __attribute__((aligned(32))) = {0};

    __m256i *f11_256 = (__m256i *)f11_8;
    __m256i *f12_256 = (__m256i *)f12_8;
    __m256i *f21_256 = (__m256i *)f21_8;


    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int jdx = 0; jdx < v_SNOVA; jdx++)
            for (int kdx = 0; kdx < o_SNOVA; kdx++)
                for (int k1 = 0; k1 < l_SNOVA; ++k1)
                {
                    uint32_t *target =
                        (uint32_t *)&f12_8[jdx * l_SNOVA * mol_SNOVA + k1 * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA];
                    uint32_t *source = (uint32_t *)&F12[mi][jdx][kdx][k1 * l_SNOVA];
                    *target = *source;
                }

    for (int mi = 0; mi < m_SNOVA; mi++)
        for (int jdx = 0; jdx < v_SNOVA; jdx++)
            for (int kdx = 0; kdx < o_SNOVA; kdx++)
            {
                uint32_t *pp11 = (uint32_t *)F21[mi][kdx][jdx];
                __m128i val128 = _mm_setr_epi32(pp11[0], pp11[1], pp11[2], pp11[3]);
                uint32_t tres[4];
                _mm_storeu_si128((__m128i_u *)&tres, _mm_shuffle_epi8(val128, transpose4x4));
                uint32_t *pres = (uint32_t *)&f21_8[jdx * l_SNOVA * mol_SNOVA + mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA];

                for (int idx = 0; idx < l_SNOVA; idx++)
                    pres[idx * mol_SNOVA / 4] = tres[idx];
            }

    // Try to find a solution
    do
    {
        uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1] = {0};
        uint8_t left_b[vl_SNOVA * lsq_SNOVA * l_SNOVA];
        uint8_t right_b[vl_SNOVA * lsq_SNOVA * l_SNOVA];

        // TODO: Reuse some of these variables to save on memory
        uint8_t temp3_8[mvl_SNOVA * lcube_SNOVA] __attribute__((aligned(32))) = {0};
        uint8_t temp_l8[mol_SNOVA * l_SNOVA] __attribute__((aligned(32)));
        uint8_t temp_r8[mol_SNOVA * l_SNOVA] __attribute__((aligned(32))) = {0};

        uint8_t temp_left8[mol_SNOVA * l_SNOVA] __attribute__((aligned(32))) = {0};
        uint8_t temp_right8[mol_SNOVA * l_SNOVA] __attribute__((aligned(32))) = {0};
        uint8_t gtemp8[mol_SNOVA * lcube_SNOVA] __attribute__((aligned(32))) = {0};
        uint8_t fvv_temp8[ml_SNOVA * v_SNOVA * l_SNOVA] __attribute__((aligned(32))) = {0};
        uint8_t fvv_res8[ml_SNOVA * l_SNOVA] = {0};

        __m256i res_left256[mol_SNOVA32 * l_SNOVA] = {0};
        __m256i res_right256[mol_SNOVA32 * l_SNOVA] = {0};

        __m256i *temp3_256 = (__m256i *)temp3_8;
        __m256i *temp_l256 = (__m256i *)temp_l8;
        __m256i *temp_r256 = (__m256i *)temp_r8;
        __m256i *temp_left256 = (__m256i *)temp_left8;
        __m256i *temp_right256 = (__m256i *)temp_right8;
        __m256i *gtemp256 = (__m256i *)gtemp8;
        __m256i *fvv_res256 = (__m256i *)fvv_res8;
        __m256i *fvv_temp256 = (__m256i *)fvv_temp8;
        uint32_t *fvv_temp32 = (uint32_t *)fvv_temp8;

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

        // u64
        uint64_t Left_u64[v_SNOVA][lsq_SNOVA] __attribute__((aligned(32))) = {0};
        uint64_t Right_u64[v_SNOVA][lsq_SNOVA] __attribute__((aligned(32))) = {0};
        uint64_t Aalpha_u64[lsq_SNOVA] __attribute__((aligned(32))) = {0};
        uint64_t Balpha_u64[lsq_SNOVA] __attribute__((aligned(32))) = {0};
        uint64_t Qalpha1_u64[lsq_SNOVA] __attribute__((aligned(32))) = {0};
        uint64_t Qalpha2_u64[lsq_SNOVA] __attribute__((aligned(32))) = {0};

        convert_GF16s_to_bytes((uint8_t *)Aalpha_u64, (uint8_t *)Aalpha, lsq_SNOVA * lsq_SNOVA);
        convert_GF16s_to_bytes((uint8_t *)Balpha_u64, (uint8_t *)Balpha, lsq_SNOVA * lsq_SNOVA);
        convert_GF16s_to_bytes((uint8_t *)Qalpha1_u64, (uint8_t *)Qalpha1, lsq_SNOVA * lsq_SNOVA);
        convert_GF16s_to_bytes((uint8_t *)Qalpha2_u64, (uint8_t *)Qalpha2, lsq_SNOVA * lsq_SNOVA);

        // evaluate the vinegar part of central map
        for (int index = 0; index < v_SNOVA; ++index)
        {
            gf16m_t vinegar_gf16_t;
            gf16m_transpose(vinegar_gf16[index], vinegar_gf16_t);

            uint64_t vinegar_gf16_u64;
            uint64_t vinegar_gf16_t_u64;
            convert_GF16s_to_bytes((uint8_t *)&vinegar_gf16_u64, (uint8_t *)vinegar_gf16[index], lsq_SNOVA);
            convert_GF16s_to_bytes((uint8_t *)&vinegar_gf16_t_u64, (uint8_t *)vinegar_gf16_t, lsq_SNOVA);
            for (int alpha = 0; alpha < lsq_SNOVA; alpha += 4)
            {
                uint64_t t256[4] __attribute__((aligned(32)));
                gf16m_u64_mul_4x4_4way_14(vinegar_gf16_t_u64, Qalpha1_u64 + alpha, t256);
                gf16m_u64_mul_4x4_4way_44(Aalpha_u64 + alpha, t256, Left_u64[index] + alpha);

                gf16m_u64_mul_4x4_4way_14(vinegar_gf16_u64, Balpha_u64 + alpha, t256);
                gf16m_u64_mul_4x4_4way_44(Qalpha2_u64 + alpha, t256, Right_u64[index] + alpha);
            }
        }

        convert_bytes_to_GF16s((uint8_t *)Left_u64, left_b, vl_SNOVA * lsq_SNOVA * l_SNOVA);
        convert_bytes_to_GF16s((uint8_t *)Right_u64, right_b, vl_SNOVA * lsq_SNOVA * l_SNOVA);

        // Main multiplication loop
        for (int jdx = 0; jdx < v_SNOVA; jdx++)
            for (int k1 = 0; k1 < l_SNOVA; k1++) {
                for (int mi = 0; mi < m_SNOVA; mi++)
                    for (int kdx = 0; kdx < v_SNOVA; kdx++)
                        {
                            uint32_t *target = (uint32_t *)&f11_8[mi * v_SNOVA * l_SNOVA + kdx * l_SNOVA];
                            uint32_t *source = (uint32_t *)&F11[mi][jdx][kdx][k1 * l_SNOVA];
                            *target = *source;
                        }

                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
                    for (int i1 = 0; i1 < l_SNOVA; i1 += 2)
                    {
                        int alpha_i1 = alpha * l_SNOVA + i1;

                        __m256i k_lh =
                            vtl_multtab(left_b[jdx * lsq_SNOVA * lsq_SNOVA + alpha_i1 * l_SNOVA + k1]) ^
                            _mm256_slli_epi16(
                                vtl_multtab(left_b[jdx * lsq_SNOVA * lsq_SNOVA + alpha_i1 * l_SNOVA + l_SNOVA + k1]), 4);
                        for (int mi_kdx_j1 = 0; mi_kdx_j1 < mvl_SNOVA32; mi_kdx_j1++)
                            temp3_256[alpha_i1 * mvl_SNOVA32 + mi_kdx_j1] ^=
                                _mm256_shuffle_epi8(k_lh, f11_256[mi_kdx_j1]);
                    }
            }

        // Separate twin result
        for (int alpha_i1 = 0; alpha_i1 < l_SNOVA * lsq_SNOVA; alpha_i1 += 2)
            for (int mi_kdx_j1 = 0; mi_kdx_j1 < mvl_SNOVA32; mi_kdx_j1++)
            {
                __m256i res = temp3_256[alpha_i1 * mvl_SNOVA32 + mi_kdx_j1];
                temp3_256[alpha_i1 * mvl_SNOVA32 + mi_kdx_j1] = res & l_mask;
                temp3_256[(alpha_i1 + 1) * mvl_SNOVA32 + mi_kdx_j1] = _mm256_srli_epi16(res, 4) & l_mask;
            }

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
            for (int mi = 0; mi < m_SNOVA; mi++)
                for (int kdx = 0; kdx < v_SNOVA; kdx++)
                {
                    uint32_t *res32 =
                        (uint32_t *)&temp3_8[alpha * l_SNOVA * mvl_SNOVA + mi * v_SNOVA * l_SNOVA + kdx * l_SNOVA];
                    __m128i val128 = _mm_setr_epi32(res32[0 * mvl_SNOVA / 4], res32[1 * mvl_SNOVA / 4],
                                                    res32[2 * mvl_SNOVA / 4], res32[3 * mvl_SNOVA / 4]);
                    uint32_t tres[4];
                    _mm_storeu_si128((__m128i_u *)&tres, _mm_shuffle_epi8(val128, transpose4x4));

                    for (int j1 = 0; j1 < l_SNOVA; j1++)
                        fvv_temp32[(j1 * v_SNOVA * ml_SNOVA + kdx * ml_SNOVA + mi * l_SNOVA) / 4] = tres[j1];
                }

            // Right multiplication
            for (int jdx = 0; jdx < v_SNOVA; ++jdx)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int k1 = 0; k1 < l_SNOVA; ++k1)
                    {
                        __m256i k_lh =
                            vtl_multtab(right_b[jdx * lsq_SNOVA * lsq_SNOVA + alpha * lsq_SNOVA + k1 * l_SNOVA + j1]);
                        for (int mi_i1 = 0; mi_i1 < ml_SNOVA32; mi_i1++)
                            fvv_res256[j1 * ml_SNOVA32 + mi_i1] ^=
                                _mm256_shuffle_epi8(k_lh, fvv_temp256[k1 * v_SNOVA * ml_SNOVA32 + jdx * ml_SNOVA32 + mi_i1]);
                    }
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
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        {
            memset(temp_l8, 0, mol_SNOVA * l_SNOVA);
            memset(res_left256, 0, mol_SNOVA * l_SNOVA);

            for (int i1 = 0; i1 < l_SNOVA; i1++)
                for (int jdx = 0; jdx < v_SNOVA; jdx++)
                    for (int k1 = 0; k1 < l_SNOVA; k1++)
                    {
                        __m256i k_lh = vtl_multtab(left_b[((jdx * lsq_SNOVA + alpha) * lsq_SNOVA + (i1 * l_SNOVA + k1))]);
                        for (int mi_kdx_j1 = 0; mi_kdx_j1 < mol_SNOVA32; mi_kdx_j1++)
                            temp_l256[i1 * mol_SNOVA32 + mi_kdx_j1] ^=
                                _mm256_shuffle_epi8(k_lh, f12_256[(jdx * l_SNOVA + k1) * mol_SNOVA32 + mi_kdx_j1]);
                    }

            // Shuffle Left
            for (int mi = 0; mi < m_SNOVA; ++mi)
                for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                {
                    uint32_t *res32 =
                        (uint32_t *)&temp_l8[mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA];
                    __m128i val128 = _mm_setr_epi32(res32[0 * mol_SNOVA / 4], res32[1 * mol_SNOVA / 4],
                                                    res32[2 * mol_SNOVA / 4], res32[3 * mol_SNOVA / 4]);
                    uint32_t tres[4];
                    _mm_storeu_si128((__m128i_u *)&tres, _mm_shuffle_epi8(val128, transpose4x4));
                    uint32_t *temp_left32 = (uint32_t *)&temp_left8[mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA];

                    for (int j1 = 0; j1 < l_SNOVA; j1++)
                        temp_left32[j1 * mol_SNOVA / 4] = tres[j1];
                }

            // O^2 * L^5
            for (int k1 = 0; k1 < l_SNOVA; ++k1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                {
                    __m256i k_lh = mtk2_16[Qalpha2[alpha][k1 * l_SNOVA + j1]];
                    for (int mi_kdx_i1 = 0; mi_kdx_i1 < mol_SNOVA32; mi_kdx_i1++)
                        res_left256[j1 * mol_SNOVA32 + mi_kdx_i1] ^=
                            _mm256_shuffle_epi8(k_lh, temp_left256[k1 * mol_SNOVA32 + mi_kdx_i1]);
                }

            // Outer product, O^2 * L^6
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                for (int j2 = 0; j2 < l_SNOVA; ++j2)
                {
                    __m256i k_lh = mtk2_16[Balpha[alpha][j2 * l_SNOVA + j1]];
                    for (int i2_mi_kdx_i1 = 0; i2_mi_kdx_i1 < l_SNOVA * mol_SNOVA32; i2_mi_kdx_i1++)
                        gtemp256[(j2 * l_SNOVA + j1) * l_SNOVA * mol_SNOVA32 + i2_mi_kdx_i1] ^=
                            _mm256_shuffle_epi8(k_lh, res_left256[i2_mi_kdx_i1]);
                }
        }

        // Shuffle to Gauss matrix
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int i2 = 0; i2 < l_SNOVA; ++i2)
                    {
                        uint32_t *res32 = (uint32_t *)&gtemp8[j1 * l_SNOVA * mol_SNOVA + i2 * mol_SNOVA +
                                                              mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA];
                        __m128i val128 =
                            _mm_setr_epi32(res32[0 * lsq_SNOVA * mol_SNOVA / 4], res32[1 * lsq_SNOVA * mol_SNOVA / 4],
                                           res32[2 * lsq_SNOVA * mol_SNOVA / 4], res32[3 * lsq_SNOVA * mol_SNOVA / 4]);
                        uint32_t tres[4];
                        _mm_storeu_si128((__m128i_u *)&tres, _mm_shuffle_epi8(val128, transpose4x4));

                        for (int i1 = 0; i1 < l_SNOVA; i1++)
                        {
                            uint32_t *temp_left32 =
                                (uint32_t *)&Gauss[mi * lsq_SNOVA + i1 * l_SNOVA + j1][kdx * lsq_SNOVA + i2 * l_SNOVA];
                            *temp_left32 = tres[i1];
                        }
                    }

        // Same for Right
        memset(gtemp256, 0, sizeof(gtemp8));

        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha)
        {
            memset(temp_r8, 0, mol_SNOVA * l_SNOVA);
            memset(res_right256, 0, mol_SNOVA * l_SNOVA);

            for (int i1 = 0; i1 < l_SNOVA; i1++)
                for (int jdx = 0; jdx < v_SNOVA; jdx++)
                    for (int k1 = 0; k1 < l_SNOVA; k1++)
                    {
                        int jdx_k1 = jdx * l_SNOVA + k1;

                        __m256i k_lh = vtl_multtab(right_b[(jdx * lsq_SNOVA + alpha) * lsq_SNOVA + k1 * l_SNOVA + i1]);
                        for (int mi_kdx_j1 = 0; mi_kdx_j1 < mol_SNOVA32; mi_kdx_j1++)
                            temp_r256[i1 * mol_SNOVA32 + mi_kdx_j1] ^=
                                _mm256_shuffle_epi8(k_lh, f21_256[jdx_k1 * mol_SNOVA32 + mi_kdx_j1]);
                    }

            // Shuffle
            for (int mi = 0; mi < m_SNOVA; ++mi)
                for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                {
                    uint32_t *res32 =
                        (uint32_t *)&temp_r8[mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA];
                    __m128i val128 = _mm_setr_epi32(res32[0 * mol_SNOVA / 4], res32[1 * mol_SNOVA / 4],
                                                    res32[2 * mol_SNOVA / 4], res32[3 * mol_SNOVA / 4]);
                    uint32_t tres[4];
                    _mm_storeu_si128((__m128i_u *)&tres, _mm_shuffle_epi8(val128, transpose4x4));

                    uint32_t *temp_left32 =
                        (uint32_t *)&temp_right8[mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA];
                    for (int j1 = 0; j1 < l_SNOVA; j1++)
                    {
                        temp_left32[j1 * mol_SNOVA / 4] = tres[j1];
                    }
                }

            // O^2 * L^5
            for (int k1 = 0; k1 < l_SNOVA; ++k1)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                {
                    __m256i k_lh = mtk2_16[Qalpha1[alpha][i1 * l_SNOVA + k1]];
                    for (int mi_kdx_j1 = 0; mi_kdx_j1 < mol_SNOVA32; mi_kdx_j1++)
                        res_right256[i1 * mol_SNOVA32 + mi_kdx_j1] ^=
                            _mm256_shuffle_epi8(k_lh, temp_right256[k1 * mol_SNOVA32 + mi_kdx_j1]);
                }

            // Outer product, O^2 * L^6
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j2 = 0; j2 < l_SNOVA; ++j2)
                {
                    __m256i k_lh = mtk2_16[Aalpha[alpha][i1 * l_SNOVA + j2]];
                    for (int i2_mi_kdx_j1 = 0; i2_mi_kdx_j1 < l_SNOVA * mol_SNOVA32; i2_mi_kdx_j1++)
                        gtemp256[(j2 * l_SNOVA + i1) * l_SNOVA * mol_SNOVA32 + i2_mi_kdx_j1] ^=
                            _mm256_shuffle_epi8(k_lh, res_right256[i2_mi_kdx_j1]);
                }
        }

        // Store in Gauss
        for (int mi = 0; mi < m_SNOVA; ++mi)
            for (int kdx = 0; kdx < o_SNOVA; ++kdx)
                for (int i1 = 0; i1 < l_SNOVA; ++i1)
                    for (int i2 = 0; i2 < l_SNOVA; ++i2)
                    {
                        uint32_t *res32 = (uint32_t *)&gtemp8[i1 * l_SNOVA * mol_SNOVA + i2 * mol_SNOVA +
                                                              mi * o_SNOVA * l_SNOVA + kdx * l_SNOVA];
                        __m128i val128 =
                            _mm_setr_epi32(res32[0 * lsq_SNOVA * mol_SNOVA / 4], res32[1 * lsq_SNOVA * mol_SNOVA / 4],
                                           res32[2 * lsq_SNOVA * mol_SNOVA / 4], res32[3 * lsq_SNOVA * mol_SNOVA / 4]);
                        uint32_t tres[4];
                        _mm_storeu_si128((__m128i_u *)&tres, _mm_shuffle_epi8(val128, transpose4x4));

                        for (int j1 = 0; j1 < l_SNOVA; j1++)
                        {
                            uint32_t *temp_left32 =
                                (uint32_t *)&Gauss[mi * lsq_SNOVA + i1 * l_SNOVA + j1][kdx * lsq_SNOVA + i2 * l_SNOVA];
                            *temp_left32 ^= tres[j1];
                        }
                    }

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
            SNOVA_CLEAR(vinegar_in_byte);
            SNOVA_CLEAR(left_b)
            SNOVA_CLEAR(right_b)
            SNOVA_CLEAR(temp3_8)
            SNOVA_CLEAR(temp_l8)
            SNOVA_CLEAR(temp_r8)
            SNOVA_CLEAR( temp_left8)
            SNOVA_CLEAR(temp_right8)
            SNOVA_CLEAR(gtemp8)
            SNOVA_CLEAR(fvv_temp8)
            SNOVA_CLEAR(fvv_res8)
            SNOVA_CLEAR(res_left256)
            SNOVA_CLEAR(res_right256)
        }
    } while (flag_redo);

    // printf("times of Gauss elimination : %d\n", num_sign);

    uint8_t t_GF16 = 0;
    uint8_t Gauss_last_col;
    for (int mil2 = m_SNOVA * lsq_SNOVA - 1; mil2 >= 0; --mil2)
    {
        Gauss_last_col = Gauss[mil2][m_SNOVA * lsq_SNOVA];
        __m256i t_GF16_256 = _mm256_setzero_si256();

        Gauss[mil2][m_SNOVA * lsq_SNOVA] = 0;

        int kstart = ((mil2 + 1) / 32) * 32;
        for (int k2 = kstart; k2 < GAUSS_ROW_mult32; k2 += 32)
            gf16_32_mul_32_add(Gauss[mil2] + k2, solution + k2, (uint8_t *)&t_GF16_256);

        uint8_t *t_GF16_256_8_ptr = (uint8_t *)(&t_GF16_256);
        for (int k2 = 0; k2 < 32; k2++)
            t_GF16 ^= t_GF16_256_8_ptr[k2];

        solution[mil2] = Gauss_last_col ^ t_GF16;
        t_GF16 = 0;
        Gauss_last_col = 0;
    }

    // Establish Signature

    for (int index = 0; index < o_SNOVA; ++index)
        for (int i1 = 0; i1 < l_SNOVA; ++i1)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                vinegar_gf16[index + v_SNOVA][i1 * l_SNOVA + j1] = solution[index * lsq_SNOVA + i1 * l_SNOVA + j1];

    uint64_t T12_u64_ov[o_SNOVA][v_SNOVA_mult4] __attribute__((aligned(32))) = {0};
    uint64_t T12_u64[v_SNOVA][o_SNOVA] __attribute__((aligned(32))) = {0};
    uint64_t sig_u64[n_SNOVA] __attribute__((aligned(32))) = {0};

    convert_GF16s_to_bytes((uint8_t *)T12_u64, (uint8_t *)T12, v_SNOVA * o_SNOVA * lsq_SNOVA);
    for (int i = 0; i < o_SNOVA; ++i)
        for (int j = 0; j < v_SNOVA; ++j)
            T12_u64_ov[i][j] = T12_u64[j][i];

    convert_GF16s_to_bytes((uint8_t *)sig_u64, (uint8_t *)vinegar_gf16, n_SNOVA * lsq_SNOVA);

    for (int i = 0; i < o_SNOVA; ++i)
        for (int index = 0; index < v_SNOVA; index += 4)
            gf16m_u64_mul_4x4_4way_41_add(T12_u64_ov[i] + index, sig_u64[v_SNOVA + i], (__m256i *)(sig_u64 + index));

    memcpy(pt_signature, sig_u64, bytes_signature);

    // output signature
    memcpy(pt_signature + bytes_signature, array_salt, bytes_salt);

    // Cleanup
    SNOVA_CLEAR(vinegar_gf16)
    SNOVA_CLEAR(solution)
    SNOVA_CLEAR(hash_in_GF16)
    SNOVA_CLEAR(signed_hash)
    SNOVA_CLEAR(Gauss)
    SNOVA_CLEAR(f11_8)
    SNOVA_CLEAR(f12_8)
    SNOVA_CLEAR(f21_8)
    SNOVA_CLEAR(T12_u64_ov)
    SNOVA_CLEAR(T12_u64)

    return 0;
}

#endif