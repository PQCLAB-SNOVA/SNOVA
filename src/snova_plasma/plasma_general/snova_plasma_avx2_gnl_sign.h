/**
 * VTL version for any l_SNOVA.
 */

#ifndef PLASMA_GNL_AVX2_SIGN_H
#define PLASMA_GNL_AVX2_SIGN_H

#include <immintrin.h>
#include <stdint.h>
#include "../../snova.h"
#include "../snova_plasma_avx2.h"

#define GAUSS_ROW (m_SNOVA * lsq_SNOVA + 1)
#define GAUSS_ROW32 ((GAUSS_ROW + 31) / 32)
#define GAUSS_ROW_mult32 (GAUSS_ROW32 * 32)


void calc_LR_J_vtl(
    uint8_t L_J_nibble[m_SNOVA][alpha_SNOVA][rank_next2 / 2][vtl_v_len * 32],
    uint8_t R_tr_J[m_SNOVA][alpha_SNOVA][rank_next2][vtl_v_len * 32], 
    uint8_t R_tr_J_nibble[m_SNOVA][alpha_SNOVA][rank_next2 / 2][vtl_v_len * 32], 
    Aalpha_t Aalpha, 
    Balpha_t Balpha, 
    Qalpha1_t Qalpha1,
    Qalpha2_t Qalpha2, 
    gf16m_t* X_in_GF16Matrix) {

    __m256i X_J[rank][vtl_v_len] = {0};
    __m256i X_tr_J[rank][vtl_v_len] = {0};

    jogressMatrix_avx2((uint8_t *)X_J, (uint8_t *)X_in_GF16Matrix, 1, v_SNOVA);
    jogressTrMatrix_avx2((uint8_t *)X_tr_J, (uint8_t *)X_in_GF16Matrix, 1, v_SNOVA);

    // calc LR
    for (int mi = 0; mi < m_SNOVA; ++mi) {
        __m256i AxS_tr_256[alpha_SNOVA_next2][rank][vtl_v_len]  = {0};
        __m256i Q2xS_tr_256[alpha_SNOVA_next2][rank][vtl_v_len] = {0};

        // A and Q2, vtl 2 alpha set end 0
        gf16m_t Aalpha_vtl[alpha_SNOVA_next2] = {0};
        gf16m_t Qalpha2_vtl[alpha_SNOVA_next2] = {0};
        memcpy((uint8_t *)Aalpha_vtl, (uint8_t *)(Aalpha[mi]), sizeof(Aalpha[mi]));
        memcpy((uint8_t *)Qalpha2_vtl, (uint8_t *)(Qalpha2[mi]), sizeof(Qalpha2[mi]));

        for (int alpha = 0; alpha < alpha_SNOVA; alpha += 2) {
            __m256i AxS_256[rank][vtl_v_len] = {0};
            __m256i Q2xS_256[rank][vtl_v_len] = {0};
            
            for (int ni = 0; ni < rank; ++ni) {
                for (int nj = 0; nj < rank; ++nj) {
                    __m256i k1_lh = mtk2_16[get_gf16m(Aalpha_vtl[alpha], ni, nj) | (get_gf16m(Aalpha_vtl[alpha + 1], ni, nj) << 4)];
                    __m256i k2_lh = mtk2_16[get_gf16m(Qalpha2_vtl[alpha], ni, nj) | (get_gf16m(Qalpha2_vtl[alpha + 1], ni, nj) << 4)];
                    for (int nk = 0; nk < vtl_v_len; ++nk) {
                        AxS_256[ni][nk] ^= _mm256_shuffle_epi8(k1_lh, X_tr_J[nj][nk]);
                        Q2xS_256[ni][nk] ^= _mm256_shuffle_epi8(k2_lh, X_J[nj][nk]);
                    }
                }
            }
            jogressMatrixTr_avx2((uint8_t *)AxS_tr_256[alpha], (uint8_t *)AxS_256, 1, v_SNOVA);
            jogressMatrixTr_avx2((uint8_t *)Q2xS_tr_256[alpha], (uint8_t *)Q2xS_256, 1, v_SNOVA);

            // nibble splite
            for (int ni = 0; ni < rank; ++ni) {
                for (int nk = 0; nk < vtl_v_len; ++nk) {
                    AxS_tr_256[alpha + 1][ni][nk] = _mm256_srli_epi16(AxS_tr_256[alpha][ni][nk], 4) & _mm256_set1_epi8(0x0f);
                    AxS_tr_256[alpha][ni][nk] &= _mm256_set1_epi8(0x0f);

                    Q2xS_tr_256[alpha + 1][ni][nk] = _mm256_srli_epi16(Q2xS_tr_256[alpha][ni][nk], 4) & _mm256_set1_epi8(0x0f);
                    Q2xS_tr_256[alpha][ni][nk] &= _mm256_set1_epi8(0x0f);
                }
            }
            SNOVA_CLEAR(AxS_256);
            SNOVA_CLEAR(Q2xS_256);
        }

        for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
            int mj = i_prime(mi, alpha);
            __m256i L_tr_J_256[rank][vtl_v_len] = {0}; 
            for (int ni = 0; ni < rank_floor2; ni += 2) {
                __m256i* R_tr_J_256 = (__m256i *)R_tr_J[mj][alpha][ni];
                for (int nj = 0; nj < rank; ++nj) {
                    __m256i k1_lh = mtk2_16[get_gf16m(Qalpha1[mi][alpha], ni, nj) ^ (get_gf16m(Qalpha1[mi][alpha], ni + 1, nj) << 4)];
                    __m256i k2_lh = mtk2_16[get_gf16m(Balpha[mi][alpha], nj, ni) ^ (get_gf16m(Balpha[mi][alpha], nj, ni + 1) << 4)];
                    for (int nk = 0; nk < vtl_v_len; ++nk) {
                        L_tr_J_256[ni][nk] ^= _mm256_shuffle_epi8(k1_lh, AxS_tr_256[alpha][nj][nk]);
                        R_tr_J_256[nk] ^= _mm256_shuffle_epi8(k2_lh, Q2xS_tr_256[alpha][nj][nk]);
                    }
                }
            }

#if rank % 2    //  rank is odd, use vl no vtl(vectorized look-up), 
            __m256i* R_tr_J_256_last = (__m256i *)R_tr_J[mj][alpha][rank_floor2];
            for (int nj = 0; nj < rank; ++nj) {
                __m256i k1_lh = mtk2_16[get_gf16m(Qalpha1[mi][alpha], rank_floor2, nj)];
                __m256i k2_lh = mtk2_16[get_gf16m(Balpha[mi][alpha], nj, rank_floor2)];
                for (int nk = 0; nk < vtl_v_len; ++nk) {
                    L_tr_J_256[rank_floor2][nk] ^= _mm256_shuffle_epi8(k1_lh, AxS_tr_256[alpha][nj][nk]);
                    R_tr_J_256_last[nk] ^= _mm256_shuffle_epi8(k2_lh, Q2xS_tr_256[alpha][nj][nk]);
                }
            }
#endif

            // nibble splite
            for (int ni = 0; ni < rank_floor2; ni +=2) {
                for (int nk = 0; nk < vtl_v_len; ++nk) {
                    __m256i* R_tr_J_256_l = (__m256i *)R_tr_J[mj][alpha][ni];
                    __m256i* R_tr_J_256_h = (__m256i *)R_tr_J[mj][alpha][ni + 1];
                    L_tr_J_256[ni + 1][nk] = _mm256_srli_epi16(L_tr_J_256[ni][nk], 4) & _mm256_set1_epi8(0x0f);
                    L_tr_J_256[ni][nk] &= _mm256_set1_epi8(0x0f);

                    R_tr_J_256_h[nk] = _mm256_srli_epi16(R_tr_J_256_l[nk], 4) & _mm256_set1_epi8(0x0f);
                    R_tr_J_256_l[nk] &= _mm256_set1_epi8(0x0f);
                }
            }

            uint8_t L_J[rank_next2][vtl_v_len * 32] __attribute__((aligned(32))) = {0};
            jogressMatrixTr_avx2((uint8_t *)L_J, (uint8_t *)L_tr_J_256, 1, v_SNOVA);

            // nibble L_J.
            for (int ni = 0; ni < rank; ni += 2) {
                __m256i* L_J_l_256 = (__m256i *)L_J[ni];
                __m256i* L_J_h_256 = (__m256i *)L_J[ni + 1];
                __m256i* L_J_256 = (__m256i *)L_J_nibble[mj][alpha][ni / 2];
                for (int nk = 0; nk < vtl_v_len; ++nk) {
                    L_J_256[nk] = L_J_l_256[nk] ^ _mm256_slli_epi16(L_J_h_256[nk], 4);
                }

                __m256i* R_tr_J_l_256 = (__m256i *)R_tr_J[mj][alpha][ni];
                __m256i* R_tr_J_h_256 = (__m256i *)R_tr_J[mj][alpha][ni + 1];
                __m256i* R_tr_J_256 = (__m256i *)R_tr_J_nibble[mj][alpha][ni / 2];
                for (int nk = 0; nk < vtl_v_len; ++nk) {
                    R_tr_J_256[nk] = R_tr_J_l_256[nk] ^ _mm256_slli_epi16(R_tr_J_h_256[nk], 4);
                }
            }
            SNOVA_CLEAR(L_tr_J_256);
            SNOVA_CLEAR(L_J);
        }
        SNOVA_CLEAR(AxS_tr_256);
        SNOVA_CLEAR(Q2xS_tr_256);
    }
    SNOVA_CLEAR(X_J);
    SNOVA_CLEAR(X_tr_J);
}

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
    uint8_t Gauss[m_SNOVA * lsq_SNOVA][GAUSS_ROW_mult32] __attribute__((aligned(32)));

    gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
    gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA];
    uint8_t signed_hash[bytes_hash];
    uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1];

    int flag_redo = 1;
    uint8_t num_sign = 0;

    createSignedHash(digest, bytes_digest, pt_public_key_seed, array_salt, signed_hash);

    // put hash value in GF16 array
    convert_bytes_to_GF16s(signed_hash, hash_in_GF16, GF16s_hash);

    do {
        gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
        uint8_t L_J_nibble[m_SNOVA][alpha_SNOVA][rank_next2 / 2][vtl_v_len * 32] __attribute__((aligned(32))) = {0};
        uint8_t R_tr_J[m_SNOVA][alpha_SNOVA][rank_next2][vtl_v_len * 32] __attribute__((aligned(32))) = {0};
        uint8_t R_tr_J_nibble[m_SNOVA][alpha_SNOVA][rank_next2 / 2][vtl_v_len * 32] __attribute__((aligned(32))) = {0};
        gf16m_t F21_vo[v_SNOVA][o_SNOVA];
        uint8_t Temp1[o_SNOVA * l_SNOVA * lsq_SNOVA * vtl_mainRow_x_rank(o_SNOVA)] __attribute__((aligned(32))) = {0};
        uint8_t Temp2[o_SNOVA * l_SNOVA * lsq_SNOVA * vtl_mainRow_x_rank(o_SNOVA)] __attribute__((aligned(32))) = {0};
        __m256i Fvv_256[m_SNOVA][rank][rank] = {0};

        memset(Gauss, 0, sizeof(Gauss));
        num_sign++;
        if(num_sign == 255) {
            // Probability of getting here is about 2^{-1020}
            memset(pt_signature, 0, bytes_sig_with_salt);
            return -1;
        }
        flag_redo = 0;
        // put hash value in the last column of Gauss matrix
        for (int index = 0; index < (m_SNOVA * lsq_SNOVA); index++) {
            Gauss[index][m_SNOVA * lsq_SNOVA] = hash_in_GF16[index];
        }
        // generate the vinegar value
        Keccak_HashInstance hashInstance;
        Keccak_HashInitialize_SHAKE256(&hashInstance);
        Keccak_HashUpdate(&hashInstance, pt_private_key_seed, 8 * seed_length_private);
        Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
        Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
        Keccak_HashUpdate(&hashInstance, &num_sign, 8);
        Keccak_HashFinal(&hashInstance, NULL);
        Keccak_HashSqueeze(&hashInstance, vinegar_in_byte, 8 * ((v_SNOVA * lsq_SNOVA + 1) >> 1));

        convert_bytes_to_GF16s(vinegar_in_byte, (uint8_t *)X_in_GF16Matrix, v_SNOVA * lsq_SNOVA);
        calc_LR_J_vtl(L_J_nibble, R_tr_J, R_tr_J_nibble, Aalpha, Balpha, Qalpha1, Qalpha2, X_in_GF16Matrix);

        // ------ test END
        for (int mi = 0; mi < m_SNOVA; ++mi) {
            gf16m_set_zero(Fvv_in_GF16Matrix[mi]);
        }

        for (int mi = 0; mi < m_SNOVA; ++mi) {

            __m256i F11_J_256[v_SNOVA * rank][vtl_v_len] = {0};
            jogressMatrix_avx2((uint8_t *)F11_J_256, (uint8_t *)F11[mi], v_SNOVA, v_SNOVA);
            for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
                int mi_prime = i_prime_inv(mi, alpha);
                __m256i LJxF11J_256[rank][vtl_v_len] = {0};
                __m256i LJxF11J_256_nibble[rank_next2 / 2][vtl_v_len] = {0};

                for (int vi = 0; vi < rank_next2 / 2; ++vi) {
                    for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
                        __m256i k_lh = mtk2_16[L_J_nibble[mi][alpha][vi][vj]];
                        for (int vk = 0; vk < vtl_v_len; ++vk) {
                            LJxF11J_256_nibble[vi][vk] ^= _mm256_shuffle_epi8(k_lh, F11_J_256[vj][vk]);
                        }
                    }
                }

                for (int ni = 0; ni < rank_floor2; ni +=2) {
                    __m256i* LJxF11J_256_nibble_l = LJxF11J_256[ni];
                    __m256i* LJxF11J_256_nibble_h = LJxF11J_256[ni + 1];
                    for (int nj = 0; nj < vtl_v_len; ++nj) {
                        LJxF11J_256_nibble_l[nj] = LJxF11J_256_nibble[ni / 2][nj] & _mm256_set1_epi8(0x0f);
                        LJxF11J_256_nibble_h[nj] = (_mm256_srli_epi16(LJxF11J_256_nibble[ni / 2][nj], 4) & _mm256_set1_epi8(0x0f));
                    }
                }

#if rank % 2 // rank is odd
                for (int nj = 0; nj < vtl_v_len; ++nj) {
                    LJxF11J_256[rank_floor2][nj] = LJxF11J_256_nibble[rank_floor2 / 2][nj];
                }
#endif

                for (int vi = 0; vi < rank; ++vi) {   
                    for (int vj = 0; vj < rank; ++vj) {
                        __m256i* R_tr_J_256 = (__m256i *)R_tr_J[mi][alpha][vj];
                        __m256i tmp_256 = _mm256_setzero_si256();
                        for (int vk = 0; vk < vtl_v_len; ++vk) {
                            gf16_32_mul_32_add((uint8_t *)(LJxF11J_256[vi] + vk), (uint8_t *)(R_tr_J_256 + vk), (uint8_t *)(&tmp_256));
                        }

                        Fvv_256[mi_prime][vi][vj] ^= tmp_256;
                    }
                }
                SNOVA_CLEAR(LJxF11J_256);
                SNOVA_CLEAR(LJxF11J_256_nibble);
            }
            SNOVA_CLEAR(F11_J_256);
        }

        for (int mi = 0; mi < m_SNOVA; ++mi) {
            for (int ni = 0; ni < rank; ++ni) {   
                for (int nj = 0; nj < rank; ++nj) {
                    set_gf16m(Fvv_in_GF16Matrix[mi], ni, nj, horizontal_xor_256(Fvv_256[mi][ni][nj]));
                }
            }
        }

        // add to the last column of Gauss matrix
        for (int i = 0; i < m_SNOVA; ++i) {
            for (int j = 0; j < rank; ++j) {
                for (int k = 0; k < rank; ++k) {
                    int index1 = i * lsq_SNOVA + j * rank + k;
                    int index2 = m_SNOVA * lsq_SNOVA;
                    Gauss[index1][index2] = gf16_get_add(Gauss[index1][index2], get_gf16m(Fvv_in_GF16Matrix[i], j, k));
                }
            }
        }

        // compute the coefficients of Xo and put into Gauss matrix and compute
        // the coefficients of Xo^t and add into Gauss matrix
        for (int mi = 0; mi < m_SNOVA; ++mi) {
            __m256i F12_J_256[v_SNOVA * rank][vtl_o_len] = {0};
            __m256i F21_vo_tr_J_256[v_SNOVA * rank][vtl_o_len] = {0};
            
            // swap F21  v <-> o
            for (int oi = 0; oi < o_SNOVA; ++oi) {
                for (int vi = 0; vi < v_SNOVA; ++vi) {
                    gf16m_clone(F21_vo[vi][oi], F21[mi][oi][vi]);
                }
            }

            jogressMatrix_avx2((uint8_t *)F12_J_256, (uint8_t *)F12[mi], v_SNOVA, o_SNOVA);
            jogressTrMatrix_avx2((uint8_t *)F21_vo_tr_J_256, (uint8_t *)F21_vo, v_SNOVA, o_SNOVA);

            for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
                int mi_prime_inv = i_prime_inv(mi, alpha);

                __m256i LJxF12J_256[rank][vtl_o_len] = {0};
                __m256i LJxF12J_256_nibble[rank_next2 / 2][vtl_o_len] = {0};

                for (int vi = 0; vi < rank_next2 / 2; ++vi) {
                    for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
                        __m256i k_lh = mtk2_16[L_J_nibble[mi][alpha][vi][vj]];
                        for (int oi = 0; oi < vtl_o_len; ++oi) {
                            LJxF12J_256_nibble[vi][oi] ^= _mm256_shuffle_epi8(k_lh, F12_J_256[vj][oi]);
                        }
                    }
                }

                // nibble splite, **PS. If rank is odd, the last value is handled separately.**
                for (int vi = 0; vi < rank_floor2; vi +=2) {
                    __m256i* LJxF12J_256_nibble_l = LJxF12J_256[vi];
                    __m256i* LJxF12J_256_nibble_h = LJxF12J_256[vi + 1];
                    for (int vj = 0; vj < vtl_o_len; ++vj) {
                        LJxF12J_256_nibble_l[vj] = LJxF12J_256_nibble[vi / 2][vj] & _mm256_set1_epi8(0x0f);
                        LJxF12J_256_nibble_h[vj] = (_mm256_srli_epi16(LJxF12J_256_nibble[vi / 2][vj], 4) & _mm256_set1_epi8(0x0f));
                    }
                }

#if rank % 2 // rank is odd
                for (int nj = 0; nj < vtl_o_len; ++nj) {
                    LJxF12J_256[rank_floor2][nj] = LJxF12J_256_nibble[rank_floor2 / 2][nj];
                }
#endif


                __m256i LJxF12J_tr_256[rank][vtl_o_len] = {0};
                __m256i LJxF12JxQJ_tr_256[rank][vtl_o_len] = {0};

                jogressMatrixTr_avx2((uint8_t *)LJxF12J_tr_256, (uint8_t *)LJxF12J_256, 1, o_SNOVA);

                for (int ri = 0; ri < rank; ++ri) {
                    for (int rj = 0; rj < rank; ++rj) {
                        __m256i k = mtk2_16[get_gf16m(Qalpha2[mi_prime_inv][alpha], rj, ri)];
                        for (int oi = 0; oi < vtl_o_len; ++oi) {
                            LJxF12JxQJ_tr_256[ri][oi] ^= _mm256_shuffle_epi8(k, LJxF12J_tr_256[rj][oi]);
                        }
                    }
                }

                for (int ti2 = 0; ti2 < l_SNOVA; ++ti2)
                    for (int tj2 = 0; tj2 < l_SNOVA; ++tj2)
                    {
                        __m256i k = mtk2_16[Balpha[mi_prime_inv][alpha][tj2 * rank + ti2]];
                        __m256i *Temp1_256 = (__m256i *)Temp1;
                        __m256i *AJ_256 = (__m256i *)LJxF12JxQJ_tr_256;

                        for (int tj1 = 0; tj1 < l_SNOVA; tj1++)
                            for (int toi = 0; toi < vtl_mainRow_x_rank32(o_SNOVA); toi++)
                                {
                                    Temp1_256[((mi_prime_inv * lsq_SNOVA + ti2 * rank + tj2) * l_SNOVA + tj1) * vtl_mainRow_x_rank32(o_SNOVA) + toi] ^= 
                                        _mm256_shuffle_epi8(k, AJ_256[(tj1 * vtl_mainRow_x_rank32(o_SNOVA)) + toi]);
                                }
                    }

                // ------- ^^^ F12    F22 vvv -------
                __m256i F21JxRJ_256[rank][vtl_o_len] = {0};

                __m256i RTRJ_x_F21TRJ_256[rank][vtl_o_len] = {0};
                __m256i RTRJ_x_F21TRJ_256_nibble[rank_next2 / 2][vtl_o_len] = {0};

                for (int vi = 0; vi < rank_next2 / 2; ++vi) {
                    for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
                        __m256i k_lh = mtk2_16[R_tr_J_nibble[mi][alpha][vi][vj]];
                        for (int oi = 0; oi < vtl_o_len; ++oi) {
                            RTRJ_x_F21TRJ_256_nibble[vi][oi] ^= _mm256_shuffle_epi8(k_lh, F21_vo_tr_J_256[vj][oi]);
                        }
                    }
                }

                for (int vi = 0; vi < rank_floor2; vi +=2) {
                    __m256i* RTRJ_x_F21TRJ_256_nibble_l = RTRJ_x_F21TRJ_256[vi];
                    __m256i* RTRJ_x_F21TRJ_256_nibble_h = RTRJ_x_F21TRJ_256[vi + 1];
                    for (int vj = 0; vj < vtl_o_len; ++vj) {
                        RTRJ_x_F21TRJ_256_nibble_l[vj] = RTRJ_x_F21TRJ_256_nibble[vi / 2][vj] & _mm256_set1_epi8(0x0f);
                        RTRJ_x_F21TRJ_256_nibble_h[vj] = (_mm256_srli_epi16(RTRJ_x_F21TRJ_256_nibble[vi / 2][vj], 4) & _mm256_set1_epi8(0x0f));
                    }
                }

#if rank % 2 // rank is odd
                for (int nj = 0; nj < vtl_o_len; ++nj) {
                    RTRJ_x_F21TRJ_256[rank_floor2][nj] = RTRJ_x_F21TRJ_256_nibble[rank_floor2 / 2][nj];
                }
#endif

                jogressMatrixTr_avx2((uint8_t *)F21JxRJ_256, (uint8_t *)RTRJ_x_F21TRJ_256, 1, o_SNOVA);

                __m256i Q1xF21JxRJ_256[rank][vtl_o_len] = {0};
                for (int ri = 0; ri < rank; ++ri) {
                    for (int rj = 0; rj < rank; ++rj) {
                        __m256i k = mtk2_16[get_gf16m(Qalpha1[mi_prime_inv][alpha], ri, rj)];
                        for (int oi = 0; oi < vtl_o_len; ++oi) {
                            Q1xF21JxRJ_256[ri][oi] ^= _mm256_shuffle_epi8(k, F21JxRJ_256[rj][oi]);
                        }
                    }
                }

                for (int ti1 = 0; ti1 < l_SNOVA; ++ti1)
                    for (int tj2 = 0; tj2 < l_SNOVA; ++tj2)
                    {
                        __m256i k = mtk2_16[Aalpha[mi_prime_inv][alpha][ti1 * rank + tj2]];
                        __m256i *Temp2_256 = (__m256i *)Temp2;
                        __m256i *AJ_256 = (__m256i *)Q1xF21JxRJ_256;

                        for (int tj1 = 0; tj1 < l_SNOVA; ++tj1)
                            for (int toi = 0; toi < vtl_mainRow_x_rank32(o_SNOVA); toi++)
                                {
                                    Temp2_256[((mi_prime_inv * lsq_SNOVA + ti1 * l_SNOVA + tj2) * l_SNOVA + tj1) * vtl_mainRow_x_rank32(o_SNOVA) + toi] ^= 
                                        _mm256_shuffle_epi8(k, AJ_256[(tj1 * vtl_mainRow_x_rank32(o_SNOVA)) + toi]);
                                }
                    }
                
                SNOVA_CLEAR(LJxF12J_256);
                SNOVA_CLEAR(LJxF12J_256_nibble);
                SNOVA_CLEAR(LJxF12J_tr_256);
                SNOVA_CLEAR(LJxF12JxQJ_tr_256);

                SNOVA_CLEAR(F21JxRJ_256);
                SNOVA_CLEAR(RTRJ_x_F21TRJ_256);
                SNOVA_CLEAR(RTRJ_x_F21TRJ_256_nibble);
                SNOVA_CLEAR(Q1xF21JxRJ_256);

            }

            SNOVA_CLEAR(F12_J_256);
            SNOVA_CLEAR(F21_vo_tr_J_256);
        }

        for (int mi_prime_inv = 0; mi_prime_inv < o_SNOVA; ++mi_prime_inv)
            for (int ti1 = 0; ti1 < l_SNOVA; ++ti1)
                for (int ti2 = 0; ti2 < l_SNOVA; ++ti2)
                    for (int tj1 = 0; tj1 < l_SNOVA; ++tj1)
                        for (int oi = 0; oi < o_SNOVA; ++oi)
                            for (int tj2 = 0; tj2 < l_SNOVA; ++tj2)
                                Gauss[mi_prime_inv * lsq_SNOVA + ti1 * rank + ti2][oi * lsq_SNOVA + tj1 * rank + tj2] ^= 
                                    Temp1[((mi_prime_inv * lsq_SNOVA + ti2 * rank + tj2) * l_SNOVA + tj1) * vtl_mainRow_x_rank(o_SNOVA) + (oi * rank + ti1)];

        for (int mi_prime_inv = 0; mi_prime_inv < o_SNOVA; ++mi_prime_inv)
            for (int oi = 0; oi < o_SNOVA; ++oi)
                for (int ti1 = 0; ti1 < l_SNOVA; ++ti1)
                    for (int ti2 = 0; ti2 < l_SNOVA; ++ti2)
                        for (int tj1 = 0; tj1 < l_SNOVA; ++tj1) 
                            for (int tj2 = 0; tj2 < l_SNOVA; ++tj2)
                                Gauss[mi_prime_inv * lsq_SNOVA + ti1 * rank + ti2][oi * lsq_SNOVA + tj1 * rank + tj2] ^= 
                                    Temp2[((mi_prime_inv * lsq_SNOVA + ti1 * l_SNOVA + tj2) * l_SNOVA + tj1) * vtl_mainRow_x_rank(o_SNOVA) + (oi * rank + ti2)];

        // Gaussian elimination in constant time
        for (int mi2 = 0; mi2 < m_SNOVA * lsq_SNOVA; ++mi2)
        {
            int swap = ct_gf16_is_not_zero(Gauss[mi2][mi2]) - 1;
            for (int j2 = mi2 + 1; j2 < m_SNOVA * lsq_SNOVA; ++j2)
            {
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

        if (!flag_redo) {
            SNOVA_CLEAR(Fvv_in_GF16Matrix);
            SNOVA_CLEAR(L_J_nibble);
            SNOVA_CLEAR(R_tr_J);
            SNOVA_CLEAR(R_tr_J_nibble);
            SNOVA_CLEAR(Temp1);
            SNOVA_CLEAR(Temp2);
            SNOVA_CLEAR(F21_vo);
            SNOVA_CLEAR(Fvv_256);
        }
    } while (flag_redo);

    uint8_t t_GF16 = 0;
    uint8_t Gauss_last_col;
    uint8_t solution[GAUSS_ROW_mult32] = {0};

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

    uint32_t xSolution[m_SNOVA * lsq_SNOVA] = {0};
    uint32_t xT12[v_SNOVA][o_SNOVA][lsq_SNOVA] = {0};

    for (int dj = 0; dj < v_SNOVA; ++dj)
        for (int dk = 0; dk < o_SNOVA; ++dk)
            for (int idx = 0; idx < lsq_SNOVA; ++idx)
                xT12[dj][dk][idx] = gf16_from_nibble(T12[dj][dk][idx]);

    for (int index = 0; index < o_SNOVA; ++index)
        for (int i1 = 0; i1 < l_SNOVA; ++i1)
            for (int j1 = 0; j1 < l_SNOVA; ++j1)
                X_in_GF16Matrix[index + v_SNOVA][i1 * l_SNOVA + j1] = solution[index * lsq_SNOVA + i1 * l_SNOVA + j1];

    for (int mi2 = m_SNOVA * lsq_SNOVA - 1; mi2 >= 0; --mi2)
        xSolution[mi2] = gf16_from_nibble(solution[mi2]);

    // Establish Signature
    uint32_t xSig[lsq_SNOVA] = {0};
    for (int dj = 0; dj < v_SNOVA; ++dj)
    {
        for (int dk = 0; dk < o_SNOVA; ++dk)
            for (int i1 = 0; i1 < l_SNOVA; ++i1)
                for (int j1 = 0; j1 < l_SNOVA; ++j1)
                    for (int k1 = 0; k1 < l_SNOVA; ++k1)
                        xSig[i1 * l_SNOVA + j1] ^= xT12[dj][dk][i1 * l_SNOVA + k1] * xSolution[dk * lsq_SNOVA + k1 * l_SNOVA + j1];

        for (int idx = 0; idx < lsq_SNOVA; ++idx)
            signature_in_GF16Matrix[dj][idx] = X_in_GF16Matrix[dj][idx] ^ gf16_to_nibble(xSig[idx]);

        memset(xSig, 0, sizeof(xSig));
    }

    for (int index = 0; index < o_SNOVA; ++index)
        for (int idx = 0; idx < lsq_SNOVA; ++idx)
            signature_in_GF16Matrix[v_SNOVA + index][idx] = X_in_GF16Matrix[v_SNOVA + index][idx];

    // output signature
    convert_GF16s_to_bytes(pt_signature, (gf16_t *)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);
    for (int i1 = 0; i1 < bytes_salt; ++i1)
        pt_signature[bytes_signature + i1] = array_salt[i1];
    
    // Cleanup
    SNOVA_CLEAR(Gauss);
    SNOVA_CLEAR(X_in_GF16Matrix);
    SNOVA_CLEAR(hash_in_GF16);
    SNOVA_CLEAR(signature_in_GF16Matrix);
    SNOVA_CLEAR(signed_hash);
    SNOVA_CLEAR(vinegar_in_byte);

    SNOVA_CLEAR(xSolution);
    SNOVA_CLEAR(xT12);
    
    return 0;
}

#endif
