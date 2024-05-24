/*
Use uint64_t to store a GF16 4x4 matrix (with each half-byte storing one element of GF16)."
Optimizing GF16 4x4 matrix multiplication using AVX2.
*/

#ifndef PLASMA_4x4_AVX2_VERIFY_H
#define PLASMA_4x4_AVX2_VERIFY_H

#include <immintrin.h>
#include <stdint.h>

#include "../../deriv_params.h"
#include "../../gf16.h"
#include "../../snova_kernel.h"
#include "../snova_plasma_avx2.h"
// uint64_t to store a GF16 4x4 matrix
#include "snova_plasma_4x4_avx2.h"
#include "snova_plasma_4x4_data.h"

static inline void convert_4x4x4_to_16x4(__m256i* restrict out, __m256i* restrict in) {
    __m256i shf1 = _mm256_setr_epi32(0, 4, 2, 6, 1, 5, 3, 7);
    __m256i shf2 = _mm256_setr_epi8(0, 1, 8, 9, 4, 5, 12, 13, 2, 3, 10, 11, 6, 7, 14, 15,  //
                                    0, 1, 8, 9, 4, 5, 12, 13, 2, 3, 10, 11, 6, 7, 14, 15);

    __m256i temp[4];
    uint64_t* temp_u64 = (uint64_t*)temp;

    for (int j = 0; j < 4; j++) {
        temp[j] = _mm256_permutevar8x32_epi32(in[j], shf1);
        temp[j] = _mm256_shuffle_epi8(temp[j], shf2);
    }

    for (int j = 0; j < 4; j++) {
        out[j] = _mm256_setr_epi64x(temp_u64[j], temp_u64[j + 4], temp_u64[j + 8], temp_u64[j + 12]);
    }
}

static inline void convert_16x4_to_4x4x4(__m256i* restrict out, __m256i* restrict in) {
    __m256i shf1 = _mm256_setr_epi32(0, 4, 2, 6, 1, 5, 3, 7);
    __m256i shf2 = _mm256_setr_epi8(0, 1, 8, 9, 4, 5, 12, 13, 2, 3, 10, 11, 6, 7, 14, 15,  //
                                    0, 1, 8, 9, 4, 5, 12, 13, 2, 3, 10, 11, 6, 7, 14, 15);

    uint64_t* in_u64 = (uint64_t*)in;
    for (int j = 0; j < 4; j++) {
        out[j] = _mm256_setr_epi64x(in_u64[j], in_u64[j + 4], in_u64[j + 8], in_u64[j + 12]);
    }

    for (int j = 0; j < 4; j++) {
        out[j] = _mm256_shuffle_epi8(out[j], shf2);
        out[j] = _mm256_permutevar8x32_epi32(out[j], shf1);
    }
}

// Reserved for reference use only.
// static inline void gf16m_mul_4x4_16way_gf16m_1_u64_16_add_cross(const gf16m_t k, const __m256i* restrict bt,
//                                                                 __m256i* restrict c_lh_256) {
//     for (int i = 0; i < 2; i++) {
//         for (int j = 0; j < 4; j++) {
//             __m256i kkk = mtk2_16[(k[i * 8 + j + 4] << 4) ^ k[i * 8 + j]];
//             c_lh_256[i * 2] ^= _mm256_shuffle_epi8(kkk, bt[j * 2]);
//             c_lh_256[i * 2 + 1] ^= _mm256_shuffle_epi8(kkk, bt[j * 2 + 1]);
//         }
//     }
// }

/**
 * Parallel Operation Transposed Matrix
 */
static inline void gf16m_transpose_u64_4way(__m256i* ap, const __m256i a) {
    __m256i shf0 = _mm256_setr_epi8(0, 2, 1, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12, 14, 13, 15,  //
                                    0, 2, 1, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12, 14, 13, 15);

    __m256i shf1 = _mm256_setr_epi8(0, 4, 1, 5, 2, 6, 3, 7, 8, 12, 9, 13, 10, 14, 11, 15,  //
                                    0, 4, 1, 5, 2, 6, 3, 7, 8, 12, 9, 13, 10, 14, 11, 15);

    __m256i mask_0 = _mm256_set1_epi64x(0xf00ff00ff00ff00full);
    __m256i mask_1 = _mm256_set1_epi64x(0x0f000f000f000f00ull);

    __m256i a_0 = _mm256_shuffle_epi8(a, shf0);
    __m256i temp = (a_0 & mask_0) ^ (_mm256_slli_epi64(a_0, 4) & mask_1) ^ _mm256_srli_epi64(a_0 & mask_1, 4);
    _mm256_store_si256(ap, _mm256_shuffle_epi8(temp, shf1));
}

/**
 * To achieve four sets of parallel computations, we changed the data structure of ARRAY from X[lsq_SNOVA][n_SNOVA] to
 * X[n_SNOVA][lsq_SNOVA]. And we changed the order of the for loop for lsq_SNOVA.
 */
void evaluation_4x4_avx2_vtl(map_group1_u64* restrict map1_u64, uint64_t* restrict sig_var_u64, uint64_t* restrict P22,
                             uint64_t* restrict hash_in_GF16Matrix_u64) {
    __m256i Left_256[n_SNOVA][l_SNOVA] __attribute__((aligned(32)));
    __m256i Right_256_lh[n_SNOVA][l_SNOVA][2] __attribute__((aligned(32)));  // 0: lo 4bits, 1: hi 4bits

    __m256i* Aalpha_256 = (__m256i*)map1_u64->Aalpha;
    __m256i* Balpha_256 = (__m256i*)map1_u64->Balpha;
    __m256i* Qalpha1_256 = (__m256i*)map1_u64->Qalpha1;
    __m256i* Qalpha2_256 = (__m256i*)map1_u64->Qalpha2;

    __m256i l_mask = _mm256_set1_epi64x(0x0f0f0f0f0f0f0f0full);
    __m256i h_mask = _mm256_set1_epi64x(0xf0f0f0f0f0f0f0f0ull);

    // Parallel Operation Transposed Matrix
    uint64_t sig_var_t_u64[n_SNOVA_mult4] __attribute__((aligned(32))) = {0};
    __m256i* sig_var_256 = (__m256i*)sig_var_u64;
    __m256i* sig_var_t_256 = (__m256i*)sig_var_t_u64;
    for (int index = 0; index < n_SNOVA_mult4 / 4; ++index) {
        gf16m_transpose_u64_4way(sig_var_t_256 + index, sig_var_256[index]);
    }

    for (int index = 0; index < n_SNOVA; ++index) {
        __m256i Right_temp[4];
        __m256i temp_16x4[4];
        for (int alpha = 0; alpha < l_SNOVA; alpha++) {
            __m256i t256;
            gf16m_256_mul_4x4_4way_14(sig_var_t_u64[index], Qalpha1_256[alpha], &t256);
            gf16m_256_mul_4x4_4way_44(Aalpha_256[alpha], t256, Left_256[index] + alpha);

            gf16m_256_mul_4x4_4way_14(sig_var_u64[index], Balpha_256[alpha], &t256);
            gf16m_256_mul_4x4_4way_44(Qalpha2_256[alpha], t256, Right_temp + alpha);
        }

        convert_4x4x4_to_16x4(temp_16x4, Right_temp);
        for (int j = 0; j < 4; j++) {
            Right_256_lh[index][j][0] = temp_16x4[j] & l_mask;
            Right_256_lh[index][j][1] = _mm256_srli_epi16(temp_16x4[j] & h_mask, 4);
        }
    }

    uint64_t P_u64_temp[n_SNOVA_mult4] __attribute__((aligned(32))) = {0};
    uint8_t Public[n_SNOVA_mult4][sq_rank / 2] __attribute__((aligned(32)));

    __m256i tl[l_SNOVA] __attribute__((aligned(32)));
    __m256i tl0[l_SNOVA] __attribute__((aligned(32)));

    __m256i shf_pub = _mm256_setr_epi8(0, 2, 1, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12, 14, 13, 15,  //
                                       0, 2, 1, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12, 14, 13, 15);
    __m256i mask_0 = _mm256_set1_epi64x(0xf00ff00ff00ff00full);
    __m256i mask_1 = _mm256_set1_epi64x(0x0f000f000f000f00ull);

    for (int mi = 0; mi < m_SNOVA; ++mi) {
        __m256i sum = _mm256_setzero_si256();

        for (int nj = 0; nj < v_SNOVA; ++nj) {
            for (int nk = 0; nk < v_SNOVA; ++nk) {
                P_u64_temp[nk] = map1_u64->P11[mi][nj][nk];
            }

            for (int nk = v_SNOVA; nk < n_SNOVA; ++nk) {
                P_u64_temp[nk] = map1_u64->P12[mi][nj][nk - v_SNOVA];
            }
            for (int nk = 0; nk < n_SNOVA; nk += 4) {
                __m256i pub256 = _mm256_setr_epi64x(P_u64_temp[nk], P_u64_temp[nk + 1], P_u64_temp[nk + 2], P_u64_temp[nk + 3]);
                __m256i pub_0 = _mm256_shuffle_epi8(pub256, shf_pub);
                _mm256_store_si256(
                    (__m256i*)Public[nk],  //
                    (pub_0 & mask_0) ^ (_mm256_slli_epi64(pub_0, 4) & mask_1) ^ _mm256_srli_epi64(pub_0 & mask_1, 4));
            }

            __m256i tl0_cross[l_SNOVA] __attribute__((aligned(32))) = {0};

            for (int nk = 0; nk < n_SNOVA; ++nk) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 4; j++) {
                        __m256i k_lh = mtk2_16[Public[nk][i * 4 + j]];
                        tl0_cross[i * 2] ^= _mm256_shuffle_epi8(k_lh, Right_256_lh[nk][j][0]);
                        tl0_cross[i * 2 + 1] ^= _mm256_shuffle_epi8(k_lh, Right_256_lh[nk][j][1]);
                    }
                }
            }

            tl0[0] = (tl0_cross[0] & l_mask) ^ (_mm256_slli_epi16(tl0_cross[1], 4) & h_mask);
            tl0[1] = (_mm256_srli_epi16(tl0_cross[0], 4) & l_mask) ^ (tl0_cross[1] & h_mask);
            tl0[2] = (tl0_cross[2] & l_mask) ^ (_mm256_slli_epi16(tl0_cross[3], 4) & h_mask);
            tl0[3] = (_mm256_srli_epi16(tl0_cross[2], 4) & l_mask) ^ (tl0_cross[3] & h_mask);

            convert_16x4_to_4x4x4(tl, tl0);
            for (int alpha = 0; alpha < l_SNOVA; alpha++) {
                gf16m_256_mul_4x4_4way_44_add(Left_256[nj][alpha], tl[alpha], &sum);
            }
        }

        for (int nj = v_SNOVA; nj < n_SNOVA; ++nj) {
            for (int nk = 0; nk < v_SNOVA; ++nk) {
                P_u64_temp[nk] = map1_u64->P21[mi][nj - v_SNOVA][nk];
            }

            for (int nk = v_SNOVA; nk < n_SNOVA; ++nk) {
                P_u64_temp[nk] = (P22 + (mi * o_SNOVA + nj - v_SNOVA) * o_SNOVA)[nk - v_SNOVA];
            }
            for (int nk = 0; nk < n_SNOVA; nk += 4) {
                __m256i pub256 = _mm256_setr_epi64x(P_u64_temp[nk], P_u64_temp[nk + 1], P_u64_temp[nk + 2], P_u64_temp[nk + 3]);
                __m256i pub_0 = _mm256_shuffle_epi8(pub256, shf_pub);
                _mm256_store_si256(
                    (__m256i*)Public[nk],  //
                    (pub_0 & mask_0) ^ (_mm256_slli_epi64(pub_0, 4) & mask_1) ^ _mm256_srli_epi64(pub_0 & mask_1, 4));
            }

            __m256i tl0_cross[l_SNOVA] __attribute__((aligned(32))) = {0};

            for (int nk = 0; nk < n_SNOVA; ++nk) {
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 4; j++) {
                        __m256i k_lh = mtk2_16[Public[nk][i * 4 + j]];
                        tl0_cross[i * 2] ^= _mm256_shuffle_epi8(k_lh, Right_256_lh[nk][j][0]);
                        tl0_cross[i * 2 + 1] ^= _mm256_shuffle_epi8(k_lh, Right_256_lh[nk][j][1]);
                    }
                }
            }

            tl0[0] = (tl0_cross[0] & l_mask) ^ (_mm256_slli_epi16(tl0_cross[1], 4) & h_mask);
            tl0[1] = (_mm256_srli_epi16(tl0_cross[0], 4) & l_mask) ^ (tl0_cross[1] & h_mask);
            tl0[2] = (tl0_cross[2] & l_mask) ^ (_mm256_slli_epi16(tl0_cross[3], 4) & h_mask);
            tl0[3] = (_mm256_srli_epi16(tl0_cross[2], 4) & l_mask) ^ (tl0_cross[3] & h_mask);

            convert_16x4_to_4x4x4(tl, tl0);
            for (int alpha = 0; alpha < l_SNOVA; alpha++) {
                gf16m_256_mul_4x4_4way_44_add(Left_256[nj][alpha], tl[alpha], &sum);
            }
        }
        hash_in_GF16Matrix_u64[mi] = _mm256_extract_epi64(sum, 0) ^ _mm256_extract_epi64(sum, 1) ^  //
                                     _mm256_extract_epi64(sum, 2) ^ _mm256_extract_epi64(sum, 3);   //
    }
}

int verify_signture_4x4(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature, const uint8_t* pk) {
    uint8_t signed_hash[bytes_hash];
    const uint8_t* pt_salt = pt_signature + bytes_signature;
    public_key* pk_stru = (public_key*)pk;

    // ------- u64 new -------
    map_group1_u64 map1_u64;
    uint64_t hash_in_GF16Matrix_u64[m_SNOVA] __attribute__((aligned(32)));
    uint64_t* signed_hash_u64 = (uint64_t*)signed_hash;
    uint64_t sig_var_u64[n_SNOVA_mult4] __attribute__((aligned(32)));

    createSignedHash(pt_digest, bytes_digest, pk_stru->pt_public_key_seed, pt_salt, signed_hash);

    memcpy(sig_var_u64, pt_signature, bytes_signature);  // without salt

#if (o_SNOVA * l_SNOVA) & 0x1 == 1
    signed_hash[bytes_hash - 1] &= 0x0f;
#endif

    // generate PRNG part of public key
    gen_A_B_Q_P_4x4(&map1_u64, pk_stru->pt_public_key_seed);

    // evaluate signature GF16Matrix_u64 array. Using AVX2 for parallel computation with 16 sets.
    evaluation_4x4_avx2_vtl(&map1_u64, sig_var_u64, (uint64_t*)(pk_stru->P22), hash_in_GF16Matrix_u64);

    int result = 0;
    for (int i = 0; i < bytes_hash / 8; ++i) {
        if (hash_in_GF16Matrix_u64[i] != signed_hash_u64[i]) {
            result = -1;
            break;
        }
    }
    return result;
}

#endif