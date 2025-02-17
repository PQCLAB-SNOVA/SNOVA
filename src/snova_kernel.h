#ifndef SNOVA_KNL_H
#define SNOVA_KNL_H

#include "gf16_matrix_inline.h"
#include "snova.h"
#include "ct_functions.h"
#include "aes/snova_aes.h"

static gf16m_t S[l_SNOVA] = {0};
static uint32_t xS[l_SNOVA][lsq_SNOVA] = {0};
static int S_is_init = 0;
// GF[x]/(x^4+x+1) reduction
static inline uint32_t gf16_reduce(uint32_t idx)
{
    uint32_t res, upper;

    res = idx & 0x49249249; // Octal 0o11111111111
    upper = idx >> 12;
    res = res ^ upper ^ (upper << 3);
    upper = res >> 12;
    res = res ^ upper ^ (upper << 3);
    upper = res >> 12;
    res = res ^ upper ^ (upper << 3);

    return res & 0x249;
}

// Conversion 4 bit -> 32 bit representation
static inline uint32_t gf16_from_nibble(uint8_t idx)
{
    uint32_t middle = idx | idx << 4;
    return (middle & 0x41) | ((middle << 2) & 0x208);
}

// Conversion 32 bit -> 4 bit representation
static inline uint8_t gf16_to_nibble(uint32_t idx)
{
    uint32_t res = gf16_reduce(idx);
    res = res | (res >> 4);
    return (res & 0x5) | ((res >> 2) & 0xa);
}

// Conversion 32 bit -> 4 bit representation
static inline uint8_t xgf16_to_nibble(uint32_t res)
{
    res = res | (res >> 4);
    return (res & 0x5) | ((res >> 2) & 0xa);
}

// Constant time GF16 inverse
// x^16 - x = 0 implies x^14 = x^-1
static inline uint32_t gf16_inv(uint32_t val)
{
    val = gf16_reduce(val);
    uint32_t res2 = gf16_reduce(val * val);
    uint32_t res4 = gf16_reduce(res2 * res2);
    uint32_t res8 = gf16_reduce(res4 * res4);

    return gf16_reduce(res2 * ((res4 * res8) & 0x49249249));
}


/**
 * Generate elements of F16[S]
 */
void gen_S_array(void) {
    if (S_is_init) {
        return;
    }

    S_is_init = 1;
    be_aI(S[0], 1);
    be_the_S(S[1]);
    for (int index = 2; index < l_SNOVA; index++) {
        gf16m_mul(S[index - 1], S[1], S[index]);
    }

    for (int index = 0; index < l_SNOVA; index++)
        for (int ij = 0; ij < lsq_SNOVA; ij++)
            xS[index][ij] = gf16_from_nibble(S[index][ij]);
}

/**
 * shake256
 * @param pt_seed_array - Pointer to the hash input.
 * @param input_bytes - hash lenth.
 * @param pt_output_array - Pointer to the hash output.
 * @param output_bytes - hash input.
 */
void shake256(const uint8_t* pt_seed_array, int input_bytes, uint8_t* pt_output_array, int output_bytes) {
    Keccak_HashInstance hashInstance;
    Keccak_HashInitialize_SHAKE256(&hashInstance);
    Keccak_HashUpdate(&hashInstance, pt_seed_array, 8 * input_bytes);
    Keccak_HashFinal(&hashInstance, NULL);
    Keccak_HashSqueeze(&hashInstance, pt_output_array, 8 * output_bytes);
}

/**
 * Convert one byte of data to GF16 representation (using only half of the
 * byte). Example: <bytes 12 34 56 78 9a bc> -> <bytes 02 01 04 03 05 ..... 0c
 * 0b>
 * @param byte_array - input (bytes)
 * @param gf16_array - output (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_bytes_to_GF16s(const uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s) {
    int i;
    int pairs = num_of_GF16s >> 1;

    // Convert each byte into two GF16 values
    for (i = 0; i < pairs; ++i) {
        gf16_array[i * 2] = byte_array[i] & 0x0F;
        gf16_array[i * 2 + 1] = (byte_array[i] >> 4) & 0x0F;
    }

    // Handle the last GF16 value if num_of_GF16s is odd
    if (num_of_GF16s % 2 == 1) {
        gf16_array[num_of_GF16s - 1] = byte_array[pairs] & 0x0F;
    }
}

/**
 * Convert two GF16 values to one byte.
 * Example:
 *  <bytes 02 01 04 03 05 ..... 0c 0b> -> <bytes 12 34 56 78 9a bc>
 * @param byte_array - output (bytes)
 * @param gf16_array - input (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_GF16s_to_bytes(uint8_t* byte_array, const gf16_t* gf16_array, int num_of_GF16s) {
    int i;
    int pairs = num_of_GF16s >> 1;

    // Convert pairs of GF16 values into one byte
    for (i = 0; i < pairs; ++i) {
        byte_array[i] = gf16_array[i * 2] | (gf16_array[i * 2 + 1] << 4);
    }

    // Handle the last GF16 value if num_of_GF16s is odd
    if (num_of_GF16s % 2 == 1) {
        byte_array[pairs] = gf16_array[num_of_GF16s - 1];
    }
}

/**
 * pk expand from seed
 * 
 * Using AES-CTR encryption as a hash function
 * AES ciphertext padded with zeros.
 * The iv is also padded with zeros.
 * Using input value as the AES key.
 * The ciphertext obtained from AES encryption serves as the output of the hash
 * function.
 * @param pt_public_key_seed - Pointer to the hash input. (Fixed length of 16)
 * @param out_pk - Pointer to the hash output. (Fixed length of
 * bytes_prng_public)
 */
void pk_expand(const uint8_t* pt_public_key_seed, uint8_t* out_pk) {
#if PK_EXPAND_SHAKE
    uint64_t pk_bytes[(bytes_prng_public + 7) / 8] __attribute__((aligned(32)));
    snova_shake(pt_public_key_seed, 16, pk_bytes, 8 * ((bytes_prng_public + 7) / 8));
    memcpy(out_pk, pk_bytes, bytes_prng_public);
#else
    AES_128_CTR(out_pk, bytes_prng_public, pt_public_key_seed, 16);
#endif
}

/**
 * Convert one byte of data to GF16 representation (using only half of the
 * byte). cut_in_half Example: <bytes 12 34 56 78 9a bc> -> <bytes 02 04 06 08
 * 0a 0c 01 03 05 07 09 0b>
 * @param byte_array - input (bytes)
 * @param gf16_array - output (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_bytes_to_GF16s_cut_in_half(const uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s) {
    int half_GF16s = (num_of_GF16s + 1) >> 1;
    int i;

    // Extract the lower 4 bits of each byte to the first half of gf16_array
    for (i = 0; i < half_GF16s; ++i) {
        gf16_array[i] = byte_array[i] & 0x0F;
    }

    // Extract the upper 4 bits of each byte to the second half of gf16_array
    for (i = 0; i < (num_of_GF16s >> 1); ++i) {
        gf16_array[i + half_GF16s] = byte_array[i] >> 4;
    }
}

/**
 * Convert two GF16 values to one byte.
 * Example:
 *  <bytes 02 04 06 08 0a 0c 01 03 05 07 09 0b> -> <bytes 12 34 56 78 9a bc>
 * @param byte_array - output (bytes)
 * @param gf16_array - input (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_GF16s_to_bytes_merger_in_half(uint8_t* byte_array, gf16_t* gf16_array, int num_of_GF16s) {
    int half_GF16s = (num_of_GF16s + 1) >> 1;
    int i;

    // Combine pairs of GF16 values into one byte
    for (i = 0; i < (num_of_GF16s >> 1); ++i) {
        byte_array[i] = gf16_array[i] | (gf16_array[i + half_GF16s] << 4);
    }

    // If num_of_GF16s is odd, handle the last GF16 value separately
    if (num_of_GF16s & 1) {
        byte_array[i] = gf16_array[i];
    }
}

/**
 * @param c - output
 * @param pt_matrix - input
 */
void gen_a_FqS(gf16_t* c, gf16m_t pt_matrix) {
    gf16m_t temp;
    be_aI(pt_matrix, c[0]);
    for (int i = 1; i < rank - 1; ++i) {
        gf16m_scale(S[i], c[i], temp);
        gf16m_add(pt_matrix, temp, pt_matrix);
    }
    gf16m_scale(S[rank - 1], (c[rank - 1] != 0) ? c[rank - 1] : 16 - (c[0] + (c[0] == 0)), temp);
    gf16m_add(pt_matrix, temp, pt_matrix);
    SNOVA_CLEAR(temp);
}

// Constant time version of gen_a_FqS
void gen_a_FqS_ct(gf16_t* c, gf16m_t pt_matrix) {
    uint32_t xTemp[lsq_SNOVA] = {0};
    uint32_t cX = gf16_from_nibble(c[0]);

    for (int ij = 0; ij < l_SNOVA; ij++)
        xTemp[ij * l_SNOVA + ij] = cX;

    for (int i1 = 1; i1 < l_SNOVA - 1; i1++) {
        cX = gf16_from_nibble(c[i1]);
        for (int ij = 0; ij < lsq_SNOVA; ij++)
            xTemp[ij] ^= cX * xS[i1][ij];
    }

    uint8_t zero = ct_gf16_is_not_zero(c[rank - 1]);
    uint8_t val = zero * c[rank - 1] + (1 - zero) * (15 + ct_gf16_is_not_zero(c[0]) - c[0]);

    cX = gf16_from_nibble(val);
    for (int ij = 0; ij < lsq_SNOVA; ij++)
        xTemp[ij] ^= cX * xS[l_SNOVA - 1][ij];

    for (int ij = 0; ij < lsq_SNOVA; ij++)
        pt_matrix[ij] = gf16_to_nibble(xTemp[ij]);
    
    SNOVA_CLEAR(xTemp);
}

/**
 * Generate the linear map T12
 * @param T12 - output
 * @param seed - input
 */
void gen_seeds_and_T12(T12_t T12, const uint8_t* seed) {
    gf16_t* pt_array;
    uint8_t prng_output_private[bytes_prng_private];
    gf16_t GF16_prng_output_private[GF16s_prng_private];

    shake256(seed, seed_length_private, prng_output_private, bytes_prng_private);
    convert_bytes_to_GF16s(prng_output_private, GF16_prng_output_private, GF16s_prng_private);

    pt_array = GF16_prng_output_private;
    for (int j = 0; j < v_SNOVA; ++j) {
        for (int k = 0; k < o_SNOVA; ++k) {
            gen_a_FqS_ct(pt_array, T12[j][k]);
            pt_array += rank;
        }
    }

    // Clear Secret!
    SNOVA_CLEAR(prng_output_private);
    SNOVA_CLEAR(GF16_prng_output_private);
}

/**
 * Generate the random part of public key
 * @param map - P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param pt_public_key_seed - input
 */

void gen_A_B_Q_P(map_group1* map, const uint8_t* pt_public_key_seed) {
    uint8_t prng_output_public[bytes_prng_public];
    uint8_t Q_temp[(sizeof(Qalpha1_t) + sizeof(Qalpha2_t)) / l_SNOVA];
    // ----- pt temp -----
    pk_expand(pt_public_key_seed, prng_output_public);
#if FIXED_ABQ
    convert_bytes_to_GF16s(prng_output_public, (uint8_t*)map, GF16s_prng_public - sizeof(Q_temp));
    memcpy(map->Aalpha, fixed_abq, 4 * m_SNOVA * alpha_SNOVA * lsq_SNOVA);
#else
    convert_bytes_to_GF16s(prng_output_public, (uint8_t*)map, GF16s_prng_public - sizeof(Q_temp));
    convert_bytes_to_GF16s(prng_output_public + sizeof(prng_output_public) - ((sizeof(Q_temp) + 1) >> 1), Q_temp, sizeof(Q_temp));

    for (int pi = 0; pi < m_SNOVA; ++pi) {
        for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
            be_invertible_by_add_aS(map->Aalpha[pi][alpha]);
        }
    }
    for (int pi = 0; pi < m_SNOVA; ++pi) {
        for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
            be_invertible_by_add_aS(map->Balpha[pi][alpha]);
        }
    }

    gf16_t* pt_array = Q_temp;
    for (int pi = 0; pi < m_SNOVA; ++pi) {
        for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
            gen_a_FqS(pt_array, map->Qalpha1[pi][alpha]);
            pt_array += l_SNOVA;
        }
    }
    for (int pi = 0; pi < m_SNOVA; ++pi) {
        for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
            gen_a_FqS(pt_array, map->Qalpha2[pi][alpha]);
            pt_array += l_SNOVA;
        }
    }
#endif
}

/**
 * Generate private key (F part)
 * @param map2 - output: F11 F12 F21
 * @param map1 - input: P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param T12 - input
 */
void gen_F_ref(map_group2* map2, map_group1* map1, T12_t T12) {
    gf16m_t temp;
    memcpy(map2->F11, map1->P11, m_SNOVA * v_SNOVA * v_SNOVA * sq_rank);
    memcpy(map2->F12, map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * sq_rank);
    memcpy(map2->F21, map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * sq_rank);

    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < v_SNOVA; ++j) {
            for (int k = 0; k < o_SNOVA; ++k) {
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(map1->P11[i][j][index], T12[index][k], temp);
                    gf16m_add(map2->F12[i][j][k], temp, map2->F12[i][j][k]);
                }
            }
        }
    }

    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < o_SNOVA; ++j) {
            for (int k = 0; k < v_SNOVA; ++k) {
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(T12[index][j], map1->P11[i][index][k], temp);
                    gf16m_add(map2->F21[i][j][k], temp, map2->F21[i][j][k]);
                }
            }
        }
    }

    // Clear Secret!
    SNOVA_CLEAR(temp);
}

/**
 * Generate public key (P22 part)
 * @param outP22 - output
 * @param T12 - input
 * @param P21 - input
 * @param F12 - input
 */
void gen_P22_ref(P22_byte_t outP22, T12_t T12, P21_t P21, F12_t F12) {
    gf16m_t temp1, temp2;
    P22_t P22 = {0};
    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < o_SNOVA; ++j) {
            for (int k = 0; k < o_SNOVA; ++k) {
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(T12[index][j], F12[i][index][k], temp1);
                    gf16m_mul(P21[i][j][index], T12[index][k], temp2);
                    gf16m_add(temp1, temp2, temp1);
                    gf16m_add(P22[i][j][k], temp1, P22[i][j][k]);
                    // P22[i][j][k].printout();
                }
            }
        }
    }

    convert_GF16s_to_bytes(outP22, (uint8_t*)P22, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);

    // Clear Secret!
    SNOVA_CLEAR(temp1);
    SNOVA_CLEAR(temp2);
}

/**
 * P22 byte to GF16
 * @param P22_gf16s - output
 * @param P22_bytes - input
 */
void input_P22(uint8_t* P22_gf16s, const uint8_t* P22_bytes) {
    convert_bytes_to_GF16s(P22_bytes, P22_gf16s, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);
}

/**
 * Pack expanded private key. esk = (key_elems, pt_private_key_seed).
 * @param esk - pointer to output expanded private key.
 * @param key_elems - pointer to input snova key elements.
 * @param pt_private_key_seed - pointer to input private key seed.
 */
void sk_pack(uint8_t* esk, snova_key_elems* key_elems, const uint8_t* pt_private_key_seed) {
    uint8_t* sk_gf16_ptr = (uint8_t*)(key_elems->map1.Aalpha);
    convert_GF16s_to_bytes_merger_in_half(esk, sk_gf16_ptr, (bytes_sk - (seed_length_public + seed_length_private)) * 2);
    memcpy(esk + (bytes_sk - (seed_length_public + seed_length_private)), key_elems->pk.pt_public_key_seed, seed_length_public);
    memcpy(esk + (bytes_sk - seed_length_private), pt_private_key_seed, seed_length_private);
}

/**
 * Unpack expanded secret key. skupk = (esk).
 * @param skupk - pointer to output private key (unpack).
 * @param esk - pointer to input expanded private key.
 */
void sk_unpack(sk_gf16* skupk, const uint8_t* esk) {
    convert_bytes_to_GF16s_cut_in_half(esk, (uint8_t*)skupk, (bytes_sk - (seed_length_public + seed_length_private)) * 2);
    memcpy(skupk->pt_public_key_seed, esk + (bytes_sk - (seed_length_public + seed_length_private)),
           seed_length_public + seed_length_private);
}

/**
 * Pack public key. pk = (key_elems).
 */
void pk_pack(uint8_t* pk, snova_key_elems* key_elems) {
    memcpy(pk, &key_elems->pk, bytes_pk);
}

/**
 * Unpack expend public key.
 */
void pkx_unpack(public_key_expand* pkx_unpck, public_key_expand_pack* pkx_pck) {
    convert_bytes_to_GF16s((uint8_t *)(pkx_pck) + seed_length_public, (uint8_t *)(pkx_unpck) + seed_length_public, sizeof(P22_t) + sizeof(map_group1));
    memcpy(pkx_unpck->pt_public_key_seed, pkx_pck->pt_public_key_seed, seed_length_public);
}

/**
 * Pack expend public key. 
 */
void pkx_pack(public_key_expand_pack* pkx_pck, public_key_expand* pkx_unpck) {
    convert_GF16s_to_bytes((uint8_t *)(pkx_pck) + seed_length_public, (uint8_t *)(pkx_unpck) + seed_length_public, sizeof(P22_t) + sizeof(map_group1));
    memcpy(pkx_pck->pt_public_key_seed, pkx_unpck->pt_public_key_seed, seed_length_public);
}

/**
 * expand public key
 * @param pkx - output
 * @param pk - input
 */
void expand_public_core(public_key_expand *pkx, const uint8_t *pk) {
    public_key* pk_stru = (public_key*)pk;
    memcpy(pkx->pt_public_key_seed, pk_stru->pt_public_key_seed, sizeof(pk_stru->pt_public_key_seed));
    // generate PRNG part of public key
    gen_A_B_Q_P(&(pkx->map1), pk_stru->pt_public_key_seed);
    // read  P22
    input_P22((uint8_t*)pkx->P22, (uint8_t*)pk_stru->P22);
}

/**
 * expand public key
 * @param pkx - output
 * @param pk - input
 */
void expand_public_pack_core(uint8_t* pkx_pck, const uint8_t *pk) {
    public_key_expand pkx_unpack;
    public_key* pk_stru = (public_key*)pk;
    memcpy(pkx_unpack.pt_public_key_seed, pk_stru->pt_public_key_seed, sizeof(pk_stru->pt_public_key_seed));
    // generate PRNG part of public key
    gen_A_B_Q_P(&(pkx_unpack.map1), pk_stru->pt_public_key_seed);
    // read  P22
    input_P22((uint8_t*)pkx_unpack.P22, (uint8_t*)pk_stru->P22);
    // pack gf16 -> bytes
    pkx_pack((public_key_expand_pack *)pkx_pck, &pkx_unpack);
}

/**
 * createHashOut
 */
void createSignedHash(const uint8_t* digest, uint64_t bytes_digest, const uint8_t* pt_public_key_seed,
                      const uint8_t* array_salt, uint8_t* signed_hash_out) {
    Keccak_HashInstance hashInstance;
    Keccak_HashInitialize_SHAKE256(&hashInstance);
    Keccak_HashUpdate(&hashInstance, pt_public_key_seed, 8 * seed_length_public);
    Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
    Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
    Keccak_HashFinal(&hashInstance, NULL);
    Keccak_HashSqueeze(&hashInstance, signed_hash_out, 8 * bytes_hash);
}

/**
 * evaluation
 * @param hash_in_GF16Matrix - pointer to output. (Evaluation Result)
 * @param pkx - pointer to expend pk.
 * @param signature_in_GF16Matrix - pointer to signature_in_GF16Matrix.
 */
void evaluation(
    gf16m_t* restrict hash_in_GF16Matrix, 
    const public_key_expand* restrict pkx, 
    gf16m_t* restrict signature_in_GF16Matrix
    ) {

    gf16m_t Left[m_SNOVA][alpha_SNOVA][n_SNOVA], Right[m_SNOVA][alpha_SNOVA][n_SNOVA];
    gf16m_t temp;

    for (int mi = 0; mi < m_SNOVA; ++mi) {
        // evaluate signature GF16Matrix array
        for (int si = 0; si < n_SNOVA; ++si) {
            gf16m_t signature_in_GF16Matrix_transpose;
            gf16m_transpose(signature_in_GF16Matrix[si], signature_in_GF16Matrix_transpose);
            for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {

                // Left[alpha][si]
                gf16m_mul(signature_in_GF16Matrix_transpose, pkx->map1.Qalpha1[mi][alpha], temp);
                gf16m_mul(pkx->map1.Aalpha[mi][alpha], temp, Left[mi][alpha][si]);
                // Right[alpha][si]
                gf16m_mul(pkx->map1.Qalpha2[mi][alpha], signature_in_GF16Matrix[si], temp);
                gf16m_mul(temp, pkx->map1.Balpha[mi][alpha], Right[mi][alpha][si]);
            }
        }
    }

    for (int mi = 0; mi < m_SNOVA; ++mi) {
        gf16m_set_zero(hash_in_GF16Matrix[mi]);
    }

    gf16m_t P[m_SNOVA][n_SNOVA][n_SNOVA] = {0};
    for (int mi = 0; mi < m_SNOVA; ++mi) {
        /*
                V        O
            +--------+--------+
            |        |        |
          V |  P11   |  P12   |
            |        |        |
            +--------+--------+   = P[n_SNOVA][n_SNOVA]
            |        |        |
          O |  P21   |  P22   |
            |        |        |
            +--------+--------+
        */
        for (int ni = 0; ni < v_SNOVA; ++ni) {
            for (int nj = 0; nj < v_SNOVA; ++nj) {
                gf16m_clone(P[mi][ni][nj], pkx->map1.P11[mi][ni][nj]);
            }

            for (int nj = v_SNOVA; nj < n_SNOVA; ++nj) {
                gf16m_clone(P[mi][ni][nj], pkx->map1.P12[mi][ni][nj - v_SNOVA]);
            }
        }
        for (int ni = v_SNOVA; ni < n_SNOVA; ++ni) {
            for (int nj = 0; nj < v_SNOVA; ++nj) {
                gf16m_clone(P[mi][ni][nj], pkx->map1.P21[mi][ni - v_SNOVA][nj]);
            }

            for (int nj = v_SNOVA; nj < n_SNOVA; ++nj) {
                gf16m_clone(P[mi][ni][nj], pkx->P22[mi][ni - v_SNOVA][nj - v_SNOVA]);
            }
        }
    }

    for (int mi = 0; mi < m_SNOVA; ++mi) {
        // main loop
        // hash_in_GF16Matrix[i] += L[alpha][ni] * P[mi][ni]][nj] * R[alpha][nj];
        for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
            int mi_prime = i_prime(mi, alpha);

            for (int ni = 0; ni < n_SNOVA; ++ni) {   
                gf16m_t sum_t0 = { 0 };
                for (int nj = 0; nj < n_SNOVA; ++nj) {
                    gf16m_mul(P[mi_prime][ni][nj], Right[mi][alpha][nj], temp);
                    gf16m_add(sum_t0, temp, sum_t0);
                }
                gf16m_mul(Left[mi][alpha][ni], sum_t0, temp);
                gf16m_add(hash_in_GF16Matrix[mi], temp, hash_in_GF16Matrix[mi]);
            }
        }
    }
}

/**
 * Verifies signature.
 * @param pt_digest - pointer to input digest.
 * @param pt_signature - pointer to output signature.
 * @param pk - pointer to output public key.
 * @returns - 0 if signature could be verified correctly and -1 otherwise
 */
int verify_signture_ref_core(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature, const public_key_expand* pkx) {
    uint8_t hash_in_bytes[bytes_hash];
    uint8_t signed_hash[bytes_hash];
    const uint8_t* pt_salt = pt_signature + bytes_signature;

    gf16m_t hash_in_GF16Matrix[m_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA];

    Keccak_HashInstance hashInstance;
    Keccak_HashInitialize_SHAKE256(&hashInstance);
    Keccak_HashUpdate(&hashInstance, pkx->pt_public_key_seed, 8 * seed_length_public);
    Keccak_HashUpdate(&hashInstance, pt_digest, 8 * bytes_digest);
    Keccak_HashUpdate(&hashInstance, pt_salt, 8 * bytes_salt);
    Keccak_HashFinal(&hashInstance, NULL);
    Keccak_HashSqueeze(&hashInstance, signed_hash, 8 * bytes_hash);

#if (o_SNOVA * l_SNOVA) & 0x1 == 1
    signed_hash[bytes_hash - 1] &= 0x0f;
#endif

    convert_bytes_to_GF16s(pt_signature, (gf16_t*)signature_in_GF16Matrix, GF16s_signature);
    evaluation(hash_in_GF16Matrix, pkx, signature_in_GF16Matrix);
    convert_GF16s_to_bytes(hash_in_bytes, (gf16_t*)hash_in_GF16Matrix, m_SNOVA * lsq_SNOVA);

    int result = 0;
    for (int i = 0; i < bytes_hash; ++i) {
        if (hash_in_bytes[i] != signed_hash[i]) {
            result = -1;
            break;
        }
    }

    return result;
}

/**
 * Regular entry point for verification
 */
int verify_signture_ref(const uint8_t *pt_digest, uint64_t bytes_digest, const uint8_t *pt_signature, const uint8_t *pk) {
    public_key_expand pkx;
    expand_public_core(&pkx, pk);
    return verify_signture_ref_core(pt_digest, bytes_digest, pt_signature, &pkx);
}

int verify_signture_pkx_ref(const uint8_t *pt_digest, uint64_t bytes_digest, const uint8_t *pt_signature, const uint8_t *pkx_pck) {
    public_key_expand pkx_unpck;
    pkx_unpack(&pkx_unpck, (public_key_expand_pack *)pkx_pck);
    return verify_signture_ref_core(pt_digest, bytes_digest, pt_signature, &pkx_unpck);
}

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
void sign_digest_core_ref(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, uint8_t* array_salt, Aalpha_t Aalpha,
                          Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, T12_t T12, F11_t F11, F12_t F12, F21_t F21,
                          const uint8_t pt_public_key_seed[seed_length_public], const uint8_t pt_private_key_seed[seed_length_private]) {
    gf16_t Gauss[m_SNOVA * lsq_SNOVA][m_SNOVA * lsq_SNOVA + 1];
    gf16_t Temp[lsq_SNOVA][lsq_SNOVA];
    gf16_t t_GF16, solution[m_SNOVA * lsq_SNOVA];

    gf16m_t Left_X_tmp, Right_X_tmp;
    gf16_t *Left_X, *Right_X;

    gf16m_t Left[m_SNOVA][alpha_SNOVA][v_SNOVA], Right[m_SNOVA][alpha_SNOVA][v_SNOVA];
    gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
    gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
    gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA];

    uint8_t signed_hash[bytes_hash];
    uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1];

    // temp
    gf16m_t gf16m_temp0;
    gf16m_t gf16m_temp1;
    gf16m_t gf16m_secret_temp0;

    Left_X = Left_X_tmp;
    Right_X = Right_X_tmp;
    int flag_redo = 1;
    uint8_t num_sign = 0;

    createSignedHash(digest, bytes_digest, pt_public_key_seed, array_salt, signed_hash);

    // put hash value in GF16 array
    convert_bytes_to_GF16s(signed_hash, hash_in_GF16, GF16s_hash);

    do {
        memset(Gauss, 0, sizeof(Gauss));
        num_sign++;
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

        // evaluate the vinegar part of central map
        for (int mi = 0; mi < m_SNOVA; ++mi) {
            for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_transpose(X_in_GF16Matrix[index], gf16m_temp0);
                    gf16m_mul(gf16m_temp0, Qalpha1[mi][alpha], gf16m_temp1);
                    gf16m_mul(Aalpha[mi][alpha], gf16m_temp1, Left[mi][alpha][index]);
                    gf16m_mul(Qalpha2[mi][alpha], X_in_GF16Matrix[index], gf16m_temp1);
                    gf16m_mul(gf16m_temp1, Balpha[mi][alpha], Right[mi][alpha][index]);
                }
            }
        }
        for (int mi = 0; mi < m_SNOVA; ++mi) {
            gf16m_set_zero(Fvv_in_GF16Matrix[mi]);
        }
        for (int mi = 0; mi < m_SNOVA; ++mi) {
            for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
                int mi_prime = i_prime(mi, alpha);
                for (int j = 0; j < v_SNOVA; ++j) {
                    for (int k = 0; k < v_SNOVA; ++k) {
                        gf16m_mul(Left[mi][alpha][j], F11[mi_prime][j][k], gf16m_temp0);
                        gf16m_mul(gf16m_temp0, Right[mi][alpha][k], gf16m_temp1);
                        gf16m_add(Fvv_in_GF16Matrix[mi], gf16m_temp1, Fvv_in_GF16Matrix[mi]);
                    }
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
            for (int index = 0; index < o_SNOVA; ++index) {
                for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
                    int mi_prime = i_prime(mi, alpha);
                    for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                        for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                            Temp[ti][tj] = 0;
                        }
                    }
                    for (int j = 0; j < v_SNOVA; ++j) {
                        gf16m_mul(Left[mi][alpha][j], F12[mi_prime][j][index], gf16m_temp0);
                        gf16m_mul(gf16m_temp0, Qalpha2[mi][alpha], Left_X_tmp);
                        Left_X = Left_X_tmp;
                        Right_X = Balpha[mi][alpha];
                        for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                            for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                                gf16_t temp3 = 0;
                                temp3 = gf16_get_mul(get_gf16m(Left_X, ti / rank, tj / rank),
                                                     get_gf16m(Right_X, tj % rank, ti % rank));
                                Temp[ti][tj] = gf16_get_add(Temp[ti][tj], temp3);
                            }
                        }
                    }

                    for (int j = 0; j < v_SNOVA; ++j) {
                        Left_X = Aalpha[mi][alpha];
                        gf16m_mul(Qalpha1[mi][alpha], F21[mi_prime][index][j], gf16m_temp0);
                        gf16m_mul(gf16m_temp0, Right[mi][alpha][j], Right_X_tmp);
                        Right_X = Right_X_tmp;
                        for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                            for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                                gf16_t temp2 = 0;
                                temp2 = gf16_get_mul(get_gf16m(Left_X, ti / rank, tj % rank),
                                                     get_gf16m(Right_X, tj / rank, ti % rank));
                                Temp[ti][tj] = gf16_get_add(Temp[ti][tj], temp2);
                            }
                        }
                    }
                    for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                        for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                            Gauss[mi * lsq_SNOVA + ti][index * lsq_SNOVA + tj] ^= Temp[ti][tj];
                        }
                    }
                }
            }
        }
        //
        // Gauss elimination
        for (int i = 0; i < m_SNOVA * lsq_SNOVA; ++i) {
            if (Gauss[i][i] == 0) {
                for (int j = i + 1; j < m_SNOVA * lsq_SNOVA; ++j) {
                    if (Gauss[j][i] != 0) {
                        for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
                            t_GF16 = Gauss[i][k];
                            Gauss[i][k] = Gauss[j][k];
                            Gauss[j][k] = t_GF16;
                        }
                        break;
                    }
                }
            }
            if (Gauss[i][i] == 0) {
                flag_redo = 1;
                break;
            }

            t_GF16 = inv(Gauss[i][i]);
            for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
                Gauss[i][k] = gf16_get_mul(Gauss[i][k], t_GF16);
            }

            for (int j = i + 1; j < m_SNOVA * lsq_SNOVA; ++j) {
                if (Gauss[j][i] != 0) {
                    t_GF16 = Gauss[j][i];
                    for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
                        Gauss[j][k] = gf16_get_add(Gauss[j][k], gf16_get_mul(Gauss[i][k], t_GF16));
                    }
                }
            }
        }

        if (!flag_redo) {
            for (int i = m_SNOVA * lsq_SNOVA - 1; i >= 0; --i) {
                t_GF16 = 0;
                for (int k = i + 1; k < m_SNOVA * lsq_SNOVA; ++k) {
                    t_GF16 = gf16_get_add(t_GF16, gf16_get_mul(Gauss[i][k], solution[k]));
                }
                solution[i] = gf16_get_add(Gauss[i][m_SNOVA * lsq_SNOVA], t_GF16);
            }
        }

    } while (flag_redo);
    // printf("times of Gauss elimination : %d\n", num_sign);
    for (int index = 0; index < o_SNOVA; ++index) {
        for (int i = 0; i < rank; ++i) {
            for (int j = 0; j < rank; ++j) {
                set_gf16m(X_in_GF16Matrix[index + v_SNOVA], i, j, solution[index * lsq_SNOVA + i * rank + j]);
            }
        }
    }

    for (int index = 0; index < v_SNOVA; ++index) {
        gf16m_clone(signature_in_GF16Matrix[index], X_in_GF16Matrix[index]);
        for (int i = 0; i < o_SNOVA; ++i) {
            gf16m_mul(T12[index][i], X_in_GF16Matrix[v_SNOVA + i], gf16m_secret_temp0);
            gf16m_add(signature_in_GF16Matrix[index], gf16m_secret_temp0, signature_in_GF16Matrix[index]);
        }
    }
    for (int index = 0; index < o_SNOVA; ++index) {
        gf16m_clone(signature_in_GF16Matrix[v_SNOVA + index], X_in_GF16Matrix[v_SNOVA + index]);
    }
    // output signature
    for (int index = 0; index < n_SNOVA * lsq_SNOVA; ++index) {
        ((gf16_t*)signature_in_GF16Matrix)[index] =
            get_gf16m(signature_in_GF16Matrix[index / lsq_SNOVA], (index % lsq_SNOVA) / l_SNOVA, (index % lsq_SNOVA) % l_SNOVA);
    }
    convert_GF16s_to_bytes(pt_signature, (gf16_t*)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);
    for (int i = 0; i < bytes_salt; ++i) {
        pt_signature[bytes_signature + i] = array_salt[i];
    }

    SNOVA_CLEAR(gf16m_secret_temp0);
}

#endif
