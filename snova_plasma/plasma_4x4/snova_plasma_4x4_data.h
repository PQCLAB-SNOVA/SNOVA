/*
Use uint64_t to store a GF16 4x4 matrix (with each half-byte storing one element of GF16)."
*/

#ifndef PLASMA_4x4_DATA_H
#define PLASMA_4x4_DATA_H

#include "../../deriv_params.h"
#include "../../gf16.h"
#include "../../snova_kernel.h"

#define v_SNOVA_mult4 ((v_SNOVA + 3) / 4 * 4)
#define n_SNOVA_mult4 ((n_SNOVA + 3) / 4 * 4)

#define GAUSS_ROW (m_SNOVA * lsq_SNOVA + 1)
#define GAUSS_ROW32 ((GAUSS_ROW + 31) / 32)
#define GAUSS_ROW_mult32 (GAUSS_ROW32 * 32)
#define GAUSS_COL (m_SNOVA * lsq_SNOVA)
#define GAUSS_COL_mult32 ((GAUSS_COL + 31) / 32 * 32)

#define mvl_SNOVA32 ((m_SNOVA * v_SNOVA * l_SNOVA + 31) / 32)
#define mvl_SNOVA (mvl_SNOVA32 * 32)

#define ml_SNOVA32 ((m_SNOVA * l_SNOVA + 31) / 32)
#define ml_SNOVA (ml_SNOVA32 * 32)

#define vl_SNOVA (v_SNOVA * l_SNOVA)
#define lcube_SNOVA (((lsq_SNOVA * l_SNOVA + 31) / 32) * 32)

#define vl4_SNOVA32 ((vl_SNOVA * lcube_SNOVA + 31) / 32)


#define map_group1_u64_byte                                                                                            \
    (m_SNOVA * v_SNOVA * v_SNOVA + m_SNOVA * v_SNOVA * o_SNOVA + m_SNOVA * o_SNOVA * v_SNOVA + lsq_SNOVA + lsq_SNOVA + \
     lsq_SNOVA + lsq_SNOVA) *                                                                                          \
        8

#define POD_u64(entry, a, b, c, d, e, f, g, h, i, j)                                                                \
    gf16_get_mul(get_gf16m_u64(entry, a, b), gf16_get_mul(get_gf16m_u64(entry, c, d), get_gf16m_u64(entry, e, f)) ^ \
                                                 gf16_get_mul(get_gf16m_u64(entry, g, h), get_gf16m_u64(entry, i, j)))

typedef uint64_t P11_u64_t[m_SNOVA][v_SNOVA][v_SNOVA] __attribute__((aligned(32)));
typedef uint64_t P12_u64_t[m_SNOVA][v_SNOVA][o_SNOVA] __attribute__((aligned(32)));
typedef uint64_t P21_u64_t[m_SNOVA][o_SNOVA][v_SNOVA] __attribute__((aligned(32)));
typedef uint64_t Aalpha_u64_t[lsq_SNOVA] __attribute__((aligned(32)));
typedef uint64_t Balpha_u64_t[lsq_SNOVA] __attribute__((aligned(32)));
typedef uint64_t Qalpha1_u64_t[lsq_SNOVA] __attribute__((aligned(32)));
typedef uint64_t Qalpha2_u64_t[lsq_SNOVA] __attribute__((aligned(32)));

typedef struct {
    P11_u64_t P11;
    P12_u64_t P12;
    P21_u64_t P21;
    Aalpha_u64_t Aalpha;
    Balpha_u64_t Balpha;
    Qalpha1_u64_t Qalpha1;
    Qalpha2_u64_t Qalpha2;
} map_group1_u64;

static uint64_t S_u64[l_SNOVA] = {0};

typedef uint64_t T12_u64_t[v_SNOVA][o_SNOVA] __attribute__((aligned(32)));
typedef uint64_t F11_u64_t[m_SNOVA][v_SNOVA][v_SNOVA] __attribute__((aligned(32)));
typedef uint64_t F12_u64_t[m_SNOVA][v_SNOVA][o_SNOVA] __attribute__((aligned(32)));
typedef uint64_t F21_u64_t[m_SNOVA][o_SNOVA][v_SNOVA] __attribute__((aligned(32)));

typedef struct {
    F11_u64_t F11;
    F12_u64_t F12;
    F21_u64_t F21;
} map_group2_u64;

uint32_t mul_count = 0;

static inline uint64_t gf16m_scale_4x4(uint64_t a, gf16_t k) {
    uint64_t mask = 0;
    uint64_t t = 0;
    for (int j = 0; j < 4; j++) {
        mask = (k & 0x0000000000000001ull) * 0xffffffffffffffffull;
        t ^= (a & mask);
        mask = ((a >> 3) & 0x1111111111111111ull);
        a = ((a ^ (mask * 0x9)) << 1) ^ mask;
        k >>= 1;
    }
    return t;
}

static inline uint8_t get_gf16m_u64(uint64_t a, uint8_t x, uint8_t y) {
    uint8_t i = (x << 2) ^ y;
    return (uint8_t)((a >> (i << 2)) & 0xfull);
}

static inline void set_gf16m_u64(uint64_t* a, uint8_t x, uint8_t y, uint64_t value) {
    uint8_t i = (x << 2) ^ y;
    uint64_t mask = ~(0xfull << (i << 2));
    *a = (*a & mask) ^ (value << (i << 2));
}

static inline void set_gf16m_u64_init(uint64_t* a, uint8_t x, uint8_t y, uint64_t value) {
    uint8_t i = (x << 2) ^ y;
    *a ^= (value << (i << 2));
}

static inline uint64_t gf16m_transpose_u64(uint64_t a) {
    uint64_t ap = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            set_gf16m_u64_init(&ap, i, j, (uint64_t)get_gf16m_u64(a, j, i));
        }
    }
    return ap;
}

static inline uint64_t be_aS_u64(uint8_t a) {
    uint64_t target = 0;
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m_u64_init(&target, i, j, (uint64_t)gf16_get_mul((8 - (i + j)), a));
        }
    }
    return target;
}

static inline gf16_t gf16m_det_u64(uint64_t entry) {
    gf16_t d0 = gf16_get_mul(get_gf16m_u64(entry, 0, 0), POD_u64(entry, 1, 1, 2, 2, 3, 3, 2, 3, 3, 2) ^
                                                             POD_u64(entry, 1, 2, 2, 1, 3, 3, 2, 3, 3, 1) ^
                                                             POD_u64(entry, 1, 3, 2, 1, 3, 2, 2, 2, 3, 1));

    gf16_t d1 = gf16_get_mul(get_gf16m_u64(entry, 0, 1), POD_u64(entry, 1, 0, 2, 2, 3, 3, 2, 3, 3, 2) ^
                                                             POD_u64(entry, 1, 2, 2, 0, 3, 3, 2, 3, 3, 0) ^
                                                             POD_u64(entry, 1, 3, 2, 0, 3, 2, 2, 2, 3, 0));

    gf16_t d2 = gf16_get_mul(get_gf16m_u64(entry, 0, 2), POD_u64(entry, 1, 0, 2, 1, 3, 3, 2, 3, 3, 1) ^
                                                             POD_u64(entry, 1, 1, 2, 0, 3, 3, 2, 3, 3, 0) ^
                                                             POD_u64(entry, 1, 3, 2, 0, 3, 1, 2, 1, 3, 0));

    gf16_t d3 = gf16_get_mul(get_gf16m_u64(entry, 0, 3), POD_u64(entry, 1, 0, 2, 1, 3, 2, 2, 2, 3, 1) ^
                                                             POD_u64(entry, 1, 1, 2, 0, 3, 2, 2, 2, 3, 0) ^
                                                             POD_u64(entry, 1, 2, 2, 0, 3, 1, 2, 1, 3, 0));
    return d0 ^ d1 ^ d2 ^ d3;
}

static inline void be_invertible_by_add_aS_u64(uint64_t* source) {
    if (gf16m_det_u64(*source) == 0) {
        for (uint8_t a = 1; a < 16; ++a) {
            *source ^= be_aS_u64(a);
            if (gf16m_det_u64(*source) != 0) {
                return;
            }
        }
    }
}

static inline uint64_t be_aI_u64(uint8_t a) {
    uint64_t target = 0;
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m_u64_init(&target, i, j, (i == j) ? a : 0);
        }
    }
    return target;
}

void gen_a_FqS_u64(uint8_t* c, uint64_t* pt_matrix) {
    uint8_t c_byte[4] = {c[0] & 0xf, (c[0] >> 4), c[1] & 0xf, (c[1] >> 4)};
    *pt_matrix = be_aI_u64(c_byte[0]);
    for (int i = 1; i < rank - 1; ++i) {
        *pt_matrix ^= gf16m_scale_4x4(S_u64[i], c_byte[i]);
    }
    uint8_t zero = ct_gf16_is_not_zero(c_byte[rank - 1]);
    uint8_t val = zero * c_byte[rank - 1] + (1 - zero) * (15 + ct_gf16_is_not_zero(c_byte[0]) - c_byte[0]);
    *pt_matrix ^= gf16m_scale_4x4(S_u64[rank - 1], val);
}

/**
 * Generate the random part of public key
 * @param map - P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param pt_public_key_seed - input
 */
void gen_A_B_Q_P_4x4(map_group1_u64* restrict map, uint8_t* restrict pt_public_key_seed) {
    uint8_t temp[lsq_SNOVA * l_SNOVA];

    pk_expand(pt_public_key_seed, (uint8_t*)map);
    memcpy(temp, ((uint8_t*)map->Qalpha1), lsq_SNOVA * l_SNOVA / 2);
    // printf("Aalpha: \n");
    for (int index = 0; index < lsq_SNOVA; ++index) {
        be_invertible_by_add_aS_u64(map->Aalpha + index);
    }

    // cout << "B======" << endl;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        be_invertible_by_add_aS_u64(map->Balpha + index);
    }
    uint8_t* pt_array = (uint8_t*)(map->Qalpha1) + l_SNOVA * lsq_SNOVA / 2;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        gen_a_FqS_u64(pt_array, map->Qalpha2 + index);
        // be_invertible_by_add_aS_u64(map->Qalpha2 + index);
        pt_array += (l_SNOVA / 2);
    }

    pt_array = temp;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        gen_a_FqS_u64(pt_array, map->Qalpha1 + index);
        // be_invertible_by_add_aS_u64(map->Qalpha1 + index);
        pt_array += (l_SNOVA / 2);
    }
}

#endif