// SPDX-License-Identifier: MIT

/**
 * AVX2 Optimized implementation for q=16.
 *
 * Copyright (c) 2025 SNOVA TEAM
 */

#include <assert.h>
#include <stdalign.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "snova.h"

#if SNOVA_q != 16
#error "SNOVA_q != 16"
#endif

#include "symmetric.h"

#define o_SNOVA SNOVA_o
#define v_SNOVA SNOVA_v
#define l_SNOVA SNOVA_l

#define seed_length_public 16
#define seed_length_private 32
#define seed_length (seed_length_public + seed_length_private)

#define n_SNOVA (v_SNOVA + o_SNOVA)
#define m_SNOVA (o_SNOVA)
#define lsq_SNOVA (l_SNOVA * l_SNOVA)
#define alpha_SNOVA (l_SNOVA * l_SNOVA + l_SNOVA)

#define GF16s_hash (o_SNOVA * lsq_SNOVA)
#define GF16s_signature (n_SNOVA * lsq_SNOVA)
#define bytes_hash ((GF16s_hash + 1) >> 1)

#define rank (l_SNOVA)
#define sq_rank (rank * rank)  // matrix size

#define bytes_signature ((GF16s_signature + 1) >> 1)
#define bytes_salt 16
#define bytes_sig_with_salt (bytes_signature + bytes_salt)

#define GF16s_prng_public                                                                          \
    (sq_rank * (2 * (m_SNOVA * alpha_SNOVA) + m_SNOVA * (n_SNOVA * n_SNOVA - m_SNOVA * m_SNOVA)) + \
     rank * 2 * m_SNOVA * alpha_SNOVA)
#define bytes_prng_public ((GF16s_prng_public + 1) >> 1)

#define GF16s_prng_private (v_SNOVA * o_SNOVA * rank)
#define bytes_prng_private ((GF16s_prng_private + 1) >> 1)

#define bytes_pk (seed_length_public + ((m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA + 1) >> 1))
#define bytes_expend_pk (seed_length_public + ((m_SNOVA * (n_SNOVA * n_SNOVA + 4 * alpha_SNOVA) * sq_rank) + 1) / 2)

#define bytes_sk                                                                                                      \
    (((sq_rank * (4 * m_SNOVA * alpha_SNOVA + m_SNOVA * (v_SNOVA * v_SNOVA + v_SNOVA * o_SNOVA + o_SNOVA * v_SNOVA) + \
                  v_SNOVA * o_SNOVA) +                                                                                \
       1) >>                                                                                                          \
      1) +                                                                                                            \
     seed_length_public + seed_length_private)

#define mt(p, q) mt4b[((p) << 4) ^ (q)]
#define inv(gf16) inv4b[(gf16)]
#define gf16_get_add(a, b) ((a) ^ (b))
#define gf16_get_mul(a, b) (mt((a), (b)))

alignas(32) static uint8_t mt4b[256] = {0};
alignas(32) static uint8_t inv4b[16] = {0};

typedef uint8_t gf16_t;

/**
 * init gf16 tables
 */
static void init_gf16_tables(void) {
	static int gf16_tables_is_init = 0;
	if (gf16_tables_is_init) {
		return;
	}
	gf16_tables_is_init = 1;
	uint8_t F_star[15] = {1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9};  // Z2[x]/(x^4+x+1)
	for (int i = 0; i < 16; ++i) {
		mt(0, i) = mt(i, 0) = 0;
	}

	for (int i = 0; i < 15; ++i)
		for (int j = 0; j < 15; ++j) {
			mt(F_star[i], F_star[j]) = F_star[(i + j) % 15];
		}
	{
		int g = F_star[1], g_inv = F_star[14], gn = 1, gn_inv = 1;
		inv4b[0] = 0;
		inv4b[1] = 1;
		for (int index = 0; index < 14; index++) {
			inv4b[(gn = mt(gn, g))] = (gn_inv = mt(gn_inv, g_inv));
		}
	}
}

#define get_gf16m(gf16m, x, y) (gf16m[(((x) * rank) + (y))])
#define set_gf16m(gf16m, x, y, value) (gf16m[(((x) * rank) + (y))] = value)

typedef gf16_t gf16m_t[sq_rank];

// POD -> entry[a][b] * (entry[c][d] * entry[e][f] + entry[g][h] * entry[i][j])
#define POD(entry, a, b, c, d, e, f, g, h, i, j)                                                                    \
    gf16_get_mul(get_gf16m(entry, a, b), gf16_get_add(gf16_get_mul(get_gf16m(entry, c, d), get_gf16m(entry, e, f)), \
                                                      gf16_get_mul(get_gf16m(entry, g, h), get_gf16m(entry, i, j))))

/**
 * Zeroing the GF16 Matrix a.
 */
static inline void gf16m_set_zero(gf16m_t a) {
	memset(a, 0, sq_rank);
}

/**
 * Adding GF16 Matrices. c = a + b
 */
static inline void gf16m_add(const gf16m_t a, const gf16m_t b, gf16m_t c) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(c, i, j, gf16_get_add(get_gf16m(a, i, j), get_gf16m(b, i, j)));
		}
	}
}

/**
 * Multiplying GF16 Matrices. c = a * b
 */
static inline void gf16m_mul(const gf16m_t a, const gf16m_t b, gf16m_t c) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(c, i, j, gf16_get_mul(get_gf16m(a, i, 0), get_gf16m(b, 0, j)));
			for (int k = 1; k < rank; ++k) {
				set_gf16m(c, i, j, gf16_get_add(get_gf16m(c, i, j), gf16_get_mul(get_gf16m(a, i, k), get_gf16m(b, k, j))));
			}
		}
	}
}

/**
 * Scaling the GF16 Matrix. c = Scaling "a" by a factor of "k"
 */
static inline void gf16m_scale(const gf16m_t a, gf16_t k, gf16m_t c) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(c, i, j, gf16_get_mul(get_gf16m(a, i, j), k));
		}
	}
}

/**
 * Transposing the GF16 Matrix. ap = aT
 */
static inline void gf16m_transpose(const gf16m_t a, gf16m_t ap) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(ap, i, j, get_gf16m(a, j, i));
		}
	}
}

/**
 * Cloning the GF16 Matrix target = source
 */
static inline void gf16m_clone(gf16m_t target, const gf16m_t source) {
	memcpy(target, source, sq_rank);
}

/**
 * be_aI
 */
static inline void be_aI(gf16m_t target, gf16_t a) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(target, i, j, (i == j) ? a : 0);
		}
	}
}

/**
 * be_the_S
 */
static inline void be_the_S(gf16m_t target) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(target, i, j, (8 - (i + j)));
		}
	}
#if rank == 5
	set_gf16m(target, 4, 4, 9);
#endif
}

/**
 * Helper for rank5 gf16m_det5
 */
static inline gf16_t gf16m_det3(gf16m_t entry, int i0, int i1, int i2, int j0, int j1, int j2) {
	return gf16_get_add(
	           gf16_get_add(gf16_get_mul(get_gf16m(entry, j0, i0),
	                                     gf16_get_add(gf16_get_mul(get_gf16m(entry, j1, i1), get_gf16m(entry, j2, i2)),
	                                             gf16_get_mul(get_gf16m(entry, j1, i2), get_gf16m(entry, j2, i1)))),
	                        gf16_get_mul(get_gf16m(entry, j0, i1),
	                                     gf16_get_add(gf16_get_mul(get_gf16m(entry, j1, i0), get_gf16m(entry, j2, i2)),
	                                             gf16_get_mul(get_gf16m(entry, j1, i2), get_gf16m(entry, j2, i0))))),
	           gf16_get_mul(get_gf16m(entry, j0, i2), gf16_get_add(gf16_get_mul(get_gf16m(entry, j1, i0), get_gf16m(entry, j2, i1)),
	                        gf16_get_mul(get_gf16m(entry, j1, i1), get_gf16m(entry, j2, i0)))));
}

static inline gf16_t gf16m_det5(gf16m_t entry) {
	uint8_t a012 = gf16m_det3(entry, 0, 1, 2, 0, 1, 2);
	uint8_t b012 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 3), get_gf16m(entry, 4, 4)),
	                            gf16_get_mul(get_gf16m(entry, 3, 4), get_gf16m(entry, 4, 3)));

	uint8_t a013 = gf16m_det3(entry, 0, 1, 3, 0, 1, 2);
	uint8_t b013 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 2), get_gf16m(entry, 4, 4)),
	                            gf16_get_mul(get_gf16m(entry, 3, 4), get_gf16m(entry, 4, 2)));

	uint8_t a014 = gf16m_det3(entry, 0, 1, 4, 0, 1, 2);
	uint8_t b014 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 2), get_gf16m(entry, 4, 3)),
	                            gf16_get_mul(get_gf16m(entry, 3, 3), get_gf16m(entry, 4, 2)));

	uint8_t a023 = gf16m_det3(entry, 0, 2, 3, 0, 1, 2);
	uint8_t b023 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 1), get_gf16m(entry, 4, 4)),
	                            gf16_get_mul(get_gf16m(entry, 3, 4), get_gf16m(entry, 4, 1)));

	uint8_t a024 = gf16m_det3(entry, 0, 2, 4, 0, 1, 2);
	uint8_t b024 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 1), get_gf16m(entry, 4, 3)),
	                            gf16_get_mul(get_gf16m(entry, 3, 3), get_gf16m(entry, 4, 1)));

	uint8_t a034 = gf16m_det3(entry, 0, 3, 4, 0, 1, 2);
	uint8_t b034 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 1), get_gf16m(entry, 4, 2)),
	                            gf16_get_mul(get_gf16m(entry, 3, 2), get_gf16m(entry, 4, 1)));

	uint8_t a123 = gf16m_det3(entry, 1, 2, 3, 0, 1, 2);
	uint8_t b123 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 0), get_gf16m(entry, 4, 4)),
	                            gf16_get_mul(get_gf16m(entry, 3, 4), get_gf16m(entry, 4, 0)));

	uint8_t a124 = gf16m_det3(entry, 1, 2, 4, 0, 1, 2);
	uint8_t b124 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 0), get_gf16m(entry, 4, 3)),
	                            gf16_get_mul(get_gf16m(entry, 3, 3), get_gf16m(entry, 4, 0)));

	uint8_t a134 = gf16m_det3(entry, 1, 3, 4, 0, 1, 2);
	uint8_t b134 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 0), get_gf16m(entry, 4, 2)),
	                            gf16_get_mul(get_gf16m(entry, 3, 2), get_gf16m(entry, 4, 0)));

	uint8_t a234 = gf16m_det3(entry, 2, 3, 4, 0, 1, 2);
	uint8_t b234 = gf16_get_add(gf16_get_mul(get_gf16m(entry, 3, 0), get_gf16m(entry, 4, 1)),
	                            gf16_get_mul(get_gf16m(entry, 3, 1), get_gf16m(entry, 4, 0)));

	return gf16_get_mul(a012, b012) ^ gf16_get_mul(a013, b013) ^ gf16_get_mul(a014, b014) ^ gf16_get_mul(a023, b023) ^
	       gf16_get_mul(a024, b024) ^ gf16_get_mul(a034, b034) ^ gf16_get_mul(a123, b123) ^ gf16_get_mul(a124, b124) ^
	       gf16_get_mul(a134, b134) ^ gf16_get_mul(a234, b234);
}

/**
 * gf16m_det
 */
static inline gf16_t gf16m_det(gf16m_t entry) {
#if rank == 2
	return gf16_get_add(gf16_get_mul(get_gf16m(entry, 0, 0), get_gf16m(entry, 1, 1)),
	                    gf16_get_mul(get_gf16m(entry, 0, 1), get_gf16m(entry, 1, 0)));
	// (entry[0][0] * entry[1][1] + entry[0][1] * entry[1][0]);
#elif rank == 3
	return gf16_get_add(
	           gf16_get_add(
	               gf16_get_mul(get_gf16m(entry, 0, 0), gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 1), get_gf16m(entry, 2, 2)),
	                            gf16_get_mul(get_gf16m(entry, 1, 2), get_gf16m(entry, 2, 1)))),
	               // AAAAA(entry, 0, 0, 1, 1, 2, 2, 1, 2, 2, 1),
	               gf16_get_mul(get_gf16m(entry, 0, 1), gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2, 2)),
	                            gf16_get_mul(get_gf16m(entry, 1, 2), get_gf16m(entry, 2, 0))))),
	           gf16_get_mul(get_gf16m(entry, 0, 2), gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2, 1)),
	                        gf16_get_mul(get_gf16m(entry, 1, 1), get_gf16m(entry, 2, 0)))));

	/*
	(
	        (entry[0][0] * (entry[1][1] * entry[2][2] + entry[1][2] *
	entry[2][1])) + (entry[0][1] * (entry[1][0] * entry[2][2] + entry[1][2] *
	entry[2][0])) + (entry[0][2] * (entry[1][0] * entry[2][1] + entry[1][1] *
	entry[2][0]))
	)

	*/

	// gf16_get_mul(gf16_get_mul(get_gf16m(entry, 0, 1),
	// gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2,
	// 2)), gf16_get_mul(get_gf16m(entry, 1, 2), get_gf16m(entry, 2, 0)))))),
	// gf16_get_mul(gf16_get_mul(get_gf16m(entry, 0, 2),
	// gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2,
	// 1)), gf16_get_mul(get_gf16m(entry, 1, 1), get_gf16m(entry, 2, 0))))));

#elif rank == 4

	gf16_t d0 = gf16_get_mul(get_gf16m(entry, 0, 0), gf16_get_add(gf16_get_add(POD(entry, 1, 1, 2, 2, 3, 3, 2, 3, 3, 2),
	                         POD(entry, 1, 2, 2, 1, 3, 3, 2, 3, 3, 1)),
	                         POD(entry, 1, 3, 2, 1, 3, 2, 2, 2, 3, 1)));

	gf16_t d1 = gf16_get_mul(get_gf16m(entry, 0, 1), gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 2, 3, 3, 2, 3, 3, 2),
	                         POD(entry, 1, 2, 2, 0, 3, 3, 2, 3, 3, 0)),
	                         POD(entry, 1, 3, 2, 0, 3, 2, 2, 2, 3, 0)));

	gf16_t d2 = gf16_get_mul(get_gf16m(entry, 0, 2), gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 1, 3, 3, 2, 3, 3, 1),
	                         POD(entry, 1, 1, 2, 0, 3, 3, 2, 3, 3, 0)),
	                         POD(entry, 1, 3, 2, 0, 3, 1, 2, 1, 3, 0)));

	gf16_t d3 = gf16_get_mul(get_gf16m(entry, 0, 3), gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 1, 3, 2, 2, 2, 3, 1),
	                         POD(entry, 1, 1, 2, 0, 3, 2, 2, 2, 3, 0)),
	                         POD(entry, 1, 2, 2, 0, 3, 1, 2, 1, 3, 0)));

	return gf16_get_add(gf16_get_add(gf16_get_add(d0, d1), d2), d3);
#elif rank == 5
	return gf16m_det5(entry);
#else
#error "rank > 5 is not supported"
#endif
}

#define Keccak_HashInstance shake_t
#define Keccak_HashInitialize_SHAKE256 shake256_init
#define Keccak_HashUpdate(a, b, c) shake_absorb(a, b, (c) / 8)
#define Keccak_HashFinal(a, b) shake_finalize(a)
#define Keccak_HashSqueeze(a, b, c) shake_squeeze(b, (c) / 8, a)

typedef gf16m_t P11_t[m_SNOVA][v_SNOVA][v_SNOVA];
typedef gf16m_t P12_t[m_SNOVA][v_SNOVA][o_SNOVA];
typedef gf16m_t P21_t[m_SNOVA][o_SNOVA][v_SNOVA];
typedef gf16m_t Aalpha_t[m_SNOVA][alpha_SNOVA];
typedef gf16m_t Balpha_t[m_SNOVA][alpha_SNOVA];
typedef gf16m_t Qalpha1_t[m_SNOVA][alpha_SNOVA];
typedef gf16m_t Qalpha2_t[m_SNOVA][alpha_SNOVA];

typedef struct {
	P11_t P11;
	P12_t P12;
	P21_t P21;
	Aalpha_t Aalpha;
	Balpha_t Balpha;
	Qalpha1_t Qalpha1;
	Qalpha2_t Qalpha2;
} map_group1;

typedef gf16m_t T12_t[v_SNOVA][o_SNOVA];
typedef gf16m_t F11_t[m_SNOVA][v_SNOVA][v_SNOVA];
typedef gf16m_t F12_t[m_SNOVA][v_SNOVA][o_SNOVA];
typedef gf16m_t F21_t[m_SNOVA][o_SNOVA][v_SNOVA];

typedef struct {
	F11_t F11;
	F12_t F12;
	F21_t F21;
} map_group2;

typedef struct {
	Aalpha_t Aalpha;
	Balpha_t Balpha;
	Qalpha1_t Qalpha1;
	Qalpha2_t Qalpha2;
	T12_t T12;
	F11_t F11;
	F12_t F12;
	F21_t F21;
	uint8_t pt_public_key_seed[seed_length_public];
	uint8_t pt_private_key_seed[seed_length_private];
} sk_gf16;

typedef gf16m_t P22_t[m_SNOVA][o_SNOVA][o_SNOVA];
typedef uint8_t P22_byte_t[(m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA + 1) >> 1];  // byte

typedef struct {
	uint8_t pt_public_key_seed[seed_length_public];
	P22_byte_t P22;
} public_key;

typedef struct {
	uint8_t pt_public_key_seed[seed_length_public];
	P22_t P22;
	map_group1 map1;
} public_key_expand;

typedef struct {
	uint8_t pt_public_key_seed[seed_length_public];
	uint8_t P22_map1[((sizeof(P22_t) + sizeof(map_group1)))];
} public_key_expand_pack;

typedef struct {
	map_group1 map1;
	T12_t T12;
	map_group2 map2;
	public_key pk;
} snova_key_elems;

/**
 * be_aS
 */
static inline void be_aS(gf16m_t target, gf16_t a) {
	for (int i = 0; i < rank; ++i) {
		for (int j = 0; j < rank; ++j) {
			set_gf16m(target, i, j, gf16_get_mul((8 - (i + j)), a));
		}
	}
#if rank == 5
	set_gf16m(target, 4, 4, gf16_get_mul(9, a));
#endif
}

/**
 * be_invertible_by_add_aS
 */
static inline void be_invertible_by_add_aS(gf16m_t source) {
	gf16m_t temp;
	if (gf16m_det(source) == 0) {
		for (uint8_t a = 1; a < 16; ++a) {
			be_aS(temp, a);
			gf16m_add(temp, source, source);
			if (gf16m_det(source) != 0) {
				return;
			}
		}
	}
}

// #define i_prime(mi, alpha) ((alpha + mi) % o_SNOVA)
#define i_prime_inv(mi, alpha) ((o_SNOVA * alpha_SNOVA - alpha + mi) % o_SNOVA)

#if FIXED_ABQ
static uint8_t fixed_abq[4 * m_SNOVA * alpha_SNOVA * lsq_SNOVA] = {0};
#endif

static gf16m_t S[l_SNOVA] = {0};
static uint32_t xS[l_SNOVA][lsq_SNOVA] = {0};
static int S_is_init = 0;
// GF[x]/(x^4+x+1) reduction
static inline uint32_t gf16_reduce(uint32_t idx) {
	uint32_t res, upper;

	res = idx & 0x49249249;  // Octal 0o11111111111
	upper = idx >> 12;
	res = res ^ upper ^ (upper << 3);
	upper = res >> 12;
	res = res ^ upper ^ (upper << 3);
	upper = res >> 12;
	res = res ^ upper ^ (upper << 3);

	return res & 0x249;
}

// Conversion 4 bit -> 32 bit representation
static inline uint32_t gf16_from_nibble(uint8_t idx) {
	uint32_t middle = idx | idx << 4;
	return (middle & 0x41) | ((middle << 2) & 0x208);
}

// Conversion 32 bit -> 4 bit representation
static inline uint8_t gf16_to_nibble(uint32_t idx) {
	uint32_t res = gf16_reduce(idx);
	res = res | (res >> 4);
	return (res & 0x5) | ((res >> 2) & 0xa);
}

// Conversion 32 bit -> 4 bit representation
static inline uint8_t xgf16_to_nibble(uint32_t res) {
	res = res | (res >> 4);
	return (res & 0x5) | ((res >> 2) & 0xa);
}

// Constant time GF16 inverse
// x^16 - x = 0 implies x^14 = x^-1
static inline uint32_t gf16_inv(uint32_t val) {
	val = gf16_reduce(val);
	uint32_t res2 = gf16_reduce(val * val);
	uint32_t res4 = gf16_reduce(res2 * res2);
	uint32_t res8 = gf16_reduce(res4 * res4);

	return gf16_reduce(res2 * ((res4 * res8) & 0x49249249));
}

/**
 * Generate elements of F16[S]
 */
static void gen_S_array(void) {
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
		for (int ij = 0; ij < lsq_SNOVA; ij++) {
			xS[index][ij] = gf16_from_nibble(S[index][ij]);
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
static void pk_expand(const uint8_t* seed, uint8_t* pk_bytes) {
	snova_pk_expand(pk_bytes, NUM_GEN_PUB_BYTES, seed, SEED_LENGTH_PUBLIC);
}

// Constant time version of: (val != 0)
static uint32_t ct_gf16_is_not_zero(uint8_t val) {
	// return (val | (val >> 1) | (val >> 2) | (val >> 3)) & 1;
	return val != 0;
}

/**
 * Convert one byte of data to GF16 representation (using only half of the
 * byte). Example: <bytes 12 34 56 78 9a bc> -> <bytes 02 01 04 03 05 ..... 0c
 * 0b>
 * @param byte_array - input (bytes)
 * @param gf16_array - output (GF16)
 * @param num_of_GF16s - GF16 amount
 */
static int convert_bytes_to_GF16s(const uint8_t* byte_array, uint8_t* gf16_array, int num_of_GF16s) {
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
		return byte_array[pairs] & 0xF0;
	}

	return 0;
}

/**
 * Convert two GF16 values to one byte.
 * Example:
 *  <bytes 02 01 04 03 05 ..... 0c 0b> -> <bytes 12 34 56 78 9a bc>
 * @param byte_array - output (bytes)
 * @param gf16_array - input (GF16)
 * @param num_of_GF16s - GF16 amount
 */
static void convert_GF16s_to_bytes(uint8_t* byte_array, const uint8_t* gf16_array, int num_of_GF16s) {
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
 * Convert one byte of data to GF16 representation (using only half of the
 * byte). cut_in_half Example: <bytes 12 34 56 78 9a bc> -> <bytes 02 04 06 08
 * 0a 0c 01 03 05 07 09 0b>
 * @param byte_array - input (bytes)
 * @param gf16_array - output (GF16)
 * @param num_of_GF16s - GF16 amount
 */
static void convert_bytes_to_GF16s_cut_in_half(const uint8_t* byte_array, uint8_t* gf16_array, int num_of_GF16s) {
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
static void convert_GF16s_to_bytes_merger_in_half(uint8_t* byte_array, uint8_t* gf16_array, int num_of_GF16s) {
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
static void gen_a_FqS(gf16_t* c, gf16m_t pt_matrix) {
	gf16m_t temp;
	be_aI(pt_matrix, c[0]);
	for (int i = 1; i < rank - 1; ++i) {
		gf16m_scale(S[i], c[i], temp);
		gf16m_add(pt_matrix, temp, pt_matrix);
	}
	gf16m_scale(S[rank - 1], (c[rank - 1] != 0) ? c[rank - 1] : 16 - (c[0] + (c[0] == 0)), temp);
	gf16m_add(pt_matrix, temp, pt_matrix);
}

// Constant time version of gen_a_FqS
static void gen_a_FqS_ct(gf16_t* c, gf16m_t pt_matrix) {
	uint32_t xTemp[lsq_SNOVA] = {0};
	uint32_t cX = gf16_from_nibble(c[0]);

	for (int ij = 0; ij < l_SNOVA; ij++) {
		xTemp[ij * l_SNOVA + ij] = cX;
	}

	for (int i1 = 1; i1 < l_SNOVA - 1; i1++) {
		cX = gf16_from_nibble(c[i1]);
		for (int ij = 0; ij < lsq_SNOVA; ij++) {
			xTemp[ij] ^= cX * xS[i1][ij];
		}
	}

	uint8_t zero = ct_gf16_is_not_zero(c[rank - 1]);
	uint8_t val = zero * c[rank - 1] + (1 - zero) * (15 + ct_gf16_is_not_zero(c[0]) - c[0]);

	cX = gf16_from_nibble(val);
	for (int ij = 0; ij < lsq_SNOVA; ij++) {
		xTemp[ij] ^= cX * xS[l_SNOVA - 1][ij];
	}

	for (int ij = 0; ij < lsq_SNOVA; ij++) {
		pt_matrix[ij] = gf16_to_nibble(xTemp[ij]);
	}
}

/**
 * Generate the linear map T12
 * @param T12 - output
 * @param seed - input
 */
static void gen_seeds_and_T12(T12_t T12, const uint8_t* seed) {
	gf16_t *pt_array;
	uint8_t prng_output_private[bytes_prng_private];
	gf16_t GF16_prng_output_private[GF16s_prng_private];

	shake256(prng_output_private, bytes_prng_private, seed, seed_length_private);
	convert_bytes_to_GF16s(prng_output_private, GF16_prng_output_private, GF16s_prng_private);

	pt_array = GF16_prng_output_private;
	for (int j = 0; j < v_SNOVA; ++j) {
		for (int k = 0; k < o_SNOVA; ++k) {
			gen_a_FqS_ct(pt_array, T12[j][k]);
			pt_array += rank;
		}
	}
}

/**
 * Generate the random part of public key
 * @param map - P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param pt_public_key_seed - input
 */

static void gen_ABQ(map_group1* map, const uint8_t* pt_public_key_seed) {
	uint8_t prng_output_public[bytes_prng_public];
	uint8_t Q_temp[(sizeof(Qalpha1_t) + sizeof(Qalpha2_t)) / l_SNOVA];
	// ----- pt temp -----
	pk_expand(pt_public_key_seed, prng_output_public);
#if FIXED_ABQ
	convert_bytes_to_GF16s(prng_output_public, (uint8_t*)map, GF16s_prng_public - sizeof(Q_temp));
	memcpy(map->Aalpha, fixed_abq, 4 * m_SNOVA * alpha_SNOVA * lsq_SNOVA);
#else
	convert_bytes_to_GF16s(prng_output_public, (uint8_t*)map, GF16s_prng_public - sizeof(Q_temp));
	convert_bytes_to_GF16s(prng_output_public + sizeof(prng_output_public) - ((sizeof(Q_temp) + 1) >> 1), Q_temp,
	                       sizeof(Q_temp));

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

	gf16_t *pt_array = Q_temp;
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
 * P22 byte to GF16
 * @param P22_gf16s - output
 * @param P22_bytes - input
 */
static void input_P22(uint8_t* P22_gf16s, const uint8_t* P22_bytes) {
	convert_bytes_to_GF16s(P22_bytes, P22_gf16s, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);
}

/**
 * Pack expanded private key. esk = (key_elems, pt_private_key_seed).
 * @param esk - pointer to output expanded private key.
 * @param key_elems - pointer to input snova key elements.
 * @param pt_private_key_seed - pointer to input private key seed.
 */
static void sk_pack(uint8_t* esk, snova_key_elems* key_elems, const uint8_t* pt_private_key_seed) {
	uint8_t *sk_gf16_ptr = (uint8_t*)(key_elems->map1.Aalpha);
	convert_GF16s_to_bytes_merger_in_half(esk, sk_gf16_ptr, (bytes_sk - (seed_length_public + seed_length_private)) * 2);
	memcpy(esk + (bytes_sk - (seed_length_public + seed_length_private)), key_elems->pk.pt_public_key_seed, seed_length_public);
	memcpy(esk + (bytes_sk - seed_length_private), pt_private_key_seed, seed_length_private);
}

/**
 * Unpack expanded secret key. skupk = (esk).
 * @param skupk - pointer to output private key (unpack).
 * @param esk - pointer to input expanded private key.
 */
static void sk_unpack(sk_gf16* skupk, const uint8_t* esk) {
	convert_bytes_to_GF16s_cut_in_half(esk, (uint8_t*)skupk, (bytes_sk - (seed_length_public + seed_length_private)) * 2);
	memcpy(skupk->pt_public_key_seed, esk + (bytes_sk - (seed_length_public + seed_length_private)),
	       seed_length_public + seed_length_private);
}

/**
 * Pack public key. pk = (key_elems).
 */
static void pk_pack(uint8_t* pk, snova_key_elems* key_elems) {
	memcpy(pk, &key_elems->pk, bytes_pk);
}

/**
 * expand public key
 * @param pkx - output
 * @param pk - input
 */
static void expand_public_pack_core(uint8_t* pkx_pck, const uint8_t* pk) {
	public_key_expand* pkx_unpack = (public_key_expand*)pkx_pck;
	public_key* pk_stru = (public_key*)pk;
	memcpy(pkx_unpack->pt_public_key_seed, pk_stru->pt_public_key_seed, sizeof(pk_stru->pt_public_key_seed));
	// generate PRNG part of public key
	gen_ABQ(&(pkx_unpack->map1), pk_stru->pt_public_key_seed);
	// read  P22
	input_P22((uint8_t*)pkx_unpack->P22, (uint8_t*)pk_stru->P22);
}

/**
 * createHashOut
 */
static void createSignedHash(const uint8_t* digest, uint64_t bytes_digest, const uint8_t* pt_public_key_seed,
                             const uint8_t *array_salt, uint8_t *signed_hash_out) {
	Keccak_HashInstance hashInstance;
	Keccak_HashInitialize_SHAKE256(&hashInstance);
	Keccak_HashUpdate(&hashInstance, pt_public_key_seed, 8 * seed_length_public);
	Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
	Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
	Keccak_HashFinal(&hashInstance, NULL);
	Keccak_HashSqueeze(&hashInstance, signed_hash_out, 8 * bytes_hash);
}

#if __AVX2__
#include <immintrin.h>
typedef __m256i VECTOR;
#define VEC_LENGTH 32
#define VEC_BYTE(vect, idx) (((uint8_t*)vect)[idx])
#define VEC_SHUFFLE _mm256_shuffle_epi8
#define VEC_CMP_GT _mm256_cmpgt_epi32
#define VEC_CMP_EQ _mm256_cmpeq_epi8

#elif __ARM_NEON
#include <arm_neon.h>
typedef uint8x16_t VECTOR;
#define VEC_LENGTH 16
#define VEC_BYTE(vect, idx) vect[(idx) / 16][(idx) % 16]
#define VEC_SHUFFLE vqtbl1q_u8
#define VEC_CMP_GT vcgtq_u8
#define VEC_CMP_EQ vceqq_u8

#else
#error "Vectorization not supported"
#endif

static alignas(16 * VEC_LENGTH) uint8_t mt4b2_16[256][VEC_LENGTH];
static VECTOR *mtk2_16 = (VECTOR*)mt4b2_16;

// inverse table, runs in constant time
static VECTOR vector_inv_table = {0};
static VECTOR l_mask = {0};

// Table used by vtl_ct_multtab
static VECTOR vtl_multmask1, vtl_multmask2, vtl_multmask4, vtl_multmask8;
static VECTOR vtl_mult_table1, vtl_mult_table2, vtl_mult_table4, vtl_mult_table8;
static VECTOR zero256 = {0};

static int init_vector_table(void) {
	static int vector_table_init_flag = 0;
	if (vector_table_init_flag) {
		return 0;
	}
	vector_table_init_flag = 1;

	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 16; ++j) {
			for (int k = 0; k < 16; ++k) {
				uint8_t temp = (mt(i, k) << 4) ^ mt(j, k);
				mt4b2_16[i * 16 + j][k] = temp;
#if VEC_LENGTH > 16
				mt4b2_16[i * 16 + j][k + 16] = temp;
#endif
			}
		}
	}
#if __ARM_NEON
	// GF16 inverse table
	uint8_t inv_table[16] = {0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8};
	vector_inv_table = vld1q_u8(inv_table);

	uint8_t numf = 0x0f;
	l_mask = vld1q_dup_u8(&numf);

	uint8_t num1 = 1;
	uint8_t num2 = 2;
	uint8_t num4 = 4;
	uint8_t num8 = 8;
	vtl_multmask1 = vld1q_dup_u8(&num1);
	vtl_multmask2 = vld1q_dup_u8(&num2);
	vtl_multmask4 = vld1q_dup_u8(&num4);
	vtl_multmask8 = vld1q_dup_u8(&num8);

#else
	// GF16 inverse table
#if VEC_LENGTH > 16
	vector_inv_table = _mm256_setr_epi8(0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8, 0, 1, 9, 14, 13, 11, 7, 6, 15, 2,
	                                    12, 5, 10, 4, 3, 8);
#else
	vector_inv_table = _mm_setr_epi8(0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8);
#endif
	l_mask = _mm256_set1_epi8(0x0f);

	vtl_multmask1 = _mm256_set1_epi8(1);
	vtl_multmask2 = _mm256_set1_epi8(2);
	vtl_multmask4 = _mm256_set1_epi8(4);
	vtl_multmask8 = _mm256_set1_epi8(8);
#endif

	vtl_mult_table1 = mtk2_16[1];
	vtl_mult_table2 = mtk2_16[2];
	vtl_mult_table4 = mtk2_16[4];
	vtl_mult_table8 = mtk2_16[8];

	return 1;
}

// Constant time VTL table
static inline VECTOR vtl_ct_multtab(uint8_t val) {
#if __ARM_NEON
	VECTOR val256 = vld1q_dup_u8(&val);
#else
	VECTOR val256 = _mm256_set1_epi8(val);
#endif

	return (vtl_mult_table1 & VEC_CMP_GT(val256 & vtl_multmask1, zero256)) ^
	       (vtl_mult_table2 & VEC_CMP_GT(val256 & vtl_multmask2, zero256)) ^
	       (vtl_mult_table4 & VEC_CMP_GT(val256 & vtl_multmask4, zero256)) ^
	       (vtl_mult_table8 & VEC_CMP_GT(val256 & vtl_multmask8, zero256));
}

static inline void gf16_32_mul_k(VECTOR* a_256, gf16_t k, VECTOR* ak) {
	*ak = VEC_SHUFFLE(vtl_ct_multtab(k), *a_256);
}

static inline void gf16_32_mul_k_add(VECTOR* a_256, gf16_t k, VECTOR* ak) {
	*ak ^= VEC_SHUFFLE(vtl_ct_multtab(k), *a_256);
}

/**
 * c[i] = a[i] * b[i]
 */
static inline void gf16_32_mul_32(VECTOR* a_256, VECTOR* b_256, VECTOR* sum) {
	memset(sum, 0, sizeof(VECTOR));

	*sum ^= *b_256 & VEC_CMP_EQ(*a_256 & vtl_multmask1, vtl_multmask1);
	*sum ^= VEC_SHUFFLE(vtl_mult_table2, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask2, vtl_multmask2);
	*sum ^= VEC_SHUFFLE(vtl_mult_table4, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask4, vtl_multmask4);
	*sum ^= VEC_SHUFFLE(vtl_mult_table8, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask8, vtl_multmask8);
}

/**
 * c[i] += a[i] * b[i]
 */
static inline void gf16_32_mul_32_add(VECTOR* a_256, VECTOR* b_256, VECTOR* sum) {
	*sum ^= *b_256 & VEC_CMP_EQ(*a_256 & vtl_multmask1, vtl_multmask1);
	*sum ^= VEC_SHUFFLE(vtl_mult_table2, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask2, vtl_multmask2);
	*sum ^= VEC_SHUFFLE(vtl_mult_table4, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask4, vtl_multmask4);
	*sum ^= VEC_SHUFFLE(vtl_mult_table8, *b_256) & VEC_CMP_EQ(*a_256 & vtl_multmask8, vtl_multmask8);
}

/**
 * Establish P22 during keygen
 */
#define mol_SNOVA32 ((m_SNOVA * o_SNOVA * l_SNOVA + VEC_LENGTH - 1) / VEC_LENGTH)
#define mol_SNOVA (mol_SNOVA32 * VEC_LENGTH)

/**
 * Generate public key (P22 part), use vector vtl
 * @param outP22 - output
 * @param T12 - input
 * @param P21 - input
 * @param F12 - input
 */
static void gen_P22_vtl(P22_byte_t outP22, T12_t T12, P21_t P21, F12_t F12) {
	P22_t P22 = {0};

	alignas(VEC_LENGTH) uint8_t temp1_8[mol_SNOVA * o_SNOVA * l_SNOVA] = {0};
	alignas(VEC_LENGTH) uint8_t F12_8[mol_SNOVA * v_SNOVA * l_SNOVA] = {0};

	VECTOR* temp1_256 = (VECTOR*)temp1_8;
	VECTOR* F12_256 = (VECTOR*)F12_8;

	for (int di = 0; di < v_SNOVA; ++di)
		for (int k1 = 0; k1 < l_SNOVA; ++k1)
			for (int mi = 0; mi < m_SNOVA; ++mi)
				for (int dk = 0; dk < o_SNOVA; ++dk)
					for (int j1 = 0; j1 < l_SNOVA; ++j1)
						F12_8[(di * l_SNOVA + k1) * mol_SNOVA + mi * o_SNOVA * l_SNOVA + dk * l_SNOVA + j1] =
						    F12[mi][di][dk][k1 * l_SNOVA + j1];

	for (int di = 0; di < v_SNOVA; ++di)
		for (int dj = 0; dj < o_SNOVA; ++dj)
			for (int i1 = 0; i1 < l_SNOVA; ++i1)
				for (int k1 = 0; k1 < l_SNOVA; ++k1) {
					VECTOR k_lh = vtl_ct_multtab(T12[di][dj][i1 * l_SNOVA + k1]);

					for (int mi_dk_j1 = 0; mi_dk_j1 < mol_SNOVA32; mi_dk_j1++)
						temp1_256[(dj * l_SNOVA + i1) * mol_SNOVA32 + mi_dk_j1] ^=
						    VEC_SHUFFLE(k_lh, F12_256[(di * l_SNOVA + k1) * mol_SNOVA32 + mi_dk_j1]);
				}

	for (int mi = 0; mi < m_SNOVA; ++mi)
		for (int dj = 0; dj < o_SNOVA; ++dj)
			for (int dk = 0; dk < o_SNOVA; ++dk)
				for (int i1 = 0; i1 < l_SNOVA; ++i1)
					for (int j1 = 0; j1 < l_SNOVA; ++j1)
						P22[mi][dj][dk][i1 * l_SNOVA + j1] ^=
						    temp1_8[(dj * l_SNOVA + i1) * mol_SNOVA + mi * o_SNOVA * l_SNOVA + dk * l_SNOVA + j1];

	alignas(VEC_LENGTH) uint8_t temp2_8[mol_SNOVA * o_SNOVA * l_SNOVA] = {0};
	alignas(VEC_LENGTH) uint8_t P21_8[mol_SNOVA * v_SNOVA * l_SNOVA] = {0};

	VECTOR* temp2_256 = (VECTOR*)temp2_8;
	VECTOR* P21_256 = (VECTOR*)P21_8;

	for (int di = 0; di < v_SNOVA; ++di)
		for (int k1 = 0; k1 < l_SNOVA; ++k1)
			for (int mi = 0; mi < m_SNOVA; ++mi)
				for (int dj = 0; dj < o_SNOVA; ++dj)
					for (int i1 = 0; i1 < l_SNOVA; ++i1)
						P21_8[(di * l_SNOVA + k1) * mol_SNOVA + mi * o_SNOVA * l_SNOVA + dj * l_SNOVA + i1] =
						    P21[mi][dj][di][i1 * l_SNOVA + k1];

	for (int di = 0; di < v_SNOVA; ++di)
		for (int dk = 0; dk < o_SNOVA; ++dk)
			for (int k1 = 0; k1 < l_SNOVA; ++k1)
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					VECTOR k_lh = vtl_ct_multtab(T12[di][dk][k1 * l_SNOVA + j1]);

					for (int mi_dj_i1 = 0; mi_dj_i1 < mol_SNOVA32; mi_dj_i1++)
						temp2_256[(dk * l_SNOVA + j1) * mol_SNOVA32 + mi_dj_i1] ^=
						    VEC_SHUFFLE(k_lh, P21_256[(di * l_SNOVA + k1) * mol_SNOVA32 + mi_dj_i1]);
				}

	for (int mi = 0; mi < m_SNOVA; ++mi)
		for (int dj = 0; dj < o_SNOVA; ++dj)
			for (int dk = 0; dk < o_SNOVA; ++dk)
				for (int i1 = 0; i1 < l_SNOVA; ++i1)
					for (int j1 = 0; j1 < l_SNOVA; ++j1)
						P22[mi][dj][dk][i1 * l_SNOVA + j1] ^=
						    temp2_8[(dk * l_SNOVA + j1) * mol_SNOVA + mi * o_SNOVA * l_SNOVA + dj * l_SNOVA + i1];

	convert_GF16s_to_bytes(outP22, (uint8_t*)P22, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);
}

/**
 * Establish F during keygen
 */
#define mvl_SNOVA32 ((m_SNOVA * v_SNOVA * l_SNOVA + VEC_LENGTH - 1) / VEC_LENGTH)
#define mvl_SNOVA8 (mvl_SNOVA32 * VEC_LENGTH)
#define vl_SNOVA (v_SNOVA * l_SNOVA)

/**
 * Generate private key (F part), use vector vtl
 * @param map2 - output: F11 F12 F21
 * @param map1 - input: P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param T12 - input
 */
static void gen_F_vtl(map_group2* map2, map_group1* map1, T12_t T12) {
	VECTOR p11_256[mvl_SNOVA32 * vl_SNOVA] = {0};
	VECTOR t12_256[mvl_SNOVA32 * l_SNOVA] = {0};
	VECTOR res256[mvl_SNOVA32 * o_SNOVA * l_SNOVA] = {0};

	uint8_t *p11_8 = (uint8_t*)p11_256;
	uint8_t *t12_8 = (uint8_t*)t12_256;
	uint8_t *res8 = (uint8_t*)res256;

	memcpy(map2->F11, map1->P11, m_SNOVA * v_SNOVA * v_SNOVA * lsq_SNOVA);
	memcpy(map2->F12, map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * lsq_SNOVA);
	memcpy(map2->F21, map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * lsq_SNOVA);

	// F12

	for (int k1 = 0; k1 < l_SNOVA; ++k1)
		for (int dk = 0; dk < v_SNOVA; ++dk)
			for (int dj = 0; dj < o_SNOVA; ++dj)
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					t12_8[(dk * l_SNOVA + k1) * o_SNOVA * l_SNOVA + dj * l_SNOVA + j1] = T12[dk][dj][k1 * l_SNOVA + j1];
				}

	for (int dk = 0; dk < v_SNOVA; dk++)
		for (int di = 0; di < v_SNOVA; di++)
			for (int j1 = 0; j1 < l_SNOVA; ++j1)
				for (int mi = 0; mi < m_SNOVA; mi++)
					for (int i1 = 0; i1 < l_SNOVA; ++i1)
						p11_8[(dk * l_SNOVA + j1) * mvl_SNOVA8 + mi * v_SNOVA * l_SNOVA + di * l_SNOVA + i1] =
						    map1->P11[mi][di][dk][i1 * l_SNOVA + j1];

	for (int dj_j1 = 0; dj_j1 < o_SNOVA * l_SNOVA; ++dj_j1)
		for (int dk_k1 = 0; dk_k1 < v_SNOVA * l_SNOVA; ++dk_k1) {
			VECTOR k_lh = vtl_ct_multtab(t12_8[dk_k1 * o_SNOVA * l_SNOVA + dj_j1]);
			for (int mi_di_i1 = 0; mi_di_i1 < mvl_SNOVA32; mi_di_i1++) {
				res256[dj_j1 * mvl_SNOVA32 + mi_di_i1] ^= VEC_SHUFFLE(k_lh, p11_256[dk_k1 * mvl_SNOVA32 + mi_di_i1]);
			}
		}

	for (int di = 0; di < v_SNOVA; ++di)
		for (int dj = 0; dj < o_SNOVA; ++dj)
			for (int j1 = 0; j1 < l_SNOVA; ++j1)
				for (int i1 = 0; i1 < l_SNOVA; ++i1)
					for (int mi = 0; mi < m_SNOVA; ++mi)
						map2->F12[mi][di][dj][i1 * l_SNOVA + j1] ^=
						    res8[dj * l_SNOVA * mvl_SNOVA8 + j1 * mvl_SNOVA8 + mi * v_SNOVA * l_SNOVA + di * l_SNOVA + i1];

	// Same for F21
	memset(res256, 0, sizeof(res256));

	for (int dj = 0; dj < v_SNOVA; ++dj)
		for (int dk = 0; dk < o_SNOVA; ++dk)
			for (int i1 = 0; i1 < l_SNOVA; ++i1)
				for (int j1 = 0; j1 < l_SNOVA; ++j1) {
					t12_8[(dj * l_SNOVA + j1) * o_SNOVA * l_SNOVA + dk * l_SNOVA + i1] = T12[dj][dk][i1 * l_SNOVA + j1];
				}

	for (int mi = 0; mi < m_SNOVA; mi++)
		for (int di = 0; di < v_SNOVA; di++)
			for (int dk = 0; dk < v_SNOVA; dk++)
				for (int i1 = 0; i1 < l_SNOVA; ++i1)
					for (int j1 = 0; j1 < l_SNOVA; ++j1)
						p11_8[(di * l_SNOVA + i1) * mvl_SNOVA8 + mi * v_SNOVA * l_SNOVA + dk * l_SNOVA + j1] =
						    map1->P11[mi][di][dk][i1 * l_SNOVA + j1];

	for (int dj_j1 = 0; dj_j1 < o_SNOVA * l_SNOVA; ++dj_j1)
		for (int dk_k1 = 0; dk_k1 < v_SNOVA * l_SNOVA; ++dk_k1) {
			VECTOR k_lh = vtl_ct_multtab(t12_8[dk_k1 * o_SNOVA * l_SNOVA + dj_j1]);
			for (int mi_di_i1 = 0; mi_di_i1 < mvl_SNOVA32; mi_di_i1++) {
				res256[dj_j1 * mvl_SNOVA32 + mi_di_i1] ^= VEC_SHUFFLE(k_lh, p11_256[dk_k1 * mvl_SNOVA32 + mi_di_i1]);
			}
		}

	// Shuffle back
	for (int mi = 0; mi < m_SNOVA; ++mi)
		for (int dj = 0; dj < o_SNOVA; ++dj)
			for (int di = 0; di < v_SNOVA; ++di)
				for (int i1 = 0; i1 < l_SNOVA; ++i1)
					for (int j1 = 0; j1 < l_SNOVA; ++j1)
						map2->F21[mi][dj][di][i1 * l_SNOVA + j1] ^=
						    res8[dj * l_SNOVA * mvl_SNOVA8 + i1 * mvl_SNOVA8 + mi * v_SNOVA * l_SNOVA + di * l_SNOVA + j1];
}

static inline uint8_t horizontal_xor_256(VECTOR vec) {
#if __ARM_NEON
	int64x2_t combined = (int64x2_t)vec;
#elif VEC_LENGTH > 16
	__m128i low = _mm256_castsi256_si128(vec);
	__m128i high = _mm256_extracti128_si256(vec, 1);
	__m128i combined = _mm_xor_si128(low, high);
#else
	__m128i combined = vec;
#endif

	combined ^= combined >> 32;
	combined ^= combined >> 16;
	combined ^= combined >> 8;

	return combined[0] ^ combined[1];
}

#define vtl_len ((n_SNOVA * rank + VEC_LENGTH - 1) / VEC_LENGTH)
#define vtl_v_len ((v_SNOVA * rank + VEC_LENGTH - 1) / VEC_LENGTH)
#define vtl_o_len ((o_SNOVA * rank + VEC_LENGTH - 1) / VEC_LENGTH)
#define vtl_mainRow_x_rank(mainRow) (((mainRow) * rank + VEC_LENGTH - 1) / VEC_LENGTH * VEC_LENGTH)
#define vtl_mainRow_x_rank32(mainRow) (((mainRow) * rank + VEC_LENGTH - 1) / VEC_LENGTH)
#define rank_next32 ((rank + 1) / 2)
#define rank_next2 (rank_next32 * 2)
#define rank_floor2 ((rank) / 2 * 2)
#define alpha_SNOVA_next2 ((alpha_SNOVA + 1) / 2 * 2)

// Optimize vector alignment. Only useful for even rank
#if rank % 2
#define MJ_MAX mainRow
#else
#define MJ_MAX (vtl_mainRow_x_rank(mainRow)) / rank
#endif

static inline void jogressMatrix_vector(VECTOR* AJ, const uint8_t* A, int mainCol, int mainRow) {
	for (int mi = 0; mi < mainCol; mi++) {
		for (int mj = 0; mj < mainRow; mj++) {
			for (int ei = 0; ei < rank; ei++) {
				for (int ej = 0; ej < rank; ej++) {
					// AJ[((mi * rank + ei) * (mainRow * rank)) + (mj * rank + ej)] =
					VEC_BYTE(AJ, ((mi * rank + ei) * vtl_mainRow_x_rank(mainRow)) + (mj * rank + ej)) =
					    A[(mi * (mainRow) * (rank * rank)) + (mj * (rank * rank)) + (ei * rank + ej)];
				}
			}
		}
	}
}

static inline void jogressTrMatrix_vector(VECTOR* AJ, const uint8_t* A, int mainCol, int mainRow) {
	for (int mi = 0; mi < mainCol; mi++) {
		for (int mj = 0; mj < mainRow; mj++) {
			for (int ei = 0; ei < rank; ei++) {
				for (int ej = 0; ej < rank; ej++) {
					// AJ[((mi * rank + ei) * (mainRow * rank)) + (mj * rank + ej)] =
					VEC_BYTE(AJ, ((mi * rank + ei) * vtl_mainRow_x_rank(mainRow)) + (mj * rank + ej)) =
					    A[(mi * (mainRow) * (rank * rank)) + (mj * (rank * rank)) + (ej * rank + ei)];
					// ei <<---->> ej
				}
			}
		}
	}
}

static inline void jogressMatrixTr_vector(VECTOR* AJ_tr, const VECTOR* AJ, int mainCol, int mainRow) {
	for (int mi = 0; mi < mainCol; mi++) {
		for (int mj = 0; mj < MJ_MAX; mj++) {
			for (int ei = 0; ei < rank; ei++) {
				for (int ej = 0; ej < rank; ej++) {
					VEC_BYTE(AJ_tr, ((mi * rank + ei) * vtl_mainRow_x_rank(mainRow)) + (mj * rank + ej)) =
					    VEC_BYTE(AJ, ((mi * rank + ej) * vtl_mainRow_x_rank(mainRow)) + (mj * rank + ei));
					//  ei <<------------------------------------------>> ej
				}
			}
		}
	}
}

/**
 * CT version for any l_SNOVA.
 */

#define GAUSS_ROW (m_SNOVA * lsq_SNOVA + 1)
#define GAUSS_ROW32 ((GAUSS_ROW + VEC_LENGTH - 1) / VEC_LENGTH)
#define GAUSS_ROW_mult32 (GAUSS_ROW32 * VEC_LENGTH)

static void calc_LR_J_vtl(uint8_t L_J[m_SNOVA][alpha_SNOVA][rank][vtl_v_len * VEC_LENGTH],
                          uint8_t R_tr_J[m_SNOVA][alpha_SNOVA][rank][vtl_v_len * VEC_LENGTH], Aalpha_t Aalpha, Balpha_t Balpha,
                          Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, gf16m_t *X_in_GF16Matrix) {
	alignas(VEC_LENGTH) VECTOR X_J[rank][vtl_v_len];
	alignas(VEC_LENGTH) VECTOR X_tr_J[rank][vtl_v_len];

	jogressMatrix_vector((VECTOR*)X_J, (uint8_t*)X_in_GF16Matrix, 1, v_SNOVA);
	jogressTrMatrix_vector((VECTOR*)X_tr_J, (uint8_t*)X_in_GF16Matrix, 1, v_SNOVA);

	// calc LR
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			int mi_prime = i_prime_inv(mi, alpha);

			alignas(VEC_LENGTH) VECTOR AxS_256[rank][vtl_v_len] = {0};
			alignas(VEC_LENGTH) VECTOR AxS_tr_256[rank][vtl_v_len] = {0};
			alignas(VEC_LENGTH) VECTOR Q2xS_256[rank][vtl_v_len] = {0};
			alignas(VEC_LENGTH) VECTOR Q2xS_tr_256[rank][vtl_v_len] = {0};
			alignas(VEC_LENGTH) VECTOR L_tr_J_256[rank][vtl_v_len] = {0};

			for (int ni = 0; ni < rank; ++ni) {
				for (int nj = 0; nj < rank; ++nj) {
					VECTOR k1 = mtk2_16[get_gf16m(Aalpha[mi_prime][alpha], ni, nj)];
					VECTOR k2 = mtk2_16[get_gf16m(Qalpha2[mi_prime][alpha], ni, nj)];
					for (int nk = 0; nk < vtl_v_len; ++nk) {
						AxS_256[ni][nk] ^= VEC_SHUFFLE(k1, X_tr_J[nj][nk]);
						Q2xS_256[ni][nk] ^= VEC_SHUFFLE(k2, X_J[nj][nk]);
					}
				}
			}
			jogressMatrixTr_vector((VECTOR*)AxS_tr_256, (VECTOR*)AxS_256, 1, v_SNOVA);
			jogressMatrixTr_vector((VECTOR*)Q2xS_tr_256, (VECTOR*)Q2xS_256, 1, v_SNOVA);

			for (int ni = 0; ni < rank; ++ni) {
				for (int nj = 0; nj < rank; ++nj) {
					VECTOR k1 = mtk2_16[get_gf16m(Qalpha1[mi_prime][alpha], ni, nj)];
					VECTOR k2 = mtk2_16[get_gf16m(Balpha[mi_prime][alpha], nj, ni)];
					for (int nk = 0; nk < vtl_v_len; ++nk) {
						VECTOR* R_tr_J_256 = (VECTOR*)R_tr_J[mi][alpha][ni];
						L_tr_J_256[ni][nk] ^= VEC_SHUFFLE(k1, AxS_tr_256[nj][nk]);
						R_tr_J_256[nk] ^= VEC_SHUFFLE(k2, Q2xS_tr_256[nj][nk]);
					}
				}
			}
			jogressMatrixTr_vector((VECTOR*)L_J[mi][alpha], (VECTOR*)L_tr_J_256, 1, v_SNOVA);
		}
	}
}

/**
 * Computes signature
 */
static int sign_digest_core_gnl_vtl(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, uint8_t* array_salt,
                                    Aalpha_t Aalpha, Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, T12_t T12,
                                    F11_t F11, F12_t F12, F21_t F21, const uint8_t pt_public_key_seed[seed_length_public],
                                    const uint8_t pt_private_key_seed[seed_length_private]) {
	alignas(VEC_LENGTH) uint8_t Gauss[m_SNOVA * lsq_SNOVA][GAUSS_ROW_mult32];

	gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
	gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
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
		memset(Gauss, 0, sizeof(Gauss));
		num_sign++;
		if (num_sign == 255) {
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

		convert_bytes_to_GF16s(vinegar_in_byte, (uint8_t*)X_in_GF16Matrix, v_SNOVA * lsq_SNOVA);

		alignas(VEC_LENGTH) uint8_t L_J[m_SNOVA][alpha_SNOVA][rank][vtl_v_len * VEC_LENGTH] = {0};
		alignas(VEC_LENGTH) uint8_t R_tr_J[m_SNOVA][alpha_SNOVA][rank][vtl_v_len * VEC_LENGTH] = {0};

		calc_LR_J_vtl(L_J, R_tr_J, Aalpha, Balpha, Qalpha1, Qalpha2, X_in_GF16Matrix);

		for (int mi = 0; mi < m_SNOVA; ++mi) {
			gf16m_set_zero(Fvv_in_GF16Matrix[mi]);
		}
		for (int mi = 0; mi < m_SNOVA; ++mi) {
			VECTOR F11_J_256[v_SNOVA * rank][vtl_v_len] = {0};
			jogressMatrix_vector((VECTOR*)F11_J_256, (uint8_t*)F11[mi], v_SNOVA, v_SNOVA);

			for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
				int mi_prime = i_prime_inv(mi, alpha);
				VECTOR LJxF11J_256[rank][vtl_v_len] = {0};
				for (int vi = 0; vi < rank; ++vi) {
					for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
						VECTOR k = vtl_ct_multtab(L_J[mi][alpha][vi][vj]);
						for (int vk = 0; vk < vtl_v_len; ++vk) {
							LJxF11J_256[vi][vk] ^= VEC_SHUFFLE(k, F11_J_256[vj][vk]);
						}
					}
				}

				uint8_t LJxF11JxRJ[rank][rank] = {0};
				for (int vi = 0; vi < rank; ++vi) {
					for (int vj = 0; vj < rank; ++vj) {
						VECTOR* R_tr_J_256 = (VECTOR*)R_tr_J[mi][alpha][vj];
						VECTOR tmp_256 = {0};
						for (int vk = 0; vk < vtl_v_len; ++vk) {
							gf16_32_mul_32_add((LJxF11J_256[vi] + vk), (R_tr_J_256 + vk), (&tmp_256));
						}

						LJxF11JxRJ[vi][vj] ^= horizontal_xor_256(tmp_256);
					}
				}

				// hash_in_GF16Matrix[mi] += LJxPJxRJ
				for (int ni = 0; ni < rank; ++ni) {
					for (int nj = 0; nj < rank; ++nj) {
						gf16_t t = get_gf16m(Fvv_in_GF16Matrix[mi_prime], ni, nj);
						set_gf16m(Fvv_in_GF16Matrix[mi_prime], ni, nj, t ^ LJxF11JxRJ[ni][nj]);
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
		alignas(VEC_LENGTH) uint8_t Temp1[o_SNOVA * l_SNOVA * lsq_SNOVA * vtl_mainRow_x_rank(o_SNOVA)] = {0};
		alignas(VEC_LENGTH) uint8_t Temp2[o_SNOVA * l_SNOVA * lsq_SNOVA * vtl_mainRow_x_rank(o_SNOVA)] = {0};

		for (int mi = 0; mi < m_SNOVA; ++mi) {
			gf16m_t F21_vo[v_SNOVA][o_SNOVA];
			VECTOR F12_J_256[v_SNOVA * rank][vtl_o_len] = {0};
			VECTOR F21_vo_tr_J_256[v_SNOVA * rank][vtl_o_len] = {0};

			// swap F21  v <-> o
			for (int oi = 0; oi < o_SNOVA; ++oi) {
				for (int vi = 0; vi < v_SNOVA; ++vi) {
					gf16m_clone(F21_vo[vi][oi], F21[mi][oi][vi]);
				}
			}

			jogressMatrix_vector((VECTOR*)F12_J_256, (uint8_t*)F12[mi], v_SNOVA, o_SNOVA);
			jogressTrMatrix_vector((VECTOR*)F21_vo_tr_J_256, (uint8_t*)F21_vo, v_SNOVA, o_SNOVA);

			for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
				int mi_prime_inv = i_prime_inv(mi, alpha);

				VECTOR LJxF12J_256[rank][vtl_o_len] = {0};

				for (int vi = 0; vi < rank; ++vi) {
					for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
						VECTOR k = vtl_ct_multtab(L_J[mi][alpha][vi][vj]);
						for (int oi = 0; oi < vtl_o_len; ++oi) {
							LJxF12J_256[vi][oi] ^= VEC_SHUFFLE(k, F12_J_256[vj][oi]);
						}
					}
				}

				VECTOR LJxF12J_tr_256[rank][vtl_o_len] = {0};
				VECTOR LJxF12JxQJ_tr_256[rank][vtl_o_len] = {0};

				jogressMatrixTr_vector((VECTOR*)LJxF12J_tr_256, (VECTOR*)LJxF12J_256, 1, o_SNOVA);

				for (int ri = 0; ri < rank; ++ri) {
					for (int rj = 0; rj < rank; ++rj) {
						VECTOR k = mtk2_16[get_gf16m(Qalpha2[mi_prime_inv][alpha], rj, ri)];
						for (int oi = 0; oi < vtl_o_len; ++oi) {
							LJxF12JxQJ_tr_256[ri][oi] ^= VEC_SHUFFLE(k, LJxF12J_tr_256[rj][oi]);
						}
					}
				}

				for (int ti2 = 0; ti2 < l_SNOVA; ++ti2)
					for (int tj2 = 0; tj2 < l_SNOVA; ++tj2) {
						VECTOR k = mtk2_16[Balpha[mi_prime_inv][alpha][tj2 * rank + ti2]];
						VECTOR* Temp1_256 = (VECTOR*)Temp1;
						VECTOR* AJ_256 = (VECTOR*)LJxF12JxQJ_tr_256;

						for (int tj1 = 0; tj1 < l_SNOVA; tj1++)
							for (int toi = 0; toi < vtl_mainRow_x_rank32(o_SNOVA); toi++) {
								Temp1_256[((mi_prime_inv * lsq_SNOVA + ti2 * rank + tj2) * l_SNOVA + tj1) *
								          vtl_mainRow_x_rank32(o_SNOVA) +
								          toi] ^= VEC_SHUFFLE(k, AJ_256[(tj1 * vtl_mainRow_x_rank32(o_SNOVA)) + toi]);
							}
					}

				// ------- ^^^ F12    F22 vvv -------
				VECTOR F21JxRJ_256[rank][vtl_o_len] = {0};
				VECTOR RTRJ_x_F21TRJ_256[rank][vtl_o_len] = {0};

				for (int vi = 0; vi < rank; ++vi) {
					for (int vj = 0; vj < v_SNOVA * rank; ++vj) {
						VECTOR k = vtl_ct_multtab(R_tr_J[mi][alpha][vi][vj]);
						for (int oi = 0; oi < vtl_o_len; ++oi) {
							RTRJ_x_F21TRJ_256[vi][oi] ^= VEC_SHUFFLE(k, F21_vo_tr_J_256[vj][oi]);
						}
					}
				}
				jogressMatrixTr_vector((VECTOR*)F21JxRJ_256, (VECTOR*)RTRJ_x_F21TRJ_256, 1, o_SNOVA);

				VECTOR Q1xF21JxRJ_256[rank][vtl_o_len] = {0};
				for (int ri = 0; ri < rank; ++ri) {
					for (int rj = 0; rj < rank; ++rj) {
						VECTOR k = mtk2_16[get_gf16m(Qalpha1[mi_prime_inv][alpha], ri, rj)];
						for (int oi = 0; oi < vtl_o_len; ++oi) {
							Q1xF21JxRJ_256[ri][oi] ^= VEC_SHUFFLE(k, F21JxRJ_256[rj][oi]);
						}
					}
				}

				for (int ti1 = 0; ti1 < l_SNOVA; ++ti1)
					for (int tj2 = 0; tj2 < l_SNOVA; ++tj2) {
						VECTOR k = mtk2_16[Aalpha[mi_prime_inv][alpha][ti1 * rank + tj2]];
						VECTOR* Temp2_256 = (VECTOR*)Temp2;
						VECTOR* AJ_256 = (VECTOR*)Q1xF21JxRJ_256;

						for (int tj1 = 0; tj1 < l_SNOVA; ++tj1)
							for (int toi = 0; toi < vtl_mainRow_x_rank32(o_SNOVA); toi++) {
								Temp2_256[((mi_prime_inv * lsq_SNOVA + ti1 * l_SNOVA + tj2) * l_SNOVA + tj1) *
								          vtl_mainRow_x_rank32(o_SNOVA) +
								          toi] ^= VEC_SHUFFLE(k, AJ_256[(tj1 * vtl_mainRow_x_rank32(o_SNOVA)) + toi]);
							}
					}
			}
		}

		for (int mi_prime_inv = 0; mi_prime_inv < o_SNOVA; ++mi_prime_inv)
			for (int ti1 = 0; ti1 < l_SNOVA; ++ti1)
				for (int ti2 = 0; ti2 < l_SNOVA; ++ti2)
					for (int tj1 = 0; tj1 < l_SNOVA; ++tj1)
						for (int oi = 0; oi < o_SNOVA; ++oi)
							for (int tj2 = 0; tj2 < l_SNOVA; ++tj2)
								Gauss[mi_prime_inv * lsq_SNOVA + ti1 * rank + ti2][oi * lsq_SNOVA + tj1 * rank + tj2] ^=
								    Temp1[((mi_prime_inv * lsq_SNOVA + ti2 * rank + tj2) * l_SNOVA + tj1) *
								          vtl_mainRow_x_rank(o_SNOVA) +
								          (oi * rank + ti1)];

		for (int mi_prime_inv = 0; mi_prime_inv < o_SNOVA; ++mi_prime_inv)
			for (int oi = 0; oi < o_SNOVA; ++oi)
				for (int ti1 = 0; ti1 < l_SNOVA; ++ti1)
					for (int ti2 = 0; ti2 < l_SNOVA; ++ti2)
						for (int tj1 = 0; tj1 < l_SNOVA; ++tj1)
							for (int tj2 = 0; tj2 < l_SNOVA; ++tj2)
								Gauss[mi_prime_inv * lsq_SNOVA + ti1 * rank + ti2][oi * lsq_SNOVA + tj1 * rank + tj2] ^=
								    Temp2[((mi_prime_inv * lsq_SNOVA + ti1 * l_SNOVA + tj2) * l_SNOVA + tj1) *
								          vtl_mainRow_x_rank(o_SNOVA) +
								          (oi * rank + ti2)];

		// Gaussian elimination in constant time
		for (int mi2 = 0; mi2 < m_SNOVA * lsq_SNOVA; ++mi2) {
			int swap = ct_gf16_is_not_zero(Gauss[mi2][mi2]) - 1;
			for (int j2 = mi2 + 1; j2 < m_SNOVA * lsq_SNOVA; ++j2) {
#if __ARM_NEON
				int8x16_t swap256 = vld1q_dup_s8((int8_t*)&swap);
#else
				VECTOR swap256 = _mm256_set1_epi32(swap);
#endif
				VECTOR* gdest = (VECTOR*)&Gauss[mi2][0];
				VECTOR* gsource = (VECTOR*)&Gauss[j2][0];
				for (int k2 = 0; k2 < GAUSS_ROW32; ++k2) {
					gdest[k2] ^= gsource[k2] & swap256;
				}

				swap = ct_gf16_is_not_zero(Gauss[mi2][mi2]) - 1;
			}
			flag_redo |= swap;

			// Constant time GF16 inverse
			VECTOR res256[1];
			VECTOR g256[1];
#if __ARM_NEON
			g256[0] = vld1q_dup_u8(&Gauss[mi2][mi2]);
#else
			g256[0] = _mm256_set1_epi8(Gauss[mi2][mi2]);
#endif
			res256[0] = VEC_SHUFFLE(vector_inv_table, g256[0]);
			uint8_t t_GF16 = VEC_BYTE(res256, 0);

			int kstart = (mi2 / VEC_LENGTH) * VEC_LENGTH;
			for (int k = kstart; k < GAUSS_ROW_mult32; k += VEC_LENGTH) {
				gf16_32_mul_k((VECTOR*)(Gauss[mi2] + k), t_GF16, (VECTOR*)(Gauss[mi2] + k));
			}

			for (int j2 = mi2 + 1; j2 < m_SNOVA * lsq_SNOVA; ++j2) {
				t_GF16 = Gauss[j2][mi2];
				for (int k2 = kstart; k2 < GAUSS_ROW_mult32; k2 += VEC_LENGTH) {
					gf16_32_mul_k_add((VECTOR*)(Gauss[mi2] + k2), t_GF16, (VECTOR*)(Gauss[j2] + k2));
				}
			}
		}
	} while (flag_redo);

	alignas(VEC_LENGTH) uint8_t solution[GAUSS_ROW_mult32] = {0};
	uint32_t xSolution[m_SNOVA * lsq_SNOVA] = {0};
	uint8_t t_GF16 = 0;
	uint8_t Gauss_last_col;

	for (int mil2 = m_SNOVA * lsq_SNOVA - 1; mil2 >= 0; --mil2) {
		Gauss_last_col = Gauss[mil2][m_SNOVA * lsq_SNOVA];
		VECTOR t_GF16_256 = {0};

		Gauss[mil2][m_SNOVA * lsq_SNOVA] = 0;

		int kstart = ((mil2 + 1) / VEC_LENGTH) * VEC_LENGTH;
		for (int k2 = kstart; k2 < GAUSS_ROW_mult32; k2 += VEC_LENGTH) {
			gf16_32_mul_32_add((VECTOR*)(&Gauss[mil2][k2]), (VECTOR*)(&solution[k2]), &t_GF16_256);
		}

		uint8_t *t_GF16_256_8_ptr = (uint8_t*)(&t_GF16_256);
		for (int k2 = 0; k2 < VEC_LENGTH; k2++) {
			t_GF16 ^= t_GF16_256_8_ptr[k2];
		}

		solution[mil2] = Gauss_last_col ^ t_GF16;
		t_GF16 = 0;
	}

	for (int mi2 = m_SNOVA * lsq_SNOVA - 1; mi2 >= 0; --mi2) {
		xSolution[mi2] = gf16_from_nibble(solution[mi2]);
	}

	uint32_t xT12[v_SNOVA][o_SNOVA][lsq_SNOVA] = {0};

	for (int dj = 0; dj < v_SNOVA; ++dj)
		for (int dk = 0; dk < o_SNOVA; ++dk)
			for (int idx = 0; idx < lsq_SNOVA; ++idx) {
				xT12[dj][dk][idx] = gf16_from_nibble(T12[dj][dk][idx]);
			}

	// Establish Signature
	uint32_t xSig[lsq_SNOVA] = {0};
	for (int dj = 0; dj < v_SNOVA; ++dj) {
		for (int dk = 0; dk < o_SNOVA; ++dk)
			for (int i1 = 0; i1 < l_SNOVA; ++i1)
				for (int j1 = 0; j1 < l_SNOVA; ++j1)
					for (int k1 = 0; k1 < l_SNOVA; ++k1) {
						xSig[i1 * l_SNOVA + j1] ^=
						    xT12[dj][dk][i1 * l_SNOVA + k1] * xSolution[dk * lsq_SNOVA + k1 * l_SNOVA + j1];
					}

		for (int idx = 0; idx < lsq_SNOVA; ++idx) {
			signature_in_GF16Matrix[dj][idx] = X_in_GF16Matrix[dj][idx] ^ gf16_to_nibble(xSig[idx]);
		}

		memset(xSig, 0, sizeof(xSig));
	}

	for (int index = 0; index < o_SNOVA; ++index)
		for (int idx = 0; idx < lsq_SNOVA; ++idx) {
			signature_in_GF16Matrix[v_SNOVA + index][idx] = solution[index * lsq_SNOVA + idx];
		}

	// output signature
	convert_GF16s_to_bytes(pt_signature, (gf16_t*)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);
	for (int i1 = 0; i1 < bytes_salt; ++i1) {
		pt_signature[bytes_signature + i1] = array_salt[i1];
	}

	return 0;
}

/**
 * VTL version for any l_SNOVA.
 */

static inline void evaluation_vector_vtl(gf16m_t* restrict hash_in_GF16Matrix, const public_key_expand* restrict pkx,
        gf16m_t *restrict signature_in_GF16Matrix) {
	VECTOR S_J_256[rank][vtl_len];
	VECTOR S_tr_J_256[rank][vtl_len];
	VECTOR L_J_nibble256[m_SNOVA * alpha_SNOVA * rank_next32 * vtl_len] = {0};
	VECTOR R_tr_J_256[m_SNOVA][alpha_SNOVA][rank][vtl_len] = {0};

	jogressMatrix_vector((VECTOR*)S_J_256, (uint8_t*)signature_in_GF16Matrix, 1, n_SNOVA);
	jogressTrMatrix_vector((VECTOR*)S_tr_J_256, (uint8_t*)signature_in_GF16Matrix, 1, n_SNOVA);

	// calc LR
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		VECTOR AxS_tr_256[alpha_SNOVA_next2][rank][vtl_len] = {0};
		VECTOR Q2xS_tr_256[alpha_SNOVA_next2][rank][vtl_len] = {0};

		// A and Q2, vtl 2 alpha set end 0
		gf16m_t Aalpha[alpha_SNOVA_next2] = {0};
		gf16m_t Qalpha2[alpha_SNOVA_next2] = {0};
		memcpy((uint8_t*)Aalpha, (uint8_t*)(pkx->map1.Aalpha[mi]), sizeof(pkx->map1.Aalpha[mi]));
		memcpy((uint8_t*)Qalpha2, (uint8_t*)(pkx->map1.Qalpha2[mi]), sizeof(pkx->map1.Qalpha2[mi]));

		for (int alpha = 0; alpha < alpha_SNOVA_next2; alpha += 2) {
			VECTOR AxS_256[rank][vtl_len] = {0};
			VECTOR Q2xS_256[rank][vtl_len] = {0};

			// use vtl (2 alpha x vtl_len)
			for (int ni = 0; ni < rank; ++ni) {
				for (int nj = 0; nj < rank; ++nj) {
					VECTOR k1_lh = mtk2_16[get_gf16m(Aalpha[alpha], ni, nj) | (get_gf16m(Aalpha[alpha + 1], ni, nj) << 4)];
					VECTOR k2_lh = mtk2_16[get_gf16m(Qalpha2[alpha], ni, nj) | (get_gf16m(Qalpha2[alpha + 1], ni, nj) << 4)];
					for (int nk = 0; nk < vtl_len; ++nk) {
						AxS_256[ni][nk] ^= VEC_SHUFFLE(k1_lh, S_tr_J_256[nj][nk]);
						Q2xS_256[ni][nk] ^= VEC_SHUFFLE(k2_lh, S_J_256[nj][nk]);
					}
				}
			}

			jogressMatrixTr_vector((VECTOR*)AxS_tr_256[alpha], (VECTOR*)AxS_256, 1, n_SNOVA);
			jogressMatrixTr_vector((VECTOR*)Q2xS_tr_256[alpha], (VECTOR*)Q2xS_256, 1, n_SNOVA);

			// nibble splite
			for (int ni = 0; ni < rank; ++ni) {
				for (int nk = 0; nk < vtl_len; ++nk) {
					AxS_tr_256[alpha + 1][ni][nk] = (AxS_tr_256[alpha][ni][nk] >> 4) & l_mask;
					AxS_tr_256[alpha][ni][nk] &= l_mask;

					Q2xS_tr_256[alpha + 1][ni][nk] = (Q2xS_tr_256[alpha][ni][nk] >> 4) & l_mask;
					Q2xS_tr_256[alpha][ni][nk] &= l_mask;
				}
			}
		}

		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			int mj = i_prime(mi, alpha);
			VECTOR L_tr_J_256[rank][vtl_len] = {0};
			// use vtl (2 ni x vtl_len)  **PS. If rank is odd, the last value is handled separately.**
			for (int ni = 0; ni < rank_floor2; ni += 2) {
				for (int nj = 0; nj < rank; ++nj) {
					VECTOR k1_lh = mtk2_16[get_gf16m(pkx->map1.Qalpha1[mi][alpha], ni, nj) ^
					                       (get_gf16m(pkx->map1.Qalpha1[mi][alpha], ni + 1, nj) << 4)];
					VECTOR k2_lh = mtk2_16[get_gf16m(pkx->map1.Balpha[mi][alpha], nj, ni) ^
					                       (get_gf16m(pkx->map1.Balpha[mi][alpha], nj, ni + 1) << 4)];
					for (int nk = 0; nk < vtl_len; ++nk) {
						L_tr_J_256[ni][nk] ^= VEC_SHUFFLE(k1_lh, AxS_tr_256[alpha][nj][nk]);
						R_tr_J_256[mj][alpha][ni][nk] ^= VEC_SHUFFLE(k2_lh, Q2xS_tr_256[alpha][nj][nk]);
					}
				}
			}

#if rank % 2  //  rank is odd, use vl no vtl(vectorized look-up),
			for (int nj = 0; nj < rank; ++nj) {
				VECTOR k1_lh = mtk2_16[get_gf16m(pkx->map1.Qalpha1[mi][alpha], rank_floor2, nj)];
				VECTOR k2_lh = mtk2_16[get_gf16m(pkx->map1.Balpha[mi][alpha], nj, rank_floor2)];
				for (int nk = 0; nk < vtl_len; ++nk) {
					L_tr_J_256[rank_floor2][nk] ^= VEC_SHUFFLE(k1_lh, AxS_tr_256[alpha][nj][nk]);
					R_tr_J_256[mj][alpha][rank_floor2][nk] ^= VEC_SHUFFLE(k2_lh, Q2xS_tr_256[alpha][nj][nk]);
				}
			}
#endif

			// nibble splite
			for (int ni = 0; ni < rank_floor2; ni += 2) {
				for (int nk = 0; nk < vtl_len; ++nk) {
					L_tr_J_256[ni + 1][nk] = (L_tr_J_256[ni][nk] >> 4) & l_mask;
					L_tr_J_256[ni][nk] &= l_mask;

					R_tr_J_256[mj][alpha][ni + 1][nk] = (R_tr_J_256[mj][alpha][ni][nk] >> 4) & l_mask;
					R_tr_J_256[mj][alpha][ni][nk] &= l_mask;
				}
			}

			alignas(VEC_LENGTH) uint8_t L_J[rank_next2][vtl_len * VEC_LENGTH] = {0};
			jogressMatrixTr_vector((VECTOR*)L_J, (VECTOR*)L_tr_J_256, 1, n_SNOVA);

			// nibble L_J
			for (int ni = 0; ni < rank; ni += 2) {
				VECTOR* L_J_l_256 = (VECTOR*)L_J[ni];
				VECTOR* L_J_h_256 = (VECTOR*)L_J[ni + 1];
				VECTOR* L_J_256 = L_J_nibble256 + ((mj * alpha_SNOVA + alpha) * rank_next32 + ni / 2) * vtl_len;
				for (int nk = 0; nk < vtl_len; ++nk) {
					L_J_256[nk] = L_J_l_256[nk] ^ (L_J_h_256[nk] << 4);
				}
			}
		}
	}

	VECTOR h_256[m_SNOVA][rank][rank] = {0};
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		// clac P
		gf16m_t P[n_SNOVA][n_SNOVA] = {0};
		VECTOR P_J_256[n_SNOVA * rank][vtl_len] = {0};
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
			memcpy((uint8_t*)P[ni], (uint8_t*)pkx->map1.P11[mi][ni], sizeof(pkx->map1.P11[mi][ni]));
			memcpy((uint8_t*)(P[ni] + v_SNOVA), (uint8_t*)pkx->map1.P12[mi][ni], sizeof(pkx->map1.P11[mi][ni]));
		}
		for (int ni = v_SNOVA; ni < n_SNOVA; ++ni) {
			memcpy((uint8_t*)P[ni], (uint8_t*)pkx->map1.P21[mi][ni - v_SNOVA], sizeof(pkx->map1.P21[mi][ni - v_SNOVA]));
			memcpy((uint8_t*)(P[ni] + v_SNOVA), (uint8_t*)pkx->P22[mi][ni - v_SNOVA], sizeof(pkx->P22[mi][ni - v_SNOVA]));
		}
		jogressMatrix_vector((VECTOR*)P_J_256, (uint8_t*)P, n_SNOVA, n_SNOVA);

		// evaluation start!!
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			int mi_prime = i_prime_inv(mi, alpha);
			VECTOR LJxPJ_256[rank][vtl_len] = {0};
			VECTOR LJxPJ_256_nibble[rank_next32][vtl_len] = {0};

			// LJ x PJ, Main loop (n_SNOVA^2 * rank^3 times)
			for (int nk = 0; nk < n_SNOVA * rank; ++nk) {
				for (int ni = 0; ni < rank_next32; ++ni) {
					VECTOR k_lh = mtk2_16[VEC_BYTE(
					                          L_J_nibble256, ((mi * alpha_SNOVA + alpha) * rank_next32 + ni) * vtl_len * VEC_LENGTH + nk)];
					for (int nj = 0; nj < vtl_len; ++nj) {
						LJxPJ_256_nibble[ni][nj] ^= VEC_SHUFFLE(k_lh, P_J_256[nk][nj]);
					}
				}
			}

			// nibble splite, **PS. If rank is odd, the last value is handled separately.**
			for (int ni = 0; ni < rank_floor2; ni += 2) {
				VECTOR* LJxPJ_256_nibble_l = LJxPJ_256[ni];
				VECTOR* LJxPJ_256_nibble_h = LJxPJ_256[ni + 1];
				for (int nj = 0; nj < vtl_len; ++nj) {
					LJxPJ_256_nibble_l[nj] = LJxPJ_256_nibble[ni / 2][nj] & l_mask;
					LJxPJ_256_nibble_h[nj] = ((LJxPJ_256_nibble[ni / 2][nj] >> 4) & l_mask);
				}
			}

#if rank % 2  // rank is odd
			for (int nj = 0; nj < vtl_len; ++nj) {
				LJxPJ_256[rank_floor2][nj] = LJxPJ_256_nibble[rank_floor2 / 2][nj];
			}
#endif

			// (LJ x PJ) x RJ, Secondary loop (n_SNOVA * rank^3 times)
			for (int nk = 0; nk < vtl_len; ++nk) {
				for (int ni = 0; ni < rank; ++ni) {
					for (int nj = 0; nj < rank; ++nj) {
						gf16_32_mul_32_add((LJxPJ_256[ni] + nk), (R_tr_J_256[mi][alpha][nj] + nk), (&h_256[mi_prime][ni][nj]));
					}
				}
			}
		}
	}

	// set hash_in_GF16Matrix
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		for (int ni = 0; ni < rank; ++ni) {
			for (int nj = 0; nj < rank; ++nj) {
				set_gf16m(hash_in_GF16Matrix[mi], ni, nj, horizontal_xor_256(h_256[mi][ni][nj]));
			}
		}
	}
}

static int verify_signture_vtl_core(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature,
                                    const public_key_expand* pkx) {
	uint8_t hash_in_bytes[bytes_hash];
	uint8_t signed_hash[bytes_hash];
	const uint8_t *pt_salt = pt_signature + bytes_signature;

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

	if (convert_bytes_to_GF16s(pt_signature, (gf16_t * )signature_in_GF16Matrix, GF16s_signature)) {
		return -1;
	}

	evaluation_vector_vtl(hash_in_GF16Matrix, pkx, signature_in_GF16Matrix);
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

static inline int verify_signture_pkx_vtl(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature,
        const uint8_t *pkx_pck) {
	return verify_signture_vtl_core(pt_digest, bytes_digest, pt_signature, (public_key_expand*)pkx_pck);
}

static void snova_plasma_init(void) {
	static int first_plasma_time = 1;
	if (first_plasma_time) {
		first_plasma_time = 0;
		init_vector_table();
	}
}

#if FIXED_ABQ
static void gen_fixedABQ(const char* abq_seed) {
	uint8_t rng_out[m_SNOVA * alpha_SNOVA * (lsq_SNOVA + l_SNOVA)];
	uint8_t q12[2 * m_SNOVA * alpha_SNOVA * l_SNOVA];

	shake256(rng_out, m_SNOVA * alpha_SNOVA * (lsq_SNOVA + l_SNOVA), (uint8_t*)abq_seed, strlen(abq_seed));

	convert_bytes_to_GF16s(rng_out, fixed_abq, 2 * m_SNOVA * alpha_SNOVA * lsq_SNOVA);
	convert_bytes_to_GF16s(&rng_out[m_SNOVA * alpha_SNOVA * lsq_SNOVA], q12, 2 * m_SNOVA * alpha_SNOVA * l_SNOVA);

	for (int pi = 0; pi < m_SNOVA; ++pi) {
		for (int index = 0; index < alpha_SNOVA; ++index) {
			be_invertible_by_add_aS(&fixed_abq[(pi * alpha_SNOVA + index) * lsq_SNOVA]);
		}
		for (int index = 0; index < alpha_SNOVA; ++index) {
			be_invertible_by_add_aS(&fixed_abq[((m_SNOVA + pi) * alpha_SNOVA + index) * lsq_SNOVA]);
		}
		for (int index = 0; index < alpha_SNOVA; ++index) {
			gen_a_FqS(&q12[(pi * alpha_SNOVA + index) * l_SNOVA],
			          &fixed_abq[((2 * m_SNOVA + pi) * alpha_SNOVA + index) * lsq_SNOVA]);
		}
		for (int index = 0; index < alpha_SNOVA; ++index) {
			gen_a_FqS(&q12[((m_SNOVA + pi) * alpha_SNOVA + index) * l_SNOVA],
			          &fixed_abq[((3 * m_SNOVA + pi) * alpha_SNOVA + index) * lsq_SNOVA]);
		}
	}
}
#endif

/**
 * SNOVA init
 */
static void snova_init(void) {
	static int first_time = 1;
	if (first_time) {
		first_time = 0;
		init_gf16_tables();
		gen_S_array();
#if FIXED_ABQ
		gen_fixedABQ("SNOVA_ABQ");
#endif
		snova_plasma_init();
	}
}

static void generate_keys_core(snova_key_elems* key_elems, const uint8_t* pk_seed, const uint8_t* sk_seed) {
	gen_seeds_and_T12(key_elems->T12, sk_seed);
	memcpy(key_elems->pk.pt_public_key_seed, pk_seed, seed_length_public);
	gen_ABQ(&(key_elems->map1), pk_seed);
	gen_F_vtl(&(key_elems->map2), &(key_elems->map1), key_elems->T12);
	gen_P22_vtl(key_elems->pk.P22, key_elems->T12, key_elems->map1.P21, key_elems->map2.F12);
}

int SNOVA_NAMESPACE(sk_expand)(expanded_SK* skx, const uint8_t* sk) {
	const uint8_t *pk_seed = sk;
	const uint8_t *sk_seed = sk + seed_length_public;
	snova_key_elems key_elems;

	generate_keys_core(&key_elems, pk_seed, sk_seed);
	sk_pack((uint8_t*)skx, &key_elems, sk_seed);

	return 0;
}

int SNOVA_NAMESPACE(genkeys)(uint8_t* pk, uint8_t* ssk, const uint8_t* seed) {
	const uint8_t *pkseed = seed;
	const uint8_t *skseed = seed + seed_length_public;

	snova_init();
	snova_key_elems key_elems;
	generate_keys_core(&key_elems, pkseed, skseed);
	pk_pack(pk, &key_elems);
	memcpy(ssk, pkseed, seed_length_public);
	memcpy(ssk + seed_length_public, skseed, seed_length_private);

	return 0;
}

int SNOVA_NAMESPACE(sign)(const expanded_SK* esk, uint8_t* pt_signature, const uint8_t* digest, const size_t bytes_digest,
                          const uint8_t *array_salt) {
	snova_init();
	sk_gf16 sk_upk;
	sk_unpack(&sk_upk, (uint8_t*)esk);
	int res = sign_digest_core_gnl_vtl(pt_signature, digest, bytes_digest, (uint8_t*)array_salt, sk_upk.Aalpha, sk_upk.Balpha,
	                                   sk_upk.Qalpha1, sk_upk.Qalpha2, sk_upk.T12, sk_upk.F11, sk_upk.F12, sk_upk.F21,
	                                   sk_upk.pt_public_key_seed, sk_upk.pt_private_key_seed);
	return res;
}

int SNOVA_NAMESPACE(pk_expand)(expanded_PK* pkx, const uint8_t* pk) {
	snova_init();
	expand_public_pack_core((uint8_t*)pkx, pk);
	return 0;
}

int SNOVA_NAMESPACE(verify)(const expanded_PK* pk, const uint8_t* pt_signature, const uint8_t* pt_digest,
                            const size_t bytes_digest) {
	snova_init();
	return verify_signture_vtl_core(pt_digest, bytes_digest, pt_signature, (public_key_expand*)pk);
}
