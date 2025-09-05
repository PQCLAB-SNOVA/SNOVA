// SPDX-License-Identifier: MIT

/**
 * Optimized implementation using AVX2 vectorization
 *
 * Copyright (c) 2025 SNOVA TEAM
 */

#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>
#include <string.h>

extern uint64_t cycles0, cycles1, cycles2, cycles3, cycles4, cycles5, cycles6, cycles7;
uint64_t get_cycles(void);

#include "snova.h"
#include "symmetric.h"

#if SNOVA_q == 16
#error "SNOVA_q == 16"
#include "stop"
#endif

#if SNOVA_l != 4
#error "SNOVA_l != 4"
#include "stop"
#endif

typedef uint8_t gf_t;

/**
 * Constant time functions. CT is according to valgrind
 */
static inline uint16_t ct_is_not_zero(uint16_t val) {
	// return (val | (val >> 1) | (val >> 2) | (val >> 3) | (val >> 4)) & 1;
	return val != 0;
}

static inline int ct_is_negative(int val) {
	// return ((val >> 31) & 1);
	return val < 0;
}

/**
 * Constant time GF(q) inverse
 *
 * Use that x^q = x and therefore x^(q-2) = x^-1
 */
static uint16_t ct_gf_inverse(uint16_t val) {
	uint32_t res = val;
	uint32_t pow = val * val;
	uint8_t bits = (SNOVA_q - 2) / 2;

	for (int j1 = 0; j1 < 3; j1++) {
		if (bits & 1) {
			res = (pow * res) % SNOVA_q;
		}
		pow = (pow * pow) % SNOVA_q;
		bits = bits >> 1;
	}
	if (bits & 1) {
		res = (pow * res) % SNOVA_q;
	}

	return res % SNOVA_q;
}

/**
 * Initialization
 */

#define gf_S SNOVA_NAMESPACE(Smat)

#if SNOVA_l == 4 && SNOVA_q == 23
uint16_t gf_S[SNOVA_l * SNOVA_l2] = {1, 0, 0, 0, 0, 1,  0,  0,  0,  0,  1,  0, 0,  0,  0, 1,  1,  2, 3,  0, 2,  3,
                                     0, 1, 3, 0, 1, 2,  0,  1,  2,  22, 14, 8, 6,  8,  8, 14, 8,  2, 6,  8, 14, 0,
                                     8, 2, 0, 6, 2, 14, 18, 12, 14, 14, 13, 5, 18, 13, 9, 13, 12, 5, 13, 19
                                    };
#define SNOVA_INIT

#elif SNOVA_l == 4 && SNOVA_q == 19
static uint16_t gf_S[SNOVA_l * SNOVA_l2] = {1, 0,  0,  0, 0,  1, 0, 0, 0, 0,  1,  0,  0, 0,  0, 1,  1, 2,  3, 0, 2,  3,
                                            0, 1,  3,  0, 1,  2, 0, 1, 2, 15, 14, 8,  6, 8,  8, 14, 8, 18, 6, 8, 14, 13,
                                            8, 18, 13, 2, 10, 3, 7, 7, 3, 0,  11, 15, 7, 11, 1, 3,  7, 15, 3, 17
                                           };
#define SNOVA_INIT

#elif SNOVA_l == 4 && SNOVA_q == 11
static uint16_t gf_S[SNOVA_l * SNOVA_l2] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 2, 3, 1, 2,
                                            3, 0, 2, 3, 0, 1, 3, 0, 1, 6, 3, 8, 6, 9, 8, 3, 8, 6, 6, 8, 3, 1,
                                            9, 6, 1, 2, 3, 4, 6, 3, 4, 5, 9, 2, 6, 9, 4, 5, 3, 2, 5, 7
                                           };
#define SNOVA_INIT

#else
uint16_t gf_S[SNOVA_l * SNOVA_l2] = {0};

static void gen_S_array(void) {
	memset(gf_S, 0, sizeof(gf_S));

	for (int i1 = 0; i1 < SNOVA_l; i1++) {
		gf_S[i1 * SNOVA_l + i1] = 1;
	}

	// Set S^1, the irreducible S matrix
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_S[SNOVA_l2 + i1 * SNOVA_l + j1] = ((Q_A + i1 + j1) & Q_B) % SNOVA_q;
		}
	gf_S[2 * SNOVA_l2 - 1] = Q_C % SNOVA_q;

	for (int si = 2; si < SNOVA_l; si++) {
		for (int i1 = 0; i1 < SNOVA_l; i1++)
			for (int j1 = 0; j1 < SNOVA_l; j1++) {
				uint16_t sum = 0;
				for (int k1 = 0; k1 < SNOVA_l; k1++) {
					sum += gf_S[SNOVA_l2 + i1 * SNOVA_l + k1] * gf_S[(si - 1) * SNOVA_l2 + k1 * SNOVA_l + j1] % SNOVA_q;
				}
				gf_S[si * SNOVA_l2 + i1 * SNOVA_l + j1] = sum % SNOVA_q;
			}
	}
}

static int first_time = 1;

#define SNOVA_INIT      \
    if (first_time) {   \
        first_time = 0; \
        gen_S_array();  \
    }
#endif

/**
 * Utilities
 */
static gf_t gf_mat_det(gf_t *a) {
#define DET_SUB(a, b) (a - b)
#define DET_MULT(a, b) (a * b)
	int32_t det = 0;
#if SNOVA_l == 2
	det = DET_SUB(DET_MULT(a[0], a[3]), DET_MULT(a[1], a[2]));
#elif SNOVA_l == 3
	det = DET_MULT(a[0], DET_SUB(DET_MULT(a[4], a[8]), DET_MULT(a[5], a[7])));
	det += DET_MULT(a[1], DET_SUB(DET_MULT(a[5], a[6]), DET_MULT(a[3], a[8])));
	det += DET_MULT(a[2], DET_SUB(DET_MULT(a[3], a[7]), DET_MULT(a[4], a[6])));
#elif SNOVA_l == 4
	int32_t DET_l;
	int32_t DET_r;
#define DET_L(x, y) DET_l = DET_SUB(DET_MULT(a[x], a[4 + y]), DET_MULT(a[y], a[4 + x]))
#define DET_R(x, y) DET_r = DET_SUB(DET_MULT(a[8 + x], a[12 + y]), DET_MULT(a[8 + y], a[12 + x]))
#define DET22(x1, y1, x2, y2) \
    DET_L(x1, y1);            \
    DET_R(x2, y2);            \
    det += DET_MULT(DET_l, DET_r)
	DET22(0, 1, 2, 3);
	DET22(0, 2, 3, 1);
	DET22(0, 3, 1, 2);
	DET22(1, 2, 0, 3);
	DET22(1, 3, 2, 0);
	DET22(2, 3, 0, 1);
#undef DET_R
#undef DET22
#undef DET_L
#elif SNOVA_l == 5
	int32_t DET_l;
	int32_t DET_r;
#define DET_L(x, y) DET_l = DET_SUB(DET_MULT(a[x], a[5 + y]), DET_MULT(a[y], a[5 + x]))
#define DET_R2(x, y, z) DET_MULT(DET_SUB(DET_MULT(a[10 + x], a[15 + y]), DET_MULT(a[10 + y], a[15 + x])), a[20 + z])
#define DET_R3(x, y, z) DET_r = DET_R2(x, y, z) + DET_R2(y, z, x) + DET_R2(z, x, y)
#define DET23(x1, y1, x2, y2, z2) \
    DET_L(x1, y1);                \
    DET_R3(x2, y2, z2);           \
    det += DET_MULT(DET_l, DET_r)
	DET23(0, 1, 2, 3, 4);
	DET23(0, 2, 3, 1, 4);
	DET23(0, 3, 1, 2, 4);
	DET23(0, 4, 1, 3, 2);
	DET23(1, 2, 0, 3, 4);
	DET23(1, 3, 2, 0, 4);
	DET23(1, 4, 2, 3, 0);
	DET23(2, 3, 0, 1, 4);
	DET23(2, 4, 0, 3, 1);
	DET23(3, 4, 2, 0, 1);
#undef DET_R2
#undef DET_R3
#undef DET23
#undef DET_L
#else
#error "Unsupported rank"
#endif
#undef DET_SUB
#undef DET_MULT
	return det % SNOVA_q;
}

static void convert_bytes_to_GF(gf_t *gf_array, const uint8_t *byte_array, size_t num) {
	for (size_t idx = 0; idx < num; idx++) {
		gf_array[idx] = byte_array[idx] % SNOVA_q;
	}
}

// Used to compress PK (genkey) and SIG(sign)
static void compress_gf(uint8_t *byte_array, const gf_t *gf_array, size_t num) {
	size_t idx = 0;
	size_t out_idx = 0;
	size_t num_bytes = BYTES_GF(num);

	do {
		uint64_t val = 0;
		uint64_t fact = 1;

		int i1 = 0;
		while (i1 < PACK_GF && idx < num) {
			val += fact * (gf_array[idx] % SNOVA_q);
			idx++;
			i1++;
			fact *= SNOVA_q;
		}

		i1 = (i1 + 1) / 2;
		int j1 = 0;
		while (j1 < PACK_BYTES && out_idx < num_bytes) {
			byte_array[out_idx] = val & 0xff;
			out_idx++;
			val = val >> 8;
			j1++;
		}
	} while (idx < num);
}

// Used to expand PK(verify) and SIG(verify)
static int expand_gf(gf_t *gf_array, const uint8_t *byte_array, size_t num) {
	size_t num_bytes = BYTES_GF(num);
	size_t idx = 0;
	size_t out_idx = 0;
	uint64_t val;
	int res = 0;

	do {
		val = 0;

		int i1 = 0;
		while (i1 < PACK_BYTES && idx < num_bytes) {
			val = val ^ ((uint64_t)(byte_array[idx]) << (8 * i1));
			idx++;
			i1++;
		}

		int j1 = 0;
		while (j1 < PACK_GF && out_idx < num) {
			gf_array[out_idx] = val % SNOVA_q;
			val = val / SNOVA_q;
			out_idx++;
			j1++;
		}

		res |= val;
	} while (out_idx < num);

	return res;
}

// Used to compress PK (genkey)
static void compress_pk(uint8_t *pk, const gf_t *P22) {
#ifdef SYMMETRIC
	gf_t P22c[NUMGF_PK] = {0};
	gf_t *curval = &P22c[0];

	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int ni = 0; ni < SNOVA_o; ++ni) {
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				for (int j1 = i1; j1 < SNOVA_l; j1++) {
					*curval = P22[((mi * SNOVA_o + ni) * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + j1];
					curval++;
				}

				for (int nj = ni + 1; nj < SNOVA_o; ++nj)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						*curval = P22[((mi * SNOVA_o + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1];
						curval++;
					}
			}
		}

	compress_gf(pk, P22c, NUMGF_PK);
#else

	compress_gf(pk, P22, NUMGF_PK);
#endif
}

/**
 * Expand the public key from a seed. Make symmetric
 */
static void expand_public(gf_t *P_matrix, const uint8_t *seed) {
	uint8_t pk_bytes[NUM_GEN_PUB_BYTES];

	snova_pk_expander_t instance;
	snova_pk_expander_init(&instance, seed, SEED_LENGTH_PUBLIC);
	snova_pk_expander(pk_bytes, NUM_GEN_PUB_BYTES, &instance);

#ifdef SYMMETRIC
	gf_t pk_gf[NUM_GEN_PUB_GF];
	convert_bytes_to_GF(pk_gf, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);

	// Make symmetric
	gf_t *P11 = P_matrix;
	gf_t *P12 = P_matrix + SNOVA_m * SNOVA_v * SNOVA_v * SNOVA_l2;
	gf_t *P21 = P_matrix + SNOVA_m * SNOVA_v * SNOVA_n * SNOVA_l2;
	gf_t *abq = P21 + SNOVA_m * SNOVA_o * SNOVA_v * SNOVA_l2;

	gf_t *curval = &pk_gf[0];

	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int ni = 0; ni < SNOVA_v; ++ni) {
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = i1; j1 < SNOVA_l; j1++) {
					P11[((mi * SNOVA_v + ni) * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
					curval++;
				}

			for (int nj = ni + 1; nj < SNOVA_v; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}

			for (int nj = 0; nj < SNOVA_o; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P12[((mi * SNOVA_v + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}
		}

	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int ni = 0; ni < SNOVA_v; ++ni)
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = i1; j1 < SNOVA_l; j1++) {
					P11[((mi * SNOVA_v + ni) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] =
					    P11[((mi * SNOVA_v + ni) * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + j1];
				}

	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int ni = 0; ni < SNOVA_v; ++ni)
			for (int nj = ni + 1; nj < SNOVA_v; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P11[((mi * SNOVA_v + nj) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] =
						    P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1];
					}

	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int ni = 0; ni < SNOVA_v; ++ni)
			for (int nj = 0; nj < SNOVA_o; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						P21[((mi * SNOVA_o + nj) * SNOVA_v + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] =
						    P12[((mi * SNOVA_v + ni) * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + j1];
					}

	for (int idx = 0; idx < SNOVA_m * SNOVA_alpha * 2 * (SNOVA_l2 + SNOVA_l); idx++) {
		abq[idx] = curval[idx];
	}
#else

	convert_bytes_to_GF(P_matrix, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);
#endif
}

static void hash_combined(uint8_t *hash_out, const uint8_t *digest, size_t len_digest, const uint8_t *pk_seed,
                          const uint8_t *salt) {
	shake_t state;
	shake256_init(&state);
	shake_absorb(&state, pk_seed, SEED_LENGTH_PUBLIC);
	shake_absorb(&state, digest, len_digest);
	shake_absorb(&state, salt, BYTES_SALT);
	shake_finalize(&state);
	shake_squeeze(hash_out, BYTES_HASH, &state);
}

/**
 * Improve q and calculate Q matrix
 *
 * Must be constant time as expand_T12 uses gen_a_FqS.
 */
static inline void gen_a_FqS(uint16_t *Qm, gf_t *q) {
	int16_t not_zero = -ct_is_not_zero(q[SNOVA_l - 1]);

	q[SNOVA_l - 1] = (not_zero & q[SNOVA_l - 1]) | ((not_zero ^ -1) & (SNOVA_q - (q[0] + 1 - ct_is_not_zero(q[0]))));

	for (int i1 = 0; i1 < SNOVA_l2; i1++) {
		uint16_t sum = 0;
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			sum += q[j1] * gf_S[j1 * SNOVA_l2 + i1];
		}
		Qm[i1] = sum % SNOVA_q;
	}
}

/**
 * Expand T12 matrix and coefficients. Shared by genkey and sign
 */
#define REJECTION_LIMIT ((256 / SNOVA_q) * SNOVA_q)
#define SK_BLOCK_SIZE 32
static void expand_T12(uint16_t *T12, const uint8_t *seed) {
	gf_t T12coef[SNOVA_o * SNOVA_v * SNOVA_l];
	gf_t sk_data[SK_BLOCK_SIZE];
	shake_t state;

	shake256_init(&state);
	shake_absorb(&state, seed, SEED_LENGTH_PRIVATE);
	shake_finalize(&state);

	size_t idx = SK_BLOCK_SIZE;
	size_t t_idx = 0;

	while (t_idx < SNOVA_o * SNOVA_v * SNOVA_l) {
		if (idx >= SK_BLOCK_SIZE) {
			shake_squeeze_keep(sk_data, SK_BLOCK_SIZE, &state);
			idx = 0;
		}

		// Rejection sampling
		T12coef[t_idx] = sk_data[idx] % SNOVA_q;
		t_idx += ct_is_negative((int)sk_data[idx] - REJECTION_LIMIT);

		idx++;
	}
	shake_release(&state);

	for (int i1 = 0; i1 < SNOVA_o * SNOVA_v; i1++) {
		gen_a_FqS(&T12[i1 * SNOVA_l2], &T12coef[i1 * SNOVA_l]);
	}
}

/**
 * Ensure that a matrix is invertible by adding multiples of S
 */
static void be_invertible_by_add_aS(gf_t *mat, const gf_t *orig) {
	for (int i1 = 0; i1 < SNOVA_l2; i1++) {
		mat[i1] = orig[i1];
	}

	if (gf_mat_det(mat) == 0) {
		for (gf_t f1 = 1; f1 < SNOVA_q; f1++) {
			for (int i1 = 0; i1 < SNOVA_l2; i1++) {
				mat[i1] = (mat[i1] + (f1 * gf_S[SNOVA_l2 + i1])) % SNOVA_q;
			}
			if (gf_mat_det(mat) != 0) {
				break;
			}
		}
	}
}

/**
 * Optimized version of genkey.
 */
int SNOVA_NAMESPACE(genkeys)(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
	SNOVA_INIT

	/**
	 * Gen T12 matrix
	 */
	uint16_t T12[SNOVA_o * SNOVA_v * SNOVA_l2];
	expand_T12(T12, seed + SEED_LENGTH_PUBLIC);

	/**
	 * Gen Public matrix but not ABQ
	 */
	gf_t P_matrix[NUM_PUB_GF];

	expand_public(P_matrix, seed);

	/**
	 * Calculate F12 matrix, use P11
	 */
	gf_t *P11gf = P_matrix;
	gf_t *P12 = P_matrix + SNOVA_o * SNOVA_v * SNOVA_v * SNOVA_l2;
	gf_t *P21 = P_matrix + SNOVA_o * SNOVA_v * SNOVA_n * SNOVA_l2;

	uint16_t P11[SNOVA_l2];
	gf_t P22gf[SNOVA_o * SNOVA_o * SNOVA_o * SNOVA_l2] = {0};

	for (int mi = 0; mi < SNOVA_m; mi++) {
		uint16_t F12[SNOVA_v * SNOVA_o * SNOVA_l2] = {0};
		uint16_t P22[SNOVA_o * SNOVA_o * SNOVA_l2] = {0};

		for (int ni = 0; ni < SNOVA_v; ni++)
			for (int nj = 0; nj < SNOVA_v; nj++) {
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					P11[i1] = P11gf[((mi * SNOVA_v + nj) * SNOVA_v + ni) * SNOVA_l2 + i1];
				}
				for (int nk = 0; nk < SNOVA_o; nk++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								F12[(nj * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    P11[i1 * SNOVA_l + k1] * T12[(ni * SNOVA_o + nk) * SNOVA_l2 + k1 * SNOVA_l + j1];
			}

		// Use P12
		for (int i1 = 0; i1 < SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
			F12[i1] = (F12[i1] + P12[mi * SNOVA_v * SNOVA_o * SNOVA_l2 + i1]) % SNOVA_q;
		}

		// Establish P22, first step
		for (int ni = 0; ni < SNOVA_v; ni++)
			for (int nk = 0; nk < SNOVA_o; nk++)
				for (int nj = 0; nj < SNOVA_o; nj++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								P22[(nj * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    T12[(ni * SNOVA_o + nj) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    F12[(ni * SNOVA_o + nk) * SNOVA_l2 + k1 * SNOVA_l + j1];

		for (int i1 = 0; i1 < SNOVA_o * SNOVA_o * SNOVA_l2; i1++) {
			P22[i1] = P22[i1] % SNOVA_q;
		}

		/**
		 * Calculate P22. Uses P21
		 */
		for (int ni = 0; ni < SNOVA_v; ni++)
			for (int nj = 0; nj < SNOVA_o; nj++)
				for (int nk = 0; nk < SNOVA_o; nk++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								P22[(nj * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    P21[((mi * SNOVA_o + nj) * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    T12[(ni * SNOVA_o + nk) * SNOVA_l2 + k1 * SNOVA_l + j1];

		// Negate P22
		for (int i1 = 0; i1 < SNOVA_o * SNOVA_o * SNOVA_l2; i1++) {
			P22gf[mi * SNOVA_o * SNOVA_o * SNOVA_l2 + i1] = (SNOVA_q - (P22[i1] % SNOVA_q)) % SNOVA_q;
		}
	}

	/**
	 * Output public and secret keys
	 */
	memcpy(pk, seed, SEED_LENGTH_PUBLIC);
	compress_pk(pk + SEED_LENGTH_PUBLIC, P22gf);
	memcpy(sk, seed, SEED_LENGTH);

	return 0;
}

/**
 * SK expansion.
 */
int SNOVA_NAMESPACE(sk_expand)(expanded_SK *skx, const uint8_t *sk) {
	SNOVA_INIT

	memcpy(skx->sk_seed, sk, SEED_LENGTH);
	expand_T12(skx->T12, sk + SEED_LENGTH_PUBLIC);

	gf_t P_matrix[NUM_PUB_GF];
	expand_public(P_matrix, skx->sk_seed);

	/**
	 * Calculate F12, F21
	 */
	gf_t *P21 = P_matrix + SNOVA_m * SNOVA_v * SNOVA_n * SNOVA_l2;

	for (int i1 = 0; i1 < SNOVA_m * SNOVA_v * SNOVA_v * SNOVA_l2; i1++) {
		skx->P11[i1] = P_matrix[i1];
	}

	alignas(32) uint16_t F21[SNOVA_m * SNOVA_o * SNOVA_v * SNOVA_l2] = {0};

	for (int nk = 0; nk < SNOVA_v; nk++)
		for (int mi = 0; mi < SNOVA_m; mi++)
			for (int ni = 0; ni < SNOVA_o; ni++)
				for (int nj = 0; nj < SNOVA_v; nj++) {
					__m256i p11_256[4];
					uint16_t *p11_16 = (uint16_t *)p11_256;
					__m256i t12_256[4];
					uint16_t *t12_16 = (uint16_t *)t12_256;
					__m256i *pF21 = (__m256i *)&F21[((mi * SNOVA_o + ni) * SNOVA_v + nj) * SNOVA_l2];

					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								p11_16[k1 * SNOVA_l2 + i1 * SNOVA_l + j1] =
								    skx->P11[((mi * SNOVA_v + nk) * SNOVA_v + nj) * SNOVA_l2 + k1 * SNOVA_l + j1];

					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								t12_16[k1 * SNOVA_l2 + i1 * SNOVA_l + j1] =
								    skx->T12[(nk * SNOVA_o + ni) * SNOVA_l2 + i1 * SNOVA_l + k1];

					for (int k1 = 0; k1 < SNOVA_l; k1++) {
						*pF21 += _mm256_mullo_epi16(p11_256[k1], t12_256[k1]);
					}
				}

	for (int i1 = 0; i1 < SNOVA_m * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
		skx->F21[i1] = (F21[i1] + P21[i1]) % SNOVA_q;
	}

#ifndef SYMMETRIC
	gf_t *P12 = P_matrix + SNOVA_m * SNOVA_v * SNOVA_v * SNOVA_l2;
	uint16_t F12[SNOVA_m * SNOVA_o * SNOVA_v * SNOVA_l2] = {0};

	for (int mi = 0; mi < SNOVA_m; mi++)
		for (int nj = 0; nj < SNOVA_v; nj++)
			for (int ni = 0; ni < SNOVA_v; ni++)
				for (int nk = 0; nk < SNOVA_o; nk++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								F12[((mi * SNOVA_v + nj) * SNOVA_o + nk) * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    skx->P11[((mi * SNOVA_v + nj) * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    skx->T12[(ni * SNOVA_o + nk) * SNOVA_l2 + k1 * SNOVA_l + j1];

	// Use P12
	for (int i1 = 0; i1 < SNOVA_o * SNOVA_v * SNOVA_o * SNOVA_l2; i1++) {
		skx->F12[i1] = (F12[i1] + P12[i1]) % SNOVA_q;
	}
#endif

	// Generate ABQ, fix q
	gf_t Am[SNOVA_m * SNOVA_alpha * SNOVA_l2];
	gf_t Bm[SNOVA_m * SNOVA_alpha * SNOVA_l2];

	gf_t *aptr = P_matrix + (SNOVA_o * (SNOVA_n * SNOVA_n - SNOVA_o * SNOVA_o)) * SNOVA_l2;
	gf_t *bptr = aptr + SNOVA_m * SNOVA_alpha * SNOVA_l2;
	gf_t *q1 = aptr + 2 * SNOVA_o * SNOVA_alpha * SNOVA_l2;
	gf_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	for (int idx = 0; idx < SNOVA_m * SNOVA_alpha; idx++) {
		be_invertible_by_add_aS(&Am[idx * SNOVA_l2], &aptr[idx * SNOVA_l2]);
		be_invertible_by_add_aS(&Bm[idx * SNOVA_l2], &bptr[idx * SNOVA_l2]);
		gen_a_FqS(&(skx->Q1[idx * SNOVA_l2]), &q1[idx * SNOVA_l]);
		gen_a_FqS(&(skx->Q2[idx * SNOVA_l2]), &q2[idx * SNOVA_l]);
	}

	for (int idx = 0; idx < SNOVA_m * SNOVA_alpha * SNOVA_l2; idx++) {
		skx->Am[idx] = Am[idx];
		skx->Bm[idx] = Bm[idx];
	}
	for (int idx = 0; idx < SNOVA_m * SNOVA_alpha * SNOVA_l; idx++) {
		skx->q1[idx] = q1[idx];
		skx->q2[idx] = q2[idx];
	}

	return 0;
}

#if __AVX2__
#define SNOVA_ml2 ((SNOVA_m + 1) * SNOVA_l2)
#else
#define SNOVA_ml2 (SNOVA_m * SNOVA_l2 + 1)
#endif

/**
 * Optimized version of Sign. Deterministic using the salt provided
 */
int SNOVA_NAMESPACE(sign)(const expanded_SK *skx, uint8_t *sig, const uint8_t *digest, size_t len_digest, const uint8_t *salt) {
	SNOVA_INIT

	// Calculate message has of size l^2o
	gf_t hash_in_GF16[GF16_HASH];

	uint8_t sign_hashb[BYTES_HASH];
	hash_combined(sign_hashb, digest, len_digest, skx->sk_seed, salt);
	expand_gf(hash_in_GF16, sign_hashb, GF16_HASH);

	// Find a solution for T.X
	alignas(32) uint16_t gauss[SNOVA_ml2][SNOVA_ml2];
	alignas(32) uint16_t gauss16[SNOVA_ml2][SNOVA_ml2];
	alignas(32) uint16_t solution[SNOVA_ml2] = {0};
	gf_t signature_in_GF[SNOVA_n * SNOVA_l2] = {0};

	int flag_redo;
	uint8_t num_sign = 0;

	do {
		memset(gauss, 0, sizeof(gauss));
		memset(gauss16, 0, sizeof(gauss16));

		num_sign++;
		if (num_sign == 255) {
			// Probability of getting here is about q^{-255}
			memset(sig, 0, BYTES_SIGNATURE);
			return -1;
		}
		flag_redo = 0;

		// generate the vinegar value
		uint8_t vinegar_in_byte[NUM_GEN_SEC_BYTES];
		shake_t v_instance;

		shake256_init(&v_instance);
		shake_absorb(&v_instance, skx->sk_seed + SEED_LENGTH_PUBLIC, SEED_LENGTH_PRIVATE);
		shake_absorb(&v_instance, sign_hashb, BYTES_HASH);
		shake_absorb(&v_instance, &num_sign, 1);
		shake_finalize(&v_instance);
		shake_squeeze(vinegar_in_byte, NUM_GEN_SEC_BYTES, &v_instance);

		expand_gf(signature_in_GF, vinegar_in_byte, SNOVA_v * SNOVA_l2);

		// Calculate Fvv
		alignas(32) uint16_t Fvv_in_GF16Matrix[SNOVA_m * SNOVA_l2] = {0};

		/**
		 * Whip signature
		 */
		alignas(32) uint16_t whipped_sig[SNOVA_l * SNOVA_v * SNOVA_l2] = {0};

		for (int ab = 0; ab < SNOVA_l; ++ab)
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++)
							whipped_sig[(ab * SNOVA_v + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] +=
							    gf_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1] * signature_in_GF[ni * SNOVA_l2 + k1 * SNOVA_l + j1];

		for (int i1 = 0; i1 < SNOVA_v * SNOVA_l * SNOVA_l2; i1++) {
			whipped_sig[i1] = whipped_sig[i1] % SNOVA_q;
		}

		/**
		 * Evaluate whipped central map
		 */
		alignas(32) uint16_t sum_t0[SNOVA_m * SNOVA_l * SNOVA_v * SNOVA_l2] = {0};
		alignas(32) uint16_t sum_t1[SNOVA_m * SNOVA_l2 * SNOVA_l2] = {0};

		// Right
		for (int mi = 0; mi < SNOVA_m; ++mi)
			for (int nj = 0; nj < SNOVA_v; ++nj)
				for (int ni = 0; ni < SNOVA_v; ++ni)
					for (int b1 = 0; b1 < SNOVA_l; ++b1) {
						__m256i p11_256[4];
						uint16_t *p11_16 = (uint16_t *)p11_256;
						__m256i ws_256[4];
						uint16_t *ws_16 = (uint16_t *)ws_256;
						__m256i *psum = (__m256i *)&sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_l2];

						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								for (int j1 = 0; j1 < SNOVA_l; j1++)
									p11_16[k1 * SNOVA_l2 + i1 * SNOVA_l + j1] =
									    skx->P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1];

						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								for (int j1 = 0; j1 < SNOVA_l; j1++)
									ws_16[k1 * SNOVA_l2 + i1 * SNOVA_l + j1] =
									    whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_l2 + k1 * SNOVA_l + j1];

						for (int k1 = 0; k1 < SNOVA_l; k1++) {
							*psum += _mm256_mullo_epi16(p11_256[k1], ws_256[k1]);
						}
					}

		for (int i1 = 0; i1 < SNOVA_m * SNOVA_v * SNOVA_l * SNOVA_l2; i1++) {
			sum_t0[i1] = sum_t0[i1] % SNOVA_q;
		}

		// Left, transposed whipped_sig
		for (int mi = 0; mi < SNOVA_m; ++mi)
			for (int ni = 0; ni < SNOVA_v; ++ni)
				for (int a1 = 0; a1 < SNOVA_l; ++a1)
					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								for (int k1 = 0; k1 < SNOVA_l; k1++)
									sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_l2 + i1 * SNOVA_l + j1] +=
									    whipped_sig[(a1 * SNOVA_v + ni) * SNOVA_l2 + k1 * SNOVA_l + i1] *
									    sum_t0[((mi * SNOVA_l + b1) * SNOVA_v + ni) * SNOVA_l2 + k1 * SNOVA_l + j1];

		for (int i1 = 0; i1 < SNOVA_m * SNOVA_l2 * SNOVA_l2; i1++) {
			sum_t1[i1] = sum_t1[i1] % SNOVA_q;
		}

		/**
		 * Apply A, B, q1 and q2, aka E matrix
		 */
		for (int mi = 0; mi < SNOVA_m; ++mi)
			for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
				int mi_prime = i_prime(mi, alpha);

				alignas(32) uint16_t gfm_temp1[SNOVA_l2] = {0};
				alignas(32) uint16_t gfm_temp2[SNOVA_l2] = {0};

				// apply q1 and q2
				for (int a1 = 0; a1 < SNOVA_l; ++a1) {
					alignas(32) uint16_t gfm_temp0[SNOVA_l2] = {0};

					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								gfm_temp0[i1 * SNOVA_l + j1] +=
								    sum_t1[(mi_prime * SNOVA_l2 + a1 * SNOVA_l + b1) * SNOVA_l2 + i1 * SNOVA_l + j1] *
								    skx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1];

					for (int i1 = 0; i1 < SNOVA_l2; i1++) {
						gfm_temp0[i1] = gfm_temp0[i1] % SNOVA_q;
					}

					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++) {
							gfm_temp1[i1 * SNOVA_l + j1] +=
							    gfm_temp0[i1 * SNOVA_l + j1] * skx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1];
						}
				}

				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					gfm_temp1[i1] = gfm_temp1[i1] % SNOVA_q;
				}

				// A and B
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++) {
							gfm_temp2[i1 * SNOVA_l + j1] += gfm_temp1[i1 * SNOVA_l + k1] *
							                                skx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + k1 * SNOVA_l + j1];
						}

				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					gfm_temp2[i1] = gfm_temp2[i1] % SNOVA_q;
				}

				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						for (int k1 = 0; k1 < SNOVA_l; k1++) {
							Fvv_in_GF16Matrix[mi * SNOVA_l2 + i1 * SNOVA_l + j1] +=
							    skx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + i1 * SNOVA_l + k1] *
							    gfm_temp2[k1 * SNOVA_l + j1];
						}

				// Set the last column of gauss matrix
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					gauss[mi * SNOVA_l2 + i1][SNOVA_o * SNOVA_l2] =
					    (hash_in_GF16[mi * SNOVA_l2 + i1] + SNOVA_q - Fvv_in_GF16Matrix[mi * SNOVA_l2 + i1] % SNOVA_q) %
					    SNOVA_q;
				}
			}

		// Whipped F21
		alignas(32) uint16_t whipped_F21[SNOVA_m * SNOVA_l * SNOVA_o * SNOVA_l2] = {0};

		for (int mi = 0; mi < SNOVA_m; mi++)
			for (int idx = 0; idx < SNOVA_o; idx++)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int nj = 0; nj < SNOVA_v; ++nj)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								for (int k1 = 0; k1 < SNOVA_l; k1++)
									whipped_F21[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_l2 + i1 * SNOVA_l + j1] +=
									    skx->F21[((mi * SNOVA_o + idx) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + k1] *
									    whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_l2 + k1 * SNOVA_l + j1];

		for (int i1 = 0; i1 < SNOVA_m * SNOVA_o * SNOVA_l * SNOVA_l2; i1++) {
			whipped_F21[i1] = whipped_F21[i1] % SNOVA_q;
		}

#ifndef SYMMETRIC
		// Whipped F12
		alignas(32) uint16_t whipped_F12[SNOVA_m * SNOVA_l * SNOVA_o * SNOVA_l2] = {0};

		for (int mi = 0; mi < SNOVA_m; mi++)
			for (int nj = 0; nj < SNOVA_v; ++nj)
				for (int b1 = 0; b1 < SNOVA_l; ++b1)
					for (int idx = 0; idx < SNOVA_o; idx++)
						for (int i1 = 0; i1 < SNOVA_l; i1++)
							for (int j1 = 0; j1 < SNOVA_l; j1++)
								for (int k1 = 0; k1 < SNOVA_l; k1++)
									whipped_F12[((mi * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_l2 + i1 * SNOVA_l + j1] +=
									    skx->F12[((mi * SNOVA_v + nj) * SNOVA_o + idx) * SNOVA_l2 + k1 * SNOVA_l + i1] *
									    whipped_sig[(b1 * SNOVA_v + nj) * SNOVA_l2 + k1 * SNOVA_l + j1];

		for (int i1 = 0; i1 < SNOVA_m * SNOVA_o * SNOVA_l * SNOVA_l2; i1++) {
			whipped_F12[i1] = whipped_F12[i1] % SNOVA_q;
		}

#else
#define whipped_F12 whipped_F21
#endif

		// compute the coefficients of Xo and put into gauss matrix and compute
		// the coefficients of Xo^t and add into gauss matrix
		for (int mi = 0; mi < SNOVA_m; mi++)
			for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
				alignas(32) uint16_t gfm_temp0[SNOVA_o * SNOVA_l2] = {0};
				alignas(32) uint16_t gfm_temp1[SNOVA_o * SNOVA_l2] = {0};
				alignas(32) uint16_t gfm_temp2[SNOVA_o * SNOVA_l2] = {0};

				int mi_prime = i_prime(mi, alpha);

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int b1 = 0; b1 < SNOVA_l; ++b1) {
						__m256i *gfm_temp0_256 = (__m256i *)&gfm_temp0[idx * SNOVA_l2];
						__m256i q2_256 = _mm256_set1_epi16(skx->q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
						__m256i *wf21 = (__m256i *)&whipped_F21[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_l2];

						*gfm_temp0_256 += _mm256_mullo_epi16(q2_256, *wf21);
					}

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
					gfm_temp0[i1] = gfm_temp0[i1] % SNOVA_q;
				}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++) {
								gfm_temp1[idx * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    gfm_temp0[idx * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    skx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + k1 * SNOVA_l + j1];
							}

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
					gfm_temp1[i1] = gfm_temp1[i1] % SNOVA_q;
				}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++) {
								gfm_temp2[idx * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    skx->Q1[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    gfm_temp1[idx * SNOVA_l2 + k1 * SNOVA_l + j1];
							}

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
					gfm_temp2[i1] = gfm_temp2[i1] % SNOVA_q;
				}

				for (int idx = 0; idx < SNOVA_o; idx++) {
					__m256i gfm[4];
					uint16_t *gfm_16 = (uint16_t *)gfm;
					__m256i am[4];
					uint16_t *am_16 = (uint16_t *)am;

					for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
						for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
							for (int tj2 = 0; tj2 < SNOVA_l; tj2++) {
								gfm_16[ti2 * SNOVA_l2 + tj1 * SNOVA_l + tj2] = gfm_temp2[idx * SNOVA_l2 + tj1 * SNOVA_l + ti2];
							}

					for (int ti1 = 0; ti1 < SNOVA_l; ti1++)
						for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
							for (int tj2 = 0; tj2 < SNOVA_l; tj2++) {
								am_16[ti1 * SNOVA_l2 + tj1 * SNOVA_l + tj2] =
								    skx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + ti1 * SNOVA_l + tj2];
							}

					for (int ti1 = 0; ti1 < SNOVA_l; ti1++)
						for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
							__m256i *gauss_256 = (__m256i *)&gauss16[mi * SNOVA_l2 + ti1 * SNOVA_l + ti2][idx * SNOVA_l2];

							*gauss_256 += _mm256_mullo_epi16(gfm[ti2], am[ti1]);
						}
				}
			}

		for (int mi = 0; mi < SNOVA_m; mi++) {
			for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
				alignas(32) uint16_t gfm_temp0[SNOVA_o * SNOVA_l2] = {0};
				alignas(32) uint16_t gfm_temp1[SNOVA_o * SNOVA_l2] = {0};
				alignas(32) uint16_t gfm_temp2[SNOVA_o * SNOVA_l2] = {0};

				int mi_prime = i_prime(mi, alpha);

				// Transpose
				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int b1 = 0; b1 < SNOVA_l; ++b1) {
						__m256i *gfm_temp0_256 = (__m256i *)&gfm_temp0[idx * SNOVA_l2];
						__m256i q1_256 = _mm256_set1_epi16(skx->q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1]);
						__m256i *wf12 = (__m256i *)&whipped_F12[((mi_prime * SNOVA_l + b1) * SNOVA_o + idx) * SNOVA_l2];

						*gfm_temp0_256 += _mm256_mullo_epi16(q1_256, *wf12);
					}

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
					gfm_temp0[i1] = gfm_temp0[i1] % SNOVA_q;
				}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								gfm_temp1[idx * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    skx->Am[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    gfm_temp0[idx * SNOVA_l2 + j1 * SNOVA_l + k1];

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
					gfm_temp1[i1] = gfm_temp1[i1] % SNOVA_q;
				}

				for (int idx = 0; idx < SNOVA_o; idx++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								gfm_temp2[idx * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    gfm_temp1[idx * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    skx->Q2[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + k1 * SNOVA_l + j1];

				for (int i1 = 0; i1 < SNOVA_o * SNOVA_l2; i1++) {
					gfm_temp2[i1] = gfm_temp2[i1] % SNOVA_q;
				}

				for (int idx = 0; idx < SNOVA_o; idx++) {
					__m256i gfm[4];
					uint16_t *gfm_16 = (uint16_t *)gfm;
					__m256i bm[4];
					uint16_t *bm_16 = (uint16_t *)bm;

					for (int ti1 = 0; ti1 < SNOVA_l; ti1++)
						for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
							for (int tj2 = 0; tj2 < SNOVA_l; tj2++) {
								gfm_16[ti1 * SNOVA_l2 + tj1 * SNOVA_l + tj2] = gfm_temp2[idx * SNOVA_l2 + ti1 * SNOVA_l + tj1];
							}

					for (int ti2 = 0; ti2 < SNOVA_l; ti2++)
						for (int tj1 = 0; tj1 < SNOVA_l; tj1++)
							for (int tj2 = 0; tj2 < SNOVA_l; tj2++) {
								bm_16[ti2 * SNOVA_l2 + tj1 * SNOVA_l + tj2] =
								    skx->Bm[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + tj2 * SNOVA_l + ti2];
							}

					for (int ti1 = 0; ti1 < SNOVA_l; ti1++)
						for (int ti2 = 0; ti2 < SNOVA_l; ti2++) {
							__m256i *gauss_256 = (__m256i *)&gauss16[mi * SNOVA_l2 + ti1 * SNOVA_l + ti2][idx * SNOVA_l2];

							*gauss_256 += _mm256_mullo_epi16(gfm[ti1], bm[ti2]);
						}
				}
			}
		}

		for (int ti = 0; ti < SNOVA_m * SNOVA_l2; ti++)
			for (int tj = 0; tj < SNOVA_m * SNOVA_l2; tj++) {
				gauss[ti][tj] = (gauss[ti][tj] + gauss16[ti][tj]) % SNOVA_q;
			}

		// Gaussian elimination in constant time
		for (int i = 0; i < SNOVA_m * SNOVA_l2; ++i) {
			for (int j = i + 1; j < SNOVA_m * SNOVA_l2; ++j) {
				int16_t mask = ct_is_not_zero(gauss[i][i]) - 1;
				for (int k = 0; k < SNOVA_ml2; ++k) {
					gauss[i][k] += mask & gauss[j][k];
				}
			}

			for (int k = 0; k < SNOVA_ml2; ++k) {
				gauss[i][k] = gauss[i][k] % SNOVA_q;
			}

			flag_redo |= 1 - ct_is_not_zero(gauss[i][i]);

			uint16_t t_GF16 = ct_gf_inverse(gauss[i][i]);
			for (int k = 0; k < SNOVA_ml2; ++k) {
				gauss[i][k] = gauss[i][k] * t_GF16;
			}

			for (int k = 0; k < SNOVA_ml2; ++k) {
				gauss[i][k] = gauss[i][k] % SNOVA_q;
			}

			for (int j = i + 1; j < SNOVA_m * SNOVA_l2; ++j) {
				uint16_t gji = SNOVA_q - gauss[j][i];
				for (int k = 0; k < SNOVA_ml2; ++k) {
					gauss[j][k] += gauss[i][k] * gji;
				}
			}

			// A periodic full cleanup is needed to prevent uint16_t overflow
			if (!(i % 64)) {
				for (int j = i + 1; j < SNOVA_m * SNOVA_l2; ++j)
					for (int k = 0; k < SNOVA_ml2; ++k) {
						gauss[j][k] = gauss[j][k] % SNOVA_q;
					}
			} else {
				for (int j = i + 1; j < SNOVA_m * SNOVA_l2; ++j) {
					gauss[j][i + 1] = gauss[j][i + 1] % SNOVA_q;
				}
			}
		}

		if (!flag_redo) {
			// Last step of Gaussian elimination
			memset(solution, 0, sizeof(solution));

			for (int i = SNOVA_m * SNOVA_l2 - 1; i >= 0; --i) {
				uint16_t sum = 0;
				for (int k = i + 1; k < SNOVA_m * SNOVA_l2; ++k) {
					sum += gauss[i][k] * solution[k];
				}
				solution[i] = (gauss[i][SNOVA_m * SNOVA_l2] + SNOVA_q - (sum % SNOVA_q)) % SNOVA_q;
			}

			for (int idx = 0; idx < SNOVA_o * SNOVA_l2; ++idx) {
				signature_in_GF[idx + SNOVA_v * SNOVA_l2] = solution[idx];
			}

			uint16_t signature_in_GF16[SNOVA_n * SNOVA_l2] = {0};

			// Establish signature using T12
			for (int index = 0; index < SNOVA_v; ++index)
				for (int mi = 0; mi < SNOVA_m; ++mi)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							for (int k1 = 0; k1 < SNOVA_l; k1++)
								signature_in_GF16[index * SNOVA_l2 + i1 * SNOVA_l + j1] +=
								    skx->T12[(index * SNOVA_m + mi) * SNOVA_l2 + i1 * SNOVA_l + k1] *
								    solution[mi * SNOVA_l2 + k1 * SNOVA_l + j1];

			for (int index = 0; index < SNOVA_v; ++index)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						signature_in_GF[index * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    (signature_in_GF16[index * SNOVA_l2 + i1 * SNOVA_l + j1] % SNOVA_q +
						     signature_in_GF[index * SNOVA_l2 + i1 * SNOVA_l + j1] + SNOVA_q) %
						    SNOVA_q;

#ifdef SYMMETRIC
			// Reject if the signature has symmetric matrices
			int num_sym = 0;
			for (int idx = 0; idx < SNOVA_n; ++idx) {
				int is_symmetric = 1;
				for (int i1 = 0; i1 < SNOVA_l - 1; i1++)
					for (int j1 = i1 + 1; j1 < SNOVA_l; j1++) {
						is_symmetric &= signature_in_GF[idx * SNOVA_l2 + i1 * SNOVA_l + j1] ==
						                signature_in_GF[idx * SNOVA_l2 + j1 * SNOVA_l + i1];
					}
				num_sym += is_symmetric;
			}
			flag_redo = num_sym > 0;
#endif
		}
	} while (flag_redo);

	compress_gf(sig, signature_in_GF, SNOVA_n * SNOVA_l2);
	memcpy(sig + BYTES_SIGNATURE - BYTES_SALT, salt, BYTES_SALT);

	return 0;
}

/**
 * PK expansion.
 */
int SNOVA_NAMESPACE(pk_expand)(expanded_PK *pkx, const uint8_t *pk) {
	SNOVA_INIT

	memcpy(pkx->pk_seed, pk, SEED_LENGTH_PUBLIC);

	/**
	 * Create P matrix
	 */

	uint8_t pk_bytes[NUM_GEN_PUB_BYTES];

	snova_pk_expander_t instance;
	snova_pk_expander_init(&instance, pk, SEED_LENGTH_PUBLIC);
	snova_pk_expander(pk_bytes, NUM_GEN_PUB_BYTES, &instance);

#ifdef SYMMETRIC
	gf_t A[SNOVA_m * SNOVA_alpha * 2 * (SNOVA_l2 + SNOVA_l)];
	gf_t pk_gf[NUM_GEN_PUB_GF];
	convert_bytes_to_GF(pk_gf, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);

	gf_t *curval = pk_gf;

	// Copy P11 and P12
	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int ni = 0; ni < SNOVA_v; ++ni) {
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = i1; j1 < SNOVA_l; j1++) {
					pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
					curval++;
				}

			for (int nj = ni + 1; nj < SNOVA_v; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}

			for (int nj = SNOVA_v; nj < SNOVA_n; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}
		}

	for (int idx = 0; idx < SNOVA_m * SNOVA_alpha * 2 * (SNOVA_l2 + SNOVA_l); idx++) {
		A[idx] = *curval;
		curval++;
	}

	// Copy P22
	gf_t P22[NUMGF_PK] = {0};
	if (expand_gf(P22, pk + SEED_LENGTH_PUBLIC, NUMGF_PK)) {
		return -1;
	}

	curval = P22;

	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int ni = SNOVA_v; ni < SNOVA_n; ++ni)
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				for (int j1 = i1; j1 < SNOVA_l; j1++) {
					pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + ni) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
					curval++;
				}

				for (int nj = ni + 1; nj < SNOVA_n; ++nj)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] = *curval;
						curval++;
					}
			}

	// Make symmetric
	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int ni = 0; ni < SNOVA_n; ++ni) {
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = i1 + 1; j1 < SNOVA_l; j1++)
					pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] =
					    pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + ni) * SNOVA_l2 + i1 * SNOVA_l + j1];

			for (int nj = ni + 1; nj < SNOVA_n; ++nj)
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						pkx->P[((mi * SNOVA_n + nj) * SNOVA_n + ni) * SNOVA_l2 + j1 * SNOVA_l + i1] =
						    pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + j1];
		}

#else

	gf_t P22[SNOVA_m * SNOVA_o * SNOVA_o * SNOVA_l2];
	if (expand_gf(P22, pk + SEED_LENGTH_PUBLIC, NUMGF_PK)) {
		return -1;
	}

	gf_t P_matrix[NUM_PUB_GF];
	gf_t *P11 = P_matrix;
	gf_t *P12 = P_matrix + SNOVA_m * SNOVA_v * SNOVA_v * SNOVA_l2;
	gf_t *P21 = P_matrix + SNOVA_m * SNOVA_v * SNOVA_n * SNOVA_l2;

	convert_bytes_to_GF(P_matrix, (uint8_t *)pk_bytes, NUM_GEN_PUB_GF);

	for (int mi = 0; mi < SNOVA_m; ++mi) {
		for (int ni = 0; ni < SNOVA_v; ++ni) {
			for (int nj = 0; nj < SNOVA_v; ++nj) {
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    P11[((mi * SNOVA_v + ni) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1];
			}

			for (int nj = SNOVA_v; nj < SNOVA_n; ++nj) {
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    P12[((mi * SNOVA_v + ni) * SNOVA_o + (nj - SNOVA_v)) * SNOVA_l2 + i1 * SNOVA_l + j1];
			}
		}
		for (int ni = SNOVA_v; ni < SNOVA_n; ++ni) {
			for (int nj = 0; nj < SNOVA_v; ++nj) {
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    P21[((mi * SNOVA_o + (ni - SNOVA_v)) * SNOVA_v + nj) * SNOVA_l2 + i1 * SNOVA_l + j1];
			}

			for (int nj = SNOVA_v; nj < SNOVA_n; ++nj) {
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    P22[((mi * SNOVA_o + (ni - SNOVA_v)) * SNOVA_o + nj - SNOVA_v) * SNOVA_l2 + i1 * SNOVA_l + j1];
			}
		}
	}

	gf_t *A = P_matrix + (SNOVA_o * (SNOVA_n * SNOVA_n - SNOVA_o * SNOVA_o)) * SNOVA_l2;
#endif

	/**
	 * Create AB matrices, improve q
	 */
	gf_t *B = A + SNOVA_o * SNOVA_alpha * SNOVA_l2;
	gf_t *q1 = A + 2 * SNOVA_o * SNOVA_alpha * SNOVA_l2;
	gf_t *q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
		be_invertible_by_add_aS(&(pkx->Am[idx * SNOVA_l2]), &A[idx * SNOVA_l2]);
		be_invertible_by_add_aS(&(pkx->Bm[idx * SNOVA_l2]), &B[idx * SNOVA_l2]);

		if (!q1[idx * SNOVA_l + SNOVA_l - 1]) {
			q1[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q1[idx * SNOVA_l] + (q1[idx * SNOVA_l] == 0));
		}
		if (!q2[idx * SNOVA_l + SNOVA_l - 1]) {
			q2[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q2[idx * SNOVA_l] + (q2[idx * SNOVA_l] == 0));
		}
	}

	for (int idx = 0; idx < SNOVA_m * SNOVA_alpha * SNOVA_l; ++idx) {
		pkx->q1[idx] = q1[idx];
		pkx->q2[idx] = q2[idx];
	}

	return 0;
}

/**
 * Optimized version of verify.
 */
int SNOVA_NAMESPACE(verify)(const expanded_PK *pkx, const uint8_t *sig, const uint8_t *digest, size_t len_digest) {
	SNOVA_INIT

	gf_t signature_in_GF[NUMGF_SIGNATURE];
	if (expand_gf(signature_in_GF, sig, NUMGF_SIGNATURE)) {
		return -1;
	}

	uint8_t Am[SNOVA_m * SNOVA_alpha * SNOVA_l2];
	uint8_t Bm[SNOVA_m * SNOVA_alpha * SNOVA_l2];
	const uint8_t *q1 = pkx->q1;
	const uint8_t *q2 = pkx->q2;

	for (int idx = 0; idx < SNOVA_m * SNOVA_alpha * SNOVA_l2; ++idx) {
		Am[idx] = pkx->Am[idx];
		Bm[idx] = pkx->Bm[idx];
	}

#ifdef SYMMETRIC
	// Reject if the signature has symmetric matrices
	int num_sym = 0;
	for (int idx = 0; idx < SNOVA_n; ++idx) {
		int is_symmetric = 1;
		for (int i1 = 0; i1 < SNOVA_l - 1; i1++)
			for (int j1 = i1 + 1; j1 < SNOVA_l; j1++) {
				is_symmetric &=
				    signature_in_GF[idx * SNOVA_l2 + i1 * SNOVA_l + j1] == signature_in_GF[idx * SNOVA_l2 + j1 * SNOVA_l + i1];
			}
		num_sym += is_symmetric;
	}
	if (num_sym > 0) {
		return -1;
	}
#endif

	/**
	 * Whip signature
	 */
	uint16_t whipped_sig[SNOVA_l * SNOVA_n * SNOVA_l2] = {0};

	for (int ab = 0; ab < SNOVA_l; ++ab)
		for (int idx = 0; idx < SNOVA_n; ++idx)
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = 0; j1 < SNOVA_l; j1++)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						whipped_sig[idx * SNOVA_l * SNOVA_l2 + i1 * SNOVA_l2 + ab * SNOVA_l + j1] +=
						    gf_S[ab * SNOVA_l2 + i1 * SNOVA_l + k1] * signature_in_GF[idx * SNOVA_l2 + k1 * SNOVA_l + j1];

	for (int i1 = 0; i1 < SNOVA_l * SNOVA_n * SNOVA_l2; i1++) {
		whipped_sig[i1] = whipped_sig[i1] % SNOVA_q;
	}

	/**
	 * Evaluate whipped central map
	 */
	uint16_t hash_in_GF[SNOVA_m * SNOVA_l2] = {0};
	uint16_t sum_t1[SNOVA_m * SNOVA_l2 * SNOVA_l2] = {0};
	uint16_t sum_t1s[SNOVA_m * SNOVA_l2 * SNOVA_l2] = {0};

	for (int mi = 0; mi < SNOVA_m; ++mi) {
		for (int ni = 0; ni < SNOVA_n; ++ni) {
			uint16_t sum_t0[SNOVA_l * SNOVA_l2] = {0};

			// Right
			for (int nj = 0; nj < SNOVA_n; ++nj)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int b1 = 0; b1 < SNOVA_l2; ++b1)
							sum_t0[i1 * SNOVA_l2 + b1] +=
							    pkx->P[((mi * SNOVA_n + ni) * SNOVA_n + nj) * SNOVA_l2 + i1 * SNOVA_l + k1] *
							    whipped_sig[nj * SNOVA_l * SNOVA_l2 + k1 * SNOVA_l2 + b1];

			for (int b1 = 0; b1 < SNOVA_l2 * SNOVA_l; ++b1) {
				sum_t0[b1] = sum_t0[b1] % SNOVA_q;
			}

			// Left, transposed whipped_sig
			for (int a1 = 0; a1 < SNOVA_l; ++a1)
				for (int k1 = 0; k1 < SNOVA_l; k1++)
					for (int i1 = 0; i1 < SNOVA_l; i1++)
						for (int b1 = 0; b1 < SNOVA_l2; ++b1)
							sum_t1[(mi * SNOVA_l2 + a1 * SNOVA_l) * SNOVA_l2 + i1 * SNOVA_l2 + b1] +=
							    whipped_sig[ni * SNOVA_l * SNOVA_l2 + k1 * SNOVA_l2 + a1 * SNOVA_l + i1] *
							    sum_t0[k1 * SNOVA_l2 + b1];
		}
	}

	for (int i1 = 0; i1 < SNOVA_m * SNOVA_l2 * SNOVA_l2; i1++) {
		sum_t1[i1] = sum_t1[i1] % SNOVA_q;
	}

	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int a1 = 0; a1 < SNOVA_l; ++a1)
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int b1 = 0; b1 < SNOVA_l2; ++b1)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						sum_t1s[(mi * SNOVA_l + a1) * SNOVA_l * SNOVA_l2 + b1 * SNOVA_l2 + i1 * SNOVA_l + j1] =
						    sum_t1[(mi * SNOVA_l + a1) * SNOVA_l * SNOVA_l2 + i1 * SNOVA_l2 + b1 * SNOVA_l + j1];

	/**
	 * Apply A, B, q1 and q2, aka E matrix
	 */
	for (int mi = 0; mi < SNOVA_m; ++mi)
		for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
			int mi_prime = i_prime(mi, alpha);

			uint16_t gfm_temp1[SNOVA_l2] = {0};
			uint16_t gfm_temp2[SNOVA_l2] = {0};

			// apply q1 and q2
			for (int a1 = 0; a1 < SNOVA_l; ++a1) {
				uint16_t gfm_temp0[SNOVA_l2] = {0};

				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int b1 = 0; b1 < SNOVA_l; ++b1)
						for (int j1 = 0; j1 < SNOVA_l; j1++)
							gfm_temp0[i1 * SNOVA_l + j1] +=
							    q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b1] *
							    sum_t1s[(mi_prime * SNOVA_l + a1) * SNOVA_l * SNOVA_l2 + b1 * SNOVA_l2 + i1 * SNOVA_l + j1];

				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					gfm_temp0[i1] = gfm_temp0[i1] % SNOVA_q;
				}

				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++) {
						gfm_temp1[i1 * SNOVA_l + j1] +=
						    gfm_temp0[i1 * SNOVA_l + j1] * q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a1];
					}
			}

			for (int i1 = 0; i1 < SNOVA_l2; i1++) {
				gfm_temp1[i1] = gfm_temp1[i1] % SNOVA_q;
			}

			// A and B
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = 0; j1 < SNOVA_l; j1++)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						gfm_temp2[i1 * SNOVA_l + j1] +=
						    gfm_temp1[i1 * SNOVA_l + k1] * Bm[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + k1 * SNOVA_l + j1];

			for (int i1 = 0; i1 < SNOVA_l2; i1++) {
				gfm_temp2[i1] = gfm_temp2[i1] % SNOVA_q;
			}

			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = 0; j1 < SNOVA_l; j1++)
					for (int k1 = 0; k1 < SNOVA_l; k1++)
						hash_in_GF[mi * SNOVA_l2 + i1 * SNOVA_l + j1] +=
						    Am[(mi * SNOVA_alpha + alpha) * SNOVA_l2 + i1 * SNOVA_l + k1] * gfm_temp2[k1 * SNOVA_l + j1];
		}

	for (int i1 = 0; i1 < SNOVA_m * SNOVA_l2; i1++) {
		hash_in_GF[i1] = hash_in_GF[i1] % SNOVA_q;
	}

	/**
	 * Check hashes
	 */
	uint8_t signed_bytes[BYTES_HASH];
	uint8_t signed_gf[GF16_HASH] = {0};
	const uint8_t *salt = sig + BYTES_SIGNATURE - BYTES_SALT;
	hash_combined(signed_bytes, digest, len_digest, pkx->pk_seed, salt);
	expand_gf(signed_gf, signed_bytes, GF16_HASH);

	for (int i = 0; i < GF16_HASH; ++i) {
		if (hash_in_GF[i] != signed_gf[i]) {
			return -1;
		}
	}

	return 0;
}
