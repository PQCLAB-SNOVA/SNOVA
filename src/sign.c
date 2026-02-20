// SPDX-License-Identifier: MIT

/**
 * Interface to NIST API.
 *
 * This file is the only point where randombytes is used by SNOVA.
 * The snova_* implementations are deterministic.
 *
 * SNOVA Team 2025
 */

#include <string.h>

#if defined(VALGRIND)
#include <valgrind/memcheck.h>
#endif

#include "api.h"
#include "rng.h"
#include "symmetric.h"

#define BYTES_DIGEST 64

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk) {
	uint8_t seed[SEED_LENGTH];

	randombytes(seed, SEED_LENGTH);

#if defined(VALGRIND)
	VALGRIND_MAKE_MEM_UNDEFINED(seed, SEED_LENGTH);
#endif

	int res = SNOVA_NAMESPACE(genkeys)(pk, sk, seed);

	return res;
}

int crypto_sign(unsigned char* sm, unsigned long long* smlen, const unsigned char* m, unsigned long long mlen,
                const unsigned char *sk) {
	uint8_t digest[BYTES_DIGEST];
	uint8_t salt[BYTES_SALT];
	expanded_SK skx_d;

	int res = SNOVA_NAMESPACE(sk_expand)(&skx_d, sk);
	if (res) {
		return res;
	}

	randombytes(salt, BYTES_SALT);

#if defined(VALGRIND)
	VALGRIND_MAKE_MEM_UNDEFINED(salt, BYTES_SALT);
#endif

	shake256(digest, BYTES_DIGEST, m, mlen);

	res = SNOVA_NAMESPACE(sign)(&skx_d, sm, digest, BYTES_DIGEST, salt);
	if (!res) {
		memcpy(sm + CRYPTO_BYTES, m, mlen);
		*smlen = mlen + CRYPTO_BYTES;
	}

	return res;
}

int crypto_sign_open(unsigned char* m, unsigned long long* mlen, const unsigned char* sm, unsigned long long smlen,
                     const unsigned char *pk) {
	uint8_t digest[BYTES_DIGEST];
	expanded_PK pkx;

	if (smlen < CRYPTO_BYTES) {
		return -1;
	}

	int res = SNOVA_NAMESPACE(pk_expand)(&pkx, pk);
	if (res) {
		return -1;
	}

	shake256(digest, BYTES_DIGEST, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);

	res = SNOVA_NAMESPACE(verify)(&pkx, sm, digest, BYTES_DIGEST);
	if (!res) {
		memcpy(m, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);
		*mlen = smlen - CRYPTO_BYTES;
	} else {
		return -1;
	}

	return 0;
}
