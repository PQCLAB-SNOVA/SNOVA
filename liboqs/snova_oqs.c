// SPDX-License-Identifier: MIT

/**
 * Glue code between SNOVA and liboqs
 *
 * Copyright (c) 2025 SNOVA TEAM
 */

#include <oqs/oqs.h>

#include "snova.h"
#include "symmetric.h"

// Size of the message digest in the hash-and-sign paragdigm.
#define BYTES_DIGEST 64

OQS_STATUS SNOVA_NAMESPACE(oqs_keypair)(uint8_t *pk, uint8_t *sk) {
	uint8_t seed_pair[SEED_LENGTH];
	uint8_t *pt_private_key_seed;
	uint8_t *pt_public_key_seed;

	OQS_randombytes(seed_pair, SEED_LENGTH);
	pt_public_key_seed = seed_pair;

	int res = SNOVA_NAMESPACE(genkeys)(pk, sk, pt_public_key_seed);

	if (res) {
		return OQS_ERROR;
	} else {
		return OQS_SUCCESS;
	}
}

OQS_STATUS SNOVA_NAMESPACE(oqs_sign)(uint8_t *signature, size_t *signature_len, const uint8_t *message, size_t message_len,
                                     const uint8_t *secret_key) {
	expanded_SK skx_d;
	uint8_t salt[BYTES_SALT];

	OQS_randombytes(salt, BYTES_SALT);

	int res = SNOVA_NAMESPACE(sk_expand)(&skx_d, secret_key);
	if (res) {
		return OQS_ERROR;
	}

#if SNOVA_q == 16
	// No change to KATs
	uint8_t digest[BYTES_DIGEST];
	shake256(digest, BYTES_DIGEST, message, message_len);
	res = SNOVA_NAMESPACE(sign)(&skx_d, signature, digest, BYTES_DIGEST, salt);
#else
	res = SNOVA_NAMESPACE(sign)(&skx_d, signature, message, message_len, salt);
#endif

	if (res) {
		return OQS_ERROR;
	} else {
		*signature_len = BYTES_SIGNATURE;
		return OQS_SUCCESS;
	}
}

OQS_STATUS SNOVA_NAMESPACE(oqs_verify)(const uint8_t *signature, size_t signature_len, const uint8_t *message,
                                       size_t message_len, const uint8_t *pk) {
	expanded_PK pkx;

	if (signature_len != BYTES_SIGNATURE) {
		return OQS_ERROR;
	}

	int res = SNOVA_NAMESPACE(pk_expand)(&pkx, pk);
	if (res) {
		return OQS_ERROR;
	}

#if SNOVA_q == 16
	// No change to KATs
	uint8_t digest[BYTES_DIGEST];
	shake256(digest, BYTES_DIGEST, message, message_len);
	res = SNOVA_NAMESPACE(verify)(&pkx, signature, digest, BYTES_DIGEST);
#else
	res = SNOVA_NAMESPACE(verify)(&pkx, signature, message, message_len);
#endif

	if (res) {
		return OQS_ERROR;
	} else {
		return OQS_SUCCESS;
	}
}
