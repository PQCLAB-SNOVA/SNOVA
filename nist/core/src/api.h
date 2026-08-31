/**
 * @file api.h
 */
#ifndef SNOVA_API_H
#define SNOVA_API_H

#include "snova.h"

#ifndef sk_is_seed
  #define sk_is_seed SNOVA_SK_IS_SEED
#endif

#define CRYPTO_PUBLICKEYBYTES SNOVA_PK_BYTES
#define CRYPTO_SECRETKEYBYTES (SNOVA_SEED_LEN + SNOVA_BYTES_PK_HASH)
#define CRYPTO_BYTES          SNOVA_SM_OVERHEAD

#ifndef CRYPTO_ALGNAME
  #define CRYPTO_ALGNAME "SNOVA"
#endif

#ifdef __cplusplus
extern "C" {
#endif

static inline void snova_init(void) {   }

static inline int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
    return snova_crypto_sign_keypair(pk, sk);
}

static inline int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                              const unsigned char *m, unsigned long long mlen,
                              const unsigned char *sk) {
    size_t s = 0;
    int rc = snova_crypto_sign(sm, &s, m, (size_t)mlen, sk);
    *smlen = (unsigned long long)s;
    return rc;
}

static inline int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                                   const unsigned char *sm, unsigned long long smlen,
                                   const unsigned char *pk) {
    size_t s = 0;
    int rc = snova_crypto_sign_open(m, &s, sm, (size_t)smlen, pk);
    *mlen = (unsigned long long)s;
    return rc;
}

#ifdef __cplusplus
}
#endif

#endif
