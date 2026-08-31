/**
 * @file wasmapi.c
 */
#include "wasmapi.h"
#include "primitives/sym_shake.h"
#include <stdlib.h>
#include <string.h>

EM_PORT_API(int) getSeedLength(void) { return SNOVA_SEED_LEN; }
EM_PORT_API(int) getSkLength(void)   { return SNOVA_SK_BYTES; }
EM_PORT_API(int) getPkLength(void)   { return SNOVA_PK_BYTES; }
EM_PORT_API(int) getSkxLength(void)  { return (int)snova_skx_bytes(); }
EM_PORT_API(int) getPkxLength(void)  { return (int)snova_pkx_bytes(); }
EM_PORT_API(int) getSaltLength(void) { return SNOVA_SALT_BYTES; }
EM_PORT_API(int) getSignLength(void) { return SNOVA_SM_OVERHEAD;   }

EM_PORT_API(void) safeFree(void *ptr, size_t size) {
    if (ptr == NULL) return;
    memset(ptr, 0, size);
    free(ptr);
}

static void wasmapi_shake256_64(const uint8_t *m, size_t mlen, uint8_t digest[64]) {
    keccak_ctx c;
    shake_init(&c, 256);
    shake_absorb(&c, m, mlen);
    shake_finalize(&c);
    shake_squeeze(&c, digest, 64);
}

EM_PORT_API(int) keygen(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
    return snova_keygen(pk, sk, seed);
}

EM_PORT_API(void) expandSkx(uint8_t *skx, const uint8_t *sk) {
    snova_sk_expand_skx(skx, sk);
}

EM_PORT_API(int) signWithSkx(uint8_t *sm, const uint8_t *m, const size_t mlen,
                             const uint8_t *salt, const uint8_t *skx) {
    uint8_t digest[64];
    wasmapi_shake256_64(m, mlen, digest);
    return snova_sign_digest_skx(sm, digest, 64, salt, skx);
}

EM_PORT_API(void) clearSkx(uint8_t *skx) {
    snova_skx_clear(skx);
}

EM_PORT_API(int) expandPk(uint8_t *pkx, const uint8_t *pk) {
    return snova_pk_expand(pkx, pk);
}

EM_PORT_API(int) verifyPkx(const uint8_t *sm, const uint8_t *m, const size_t mlen,
                           const uint8_t *pkx) {
    uint8_t digest[64];
    wasmapi_shake256_64(m, mlen, digest);
    return snova_verify_digest_pkx(digest, 64, sm, pkx);
}
