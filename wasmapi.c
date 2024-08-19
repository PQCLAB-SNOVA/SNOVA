#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "wasmapi.h"

EM_PORT_API(int) getSeedLength(void) {
    return seed_length;
}

EM_PORT_API(int) getSskLength(void) {
    return seed_length;
}

EM_PORT_API(int) getEskLength(void) {
    return bytes_sk;
}

EM_PORT_API(int) getPkLength(void) {
    return bytes_pk;
}

EM_PORT_API(int) getSaltLength(void) {
    return bytes_salt;
}

EM_PORT_API(int) getSignLength(void) {
    return bytes_signature + bytes_salt;
}

EM_PORT_API(void) safeFree(void* ptr, size_t size) {
    if (ptr == NULL) {
        return;
    }
    memset(ptr, 0, size);
    free(ptr);
}

EM_PORT_API(void) genKeyPairSsk(uint8_t* pk, uint8_t* ssk, const uint8_t* seed) {
    const uint8_t* pkseed = seed;
    const uint8_t* skseed = seed + seed_length_public;

    generate_keys_ssk(pk, ssk, pkseed, skseed);
}

EM_PORT_API(void) genKeyPairEsk(uint8_t* pk, uint8_t* esk, const uint8_t* seed) {
    const uint8_t* pkseed = seed;
    const uint8_t* skseed = seed + seed_length_public;

    generate_keys_esk(pk, esk, pkseed, skseed);
}

EM_PORT_API(void) genPkWithSsk(uint8_t* pk, uint8_t* ssk) {
    generate_pk_with_ssk(pk, ssk);
}

EM_PORT_API(void) genPkWithEsk(uint8_t* pk, uint8_t* esk) {
    generate_pk_with_esk(pk, esk);
}

EM_PORT_API(void) signWithSsk(uint8_t *sm, const uint8_t *m, const size_t mlen, uint8_t *salt, const uint8_t *ssk) {
    uint8_t digest[64];

    shake256(m, mlen, digest, 64);
    sign_digest_ssk(sm, digest, 64, salt, ssk);
}

EM_PORT_API(void) signWithEsk(uint8_t *sm, const uint8_t *m, const size_t mlen, uint8_t *salt,const uint8_t *esk) {
    uint8_t digest[64];

    shake256(m, mlen, digest, 64);
    sign_digest_esk(sm, digest, 64, salt, esk);
}

EM_PORT_API(int) verify(const uint8_t* sm, const uint8_t *m, const size_t mlen, const uint8_t* pk) {
    uint8_t digest[64];
    
    shake256(m, mlen, digest, 64);
    return verify_signture(digest, 64, sm, pk);
}
