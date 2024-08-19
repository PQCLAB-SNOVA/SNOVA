#ifndef WASM_API_H
#define WASM_API_H

#define NO_MT4B
#include "snova.h"

#define CRYPTO_PUBLICKEYBYTES bytes_pk
#define CRYPTO_BYTES (bytes_signature + bytes_salt)

#ifndef CRYPTO_ALGNAME
#define CRYPTO_ALGNAME "SNOVA_24_5_4"
#endif

#ifndef EM_PORT_API
#    if defined(__EMSCRIPTEN__)
#        include <emscripten.h>
#        if defined(__cplusplus)
#            define EM_PORT_API(rettype) extern "C" rettype EMSCRIPTEN_KEEPALIVE
#        else
#            define EM_PORT_API(rettype) rettype EMSCRIPTEN_KEEPALIVE
#        endif
#    else
#        if defined(__cplusplus)
#            define EM_PORT_API(rettype) extern "C" rettype
#        else
#            define EM_PORT_API(rettype) rettype
#        endif
#    endif
#endif

EM_PORT_API(int) getSeedLength(void);

EM_PORT_API(int) getSskLength(void);

EM_PORT_API(int) getEskLength(void);

EM_PORT_API(int) getPkLength(void);

EM_PORT_API(int) getSaltLength(void);

EM_PORT_API(int) getSignLength(void);

EM_PORT_API(void) safeFree(void* ptr, size_t size);

EM_PORT_API(void) genKeyPairSsk(uint8_t* pk, uint8_t* ssk, const uint8_t* seed);

EM_PORT_API(void) genKeyPairEsk(uint8_t* pk, uint8_t* esk, const uint8_t* seed);

EM_PORT_API(void) genPkWithSsk(uint8_t* pk, uint8_t* ssk);

EM_PORT_API(void) genPkWithEsk(uint8_t* pk, uint8_t* esk);

EM_PORT_API(void) signWithSsk(uint8_t *sm, const uint8_t *m, const size_t mlen, uint8_t *salt, const uint8_t *ssk);

EM_PORT_API(void) signWithEsk(uint8_t *sm, const uint8_t *m, const size_t mlen, uint8_t *salt,const uint8_t *esk);

EM_PORT_API(int) verify(const uint8_t* sm, const uint8_t *m, const size_t mlen, const uint8_t* pk);

#endif
