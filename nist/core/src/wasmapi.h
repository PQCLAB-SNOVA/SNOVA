/**
 * @file wasmapi.h
 */
#ifndef WASM_API_H
#define WASM_API_H

#include <stddef.h>
#include <stdint.h>
#include "snova.h"

#ifndef EM_PORT_API
#  if defined(__EMSCRIPTEN__)
#    include <emscripten.h>
#    if defined(__cplusplus)
#      define EM_PORT_API(rettype) extern "C" rettype EMSCRIPTEN_KEEPALIVE
#    else
#      define EM_PORT_API(rettype) rettype EMSCRIPTEN_KEEPALIVE
#    endif
#  else
#    if defined(__cplusplus)
#      define EM_PORT_API(rettype) extern "C" rettype
#    else
#      define EM_PORT_API(rettype) rettype
#    endif
#  endif
#endif

EM_PORT_API(int) getSeedLength(void);
EM_PORT_API(int) getSkLength(void);
EM_PORT_API(int) getPkLength(void);
EM_PORT_API(int) getSkxLength(void);
EM_PORT_API(int) getPkxLength(void);
EM_PORT_API(int) getSaltLength(void);
EM_PORT_API(int) getSignLength(void);

EM_PORT_API(void) safeFree(void *ptr, size_t size);

EM_PORT_API(int) keygen(uint8_t *pk, uint8_t *sk, const uint8_t *seed);

EM_PORT_API(void) expandSkx(uint8_t *skx, const uint8_t *sk);
EM_PORT_API(int)  signWithSkx(uint8_t *sm, const uint8_t *m, const size_t mlen,
                              const uint8_t *salt, const uint8_t *skx);
EM_PORT_API(void) clearSkx(uint8_t *skx);

EM_PORT_API(int) expandPk(uint8_t *pkx, const uint8_t *pk);
EM_PORT_API(int) verifyPkx(const uint8_t *sm, const uint8_t *m, const size_t mlen,
                           const uint8_t *pkx);

#endif
