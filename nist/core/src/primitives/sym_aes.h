/**
 * @file primitives/sym_aes.h
 */
#ifndef SNOVA_PRIMITIVES_SYM_AES_H
#define SNOVA_PRIMITIVES_SYM_AES_H

#include "platform.h"

#if defined(SNOVA_ARCH_X86_AVX2)
  #include "../platforms/x86_avx2/aes_ni.h"
#elif defined(SNOVA_ARCH_WASM_SIMD128)
  #include "../platforms/wasm_simd128/aes_wasm.h"
#elif defined(SNOVA_ARCH_PORTABLE_OPT)
  #include "../platforms/portable_opt/aes_opt.h"
#elif defined(SNOVA_ARCH_REF)
  #include "../platforms/ref/aes_ref.h"
#endif

#endif
