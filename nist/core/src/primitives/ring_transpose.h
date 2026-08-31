/**
 * @file primitives/ring_transpose.h
 */
#ifndef SNOVA_PRIMITIVES_RING_TRANSPOSE_H
#define SNOVA_PRIMITIVES_RING_TRANSPOSE_H

#include "platform.h"

#if defined(SNOVA_ARCH_X86_AVX2)
  #include "../platforms/x86_avx2/ring_transpose_avx2.h"
#elif defined(SNOVA_ARCH_WASM_SIMD128)
  #include "../platforms/wasm_simd128/ring_transpose_wasm.h"
#elif defined(SNOVA_ARCH_PORTABLE_OPT)
  #include "../platforms/portable_opt/ring_transpose_opt.h"
#elif defined(SNOVA_ARCH_REF)
  #include "../platforms/ref/ring_transpose_ref.h"
#endif

#endif
