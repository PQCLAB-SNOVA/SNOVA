/**
 * @file primitives/gf16_simd.h
 */
#ifndef SNOVA_PRIMITIVES_GF16_SIMD_H
#define SNOVA_PRIMITIVES_GF16_SIMD_H

#include "platform.h"

#if defined(SNOVA_ARCH_X86_AVX2)
  #include "../platforms/x86_avx2/gf16_simd_avx2.h"
  #include "../platforms/x86_avx2/xgf16_avx.h"
#elif defined(SNOVA_ARCH_WASM_SIMD128)
  #include "../platforms/wasm_simd128/gf16_simd_wasm.h"
#elif defined(SNOVA_ARCH_PORTABLE_OPT)
  #include "../platforms/portable_opt/gf16_simd_opt.h"
#elif defined(SNOVA_ARCH_REF)
  #include "../platforms/ref/gf16_simd_ref.h"
#endif

#endif
