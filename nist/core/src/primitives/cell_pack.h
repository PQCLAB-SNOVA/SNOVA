/**
 * @file primitives/cell_pack.h
 */
#ifndef SNOVA_PRIMITIVES_CELL_PACK_H
#define SNOVA_PRIMITIVES_CELL_PACK_H

#include "platform.h"

#if defined(SNOVA_ARCH_X86_AVX2)
  #include "../platforms/x86_avx2/cell_pack_avx2.h"
#elif defined(SNOVA_ARCH_WASM_SIMD128)
  #include "../platforms/wasm_simd128/cell_pack_wasm.h"
#elif defined(SNOVA_ARCH_PORTABLE_OPT)
  #include "../platforms/portable_opt/cell_pack_opt.h"
#elif defined(SNOVA_ARCH_REF)
  #include "../platforms/ref/cell_pack_ref.h"
#endif

#endif
