/**
 * @file primitives/platform.h
 */
#ifndef SNOVA_PRIMITIVES_PLATFORM_H
#define SNOVA_PRIMITIVES_PLATFORM_H

#if defined(SNOVA_ARCH_X86_AVX2)
#elif defined(SNOVA_ARCH_WASM_SIMD128)
#elif defined(SNOVA_ARCH_PORTABLE_OPT)
#elif defined(SNOVA_ARCH_REF)
#else
  #error "snova: No SNOVA_ARCH_* defined. Use -DSNOVA_ARCH_X86_AVX2=1, -DSNOVA_ARCH_WASM_SIMD128=1, -DSNOVA_ARCH_PORTABLE_OPT=1, or -DSNOVA_ARCH_REF=1."
#endif

#endif
