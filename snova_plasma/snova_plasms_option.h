#ifndef PLASMA_OPTION_H
#define PLASMA_OPTION_H

// Optimisation level. 0: Ref, 1: Opt, 2: VTL
#ifndef OPTIMISATION
#if __AVX2__
#define OPTIMISATION 2
#else
#define OPTIMISATION 1
#endif
#endif

#if OPTIMISATION == 2 && rank == 2
#include "plasma_2x2/snova_plasma_2x2_avx2.h"
#include "plasma_2x2/snova_plasma_2x2_avx2_verify.h"
#include "plasma_2x2/snova_plasma_2x2_avx2_sign.h"
#elif OPTIMISATION == 2 && rank == 4
#include "plasma_4x4/snova_plasma_4x4_avx2.h"
#include "plasma_4x4/snova_plasma_4x4_avx2_sign.h"
#include "plasma_4x4/snova_plasma_4x4_avx2_verify.h"
#elif OPTIMISATION == 2
#include "plasma_3x3/snova_plasma_3x3_avx2.h"
#include "plasma_3x3/snova_plasma_3x3_avx2_sign.h"
#include "plasma_3x3/snova_plasma_3x3_avx2_verify.h"
#elif OPTIMISATION == 1
#include "snova_opt.h"
#endif

// ---
#if OPTIMISATION == 2
#define gen_P22 gen_P22_vtl

#if rank == 4
#define gen_F gen_F_4x4_vtl
#define sign_digest_core sign_digest_core_4x4_avx2_vtl
#define verify_core verify_signture_4x4

#elif rank == 2
#define gen_F gen_F_2x2_vtl
#define sign_digest_core sign_digest_core_2x2_avx2_vtl
#define verify_core verify_signture_2x2

#else  // 3
#define gen_F gen_F_vtl
#define sign_digest_core sign_digest_core_gnl_vtl
#define verify_core verify_signture_vtl

#endif

// ---
#elif OPTIMISATION == 1
#define gen_F gen_F_opt
#define gen_P22 gen_P22_ref
#define sign_digest_core sign_digest_core_opt
#define verify_core verify_signture_opt
// ---
#else
#define gen_F gen_F_ref
#define gen_P22 gen_P22_ref
#define sign_digest_core sign_digest_core_ref
#define verify_core verify_signture_ref

#endif

#endif