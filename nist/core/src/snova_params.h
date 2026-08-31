/**
 * @file snova_params.h
 */
#ifndef SNOVA_PARAMS_H
#define SNOVA_PARAMS_H

#ifndef SNOVA_V
#define SNOVA_V 24
#endif
#ifndef SNOVA_O
#define SNOVA_O 5
#endif
#ifndef SNOVA_L
#define SNOVA_L 4
#endif

#ifndef SNOVA_Q
#define SNOVA_Q 16
#endif

#if SNOVA_Q == 11
#define SNOVA_Q_A 0
#define SNOVA_Q_B 3
#define SNOVA_Q_C 6
#define SNOVA_PACK_GF 16
#define SNOVA_PACK_BYTES 7
#elif SNOVA_Q == 13
#define SNOVA_Q_A 2
#define SNOVA_Q_B 11
#define SNOVA_Q_C 3
#define SNOVA_PACK_GF 15
#define SNOVA_PACK_BYTES 7
#elif SNOVA_Q == 16
#define SNOVA_PACK_GF 2
#define SNOVA_PACK_BYTES 1
#elif SNOVA_Q == 19
#define SNOVA_Q_A 1
#define SNOVA_Q_B 3
#define SNOVA_Q_C 15
#define SNOVA_PACK_GF 15
#define SNOVA_PACK_BYTES 8
#else
#error "Unsupported SNOVA_Q (supported: 11, 13, 16, 19)"
#endif
#define SNOVA_BYTES_GF(x) ((SNOVA_PACK_BYTES * (x) + SNOVA_PACK_GF - 1) / SNOVA_PACK_GF)
#define SNOVA_REJECTION_LIMIT ((256 / SNOVA_Q) * SNOVA_Q)
#ifndef SNOVA_R
#define SNOVA_R SNOVA_L
#endif
#ifndef SNOVA_M1
#define SNOVA_M1 ((SNOVA_O * SNOVA_R) / SNOVA_L)
#endif
#define SNOVA_M2 (SNOVA_O * SNOVA_L * SNOVA_R)
#ifdef SNOVA_M2_ASSERT
_Static_assert(SNOVA_M2_ASSERT == (SNOVA_O * SNOVA_L * SNOVA_R),
               "SNOVA_M2 supplied by the build disagrees with the derived o*l*r "
               "(build configuration vs snova_params.h formula divergence)");
#endif

#ifndef FIXED_ABQ
#define FIXED_ABQ 1
#endif
#ifndef HASH_PK
#define HASH_PK (SNOVA_l > 2)
#endif
#ifndef SNOVA_BYTES_PK_HASH
#define SNOVA_BYTES_PK_HASH 48
#endif
#ifndef ROUND2_T12
#define ROUND2_T12 0
#endif
#ifndef ABQ_ALG2
#define ABQ_ALG2 1
#endif

#define SNOVA_SK_IS_SEED 1
#ifndef SNOVA_PK_EXPAND_SHAKE
#define SNOVA_PK_EXPAND_SHAKE 0
#endif

#define SNOVA_N      (SNOVA_V + SNOVA_O)
#define SNOVA_M      (SNOVA_O)
#define SNOVA_L2     (SNOVA_L * SNOVA_L)
#ifndef SNOVA_ALPHA
#define SNOVA_ALPHA  (SNOVA_L * SNOVA_R + 2 * SNOVA_R)
#endif
#define SNOVA_RANK   (SNOVA_L)
#define SNOVA_SQ_RANK (SNOVA_RANK * SNOVA_RANK)

#define SNOVA_SEED_PUB   16
#define SNOVA_SEED_PRIV  32
#define SNOVA_SEED_LEN   (SNOVA_SEED_PUB + SNOVA_SEED_PRIV)
#define SNOVA_SALT_BYTES 16

#define SNOVA_GF16S_HASH   (SNOVA_O * SNOVA_L * SNOVA_R)
#define SNOVA_GF16S_SIG    (SNOVA_N * SNOVA_L * SNOVA_R)
#if SNOVA_Q != 16
#define SNOVA_BYTES_HASH   (SNOVA_BYTES_GF(SNOVA_GF16S_HASH))
#define SNOVA_BYTES_SIG    (SNOVA_BYTES_GF(SNOVA_GF16S_SIG))
#else
#define SNOVA_BYTES_HASH   ((SNOVA_GF16S_HASH + 1) >> 1)
#define SNOVA_BYTES_SIG    ((SNOVA_GF16S_SIG + 1) >> 1)
#endif
#define SNOVA_BYTES_SIG_SALT (SNOVA_BYTES_SIG + SNOVA_SALT_BYTES)

#define SNOVA__STR2(x) #x
#define SNOVA__STR(x) SNOVA__STR2(x)
#define SNOVA_NAME ("SNOVA_" SNOVA__STR(SNOVA_V) "_" SNOVA__STR(SNOVA_O) "_" SNOVA__STR(SNOVA_L))

#if SNOVA_R == SNOVA_L
_Static_assert((SNOVA_O * SNOVA_R) % SNOVA_L == 0,
    "seven-param: m1 = (o*r)/l must divide evenly for square cells");
_Static_assert(SNOVA_M1 == SNOVA_O,
    "seven-param: m1 must collapse to o for square cells");
_Static_assert(SNOVA_M2 == SNOVA_M * SNOVA_L2,
    "seven-param: m2 must equal o*l^2 (GAUSS_ROW invariant)");
#endif

#endif
