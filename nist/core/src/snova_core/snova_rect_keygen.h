
#ifdef RCT_PROFILE
#include <x86intrin.h>
#include <stdio.h>
#define RCT_PT(x) do{ x = __rdtsc(); }while(0)
static uint64_t rct_pf[8];
#define RCT_PACC(i,a,b) rct_pf[i]+=(b)-(a)
static uint64_t rct_pf6[3];
#define RCT_PACC6(i,a,b) rct_pf6[i]+=(b)-(a)
static uint64_t rct_pf6b[3];
#define RCT_PACC6B(i,a,b) rct_pf6b[i]+=(b)-(a)
#else
#define RCT_PT(x) do{ (void)(x); }while(0)
#define RCT_PACC(i,a,b)
#define RCT_PACC6(i,a,b)
#define RCT_PACC6B(i,a,b)
#endif

#ifdef RCT_VPROF
#include <x86intrin.h>
#include <stdio.h>
#define VPT(x) do{ x = __rdtsc(); }while(0)
#define VPA(i,a,b) vpf[i]+=(b)-(a)
#ifndef VPROF_AT
#define VPROF_AT 2000
#endif
#else
#define VPT(x)
#define VPA(i,a,b)
#endif

static int rct_genkeys(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
    rct_init();
#if RCT_SIGN_JOG
    rct_sj_ensure();
#endif
    uint64_t _pa,_pb;

#if SNOVA_WRAPPER_STACK
    _Alignas(32) gf_t T12[SNOVA_o * SNOVA_v * SNOVA_l2 + 16];
    _Alignas(32) gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
#else
    static gf_t T12[SNOVA_o * SNOVA_v * SNOVA_l2 + 16];
    static gf_t P22[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
#endif
    memset(T12 + SNOVA_o * SNOVA_v * SNOVA_l2, 0, 16);
    memset(P22, 0, sizeof(P22));

    expand_T12(T12, seed + SEED_LENGTH_PUBLIC);

#if RCT_KG_STREAM
    (void)_pa; (void)_pb;
#if SNOVA_WRAPPER_STACK
    _Alignas(32) gf_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2 + 16];
#else
    static gf_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2 + 16];
#endif
    memset(F12, 0, sizeof(F12));
    {
        snova_pgen_t pg;
        snova_pgen_init(&pg, seed);
        snova_prow_t rw;
#if RCT_Q_SIMD
        _Static_assert((uint32_t)SNOVA_v * (SNOVA_q - 1) * (SNOVA_q - 1)
                           + (SNOVA_q - 1) < 65536u,
                       "u16 accumulation overflow guard");
        _Static_assert(2u * (uint32_t)SNOVA_v * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                       "u16 accumulation overflow guard");
        _Alignas(32) uint16_t rowacc[SNOVA_o * SNOVA_l2];
        while (snova_pgen_next_row(&pg, &rw)) {
            const int i1 = rw.mi, j1 = rw.ni;
            if (rw.block == 0) {
                memset(rowacc, 0, sizeof(rowacc));
                for (int idx = 0; idx < SNOVA_v; idx++)
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        rct_q_matmul4_add(&rowacc[k1 * SNOVA_l2],
                                          &rw.cells[(size_t)idx * SNOVA_l2],
                                          &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
                gf_t *frow = &F12[((size_t)(i1 * SNOVA_v + j1) * SNOVA_o) * SNOVA_l2];
                for (int t = 0; t < SNOVA_o * SNOVA_l2; t++)
                    frow[t] = (gf_t)(rowacc[t] % SNOVA_q);
            } else if (rw.block == 1) {
                gf_t *frow = &F12[((size_t)(i1 * SNOVA_v + j1) * SNOVA_o) * SNOVA_l2];
                for (int t = 0; t < SNOVA_o * SNOVA_l2; t++)
                    frow[t] = (gf_t)((frow[t] + rw.cells[t]) % SNOVA_q);
            } else {
                memset(rowacc, 0, sizeof(rowacc));
                for (int idx = 0; idx < SNOVA_v; idx++)
                    for (int k1 = 0; k1 < SNOVA_o; k1++) {
                        rct_q_matmul4_add(&rowacc[k1 * SNOVA_l2],
                                          &T12[(idx * SNOVA_o + j1) * SNOVA_l2],
                                          &F12[((size_t)(i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]);
                        rct_q_matmul4_add(&rowacc[k1 * SNOVA_l2],
                                          &rw.cells[(size_t)idx * SNOVA_l2],
                                          &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
                    }
                gf_t *prow = &P22[((size_t)(i1 * SNOVA_o + j1) * SNOVA_o) * SNOVA_l2];
                for (int t = 0; t < SNOVA_o * SNOVA_l2; t++)
                    prow[t] = (gf_t)((SNOVA_q - (rowacc[t] % SNOVA_q)) % SNOVA_q);
            }
        }
        SNOVA_CLEAR_OBJ(rowacc);
#elif RCT_KGM_ACTIVE
        while (snova_pgen_next_row(&pg, &rw)) {
            const int i1 = rw.mi, j1 = rw.ni;
            if (rw.block == 0) {
                __m256i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm256_setzero_si256();
                for (int idx = 0; idx < SNOVA_v; idx++) {
                    __m256i bsh[4];
                    rct_gf16_bshuf_exp(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&rw.cells[(size_t)idx * SNOVA_l2])), bsh);
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        acc[k1] = _mm256_xor_si256(acc[k1], rct_gf16_mm4_bs(bsh,
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&T12[(idx * SNOVA_o + k1) * SNOVA_l2]))));
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&F12[((size_t)(i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
                                     gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc[k1])));
            } else if (rw.block == 1) {
                gf_t *frow = &F12[((size_t)(i1 * SNOVA_v + j1) * SNOVA_o) * SNOVA_l2];
                for (int t = 0; t < SNOVA_o * SNOVA_l2; t++) gf_set_add(&frow[t], rw.cells[t]);
            } else {
                __m256i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm256_setzero_si256();
                for (int idx = 0; idx < SNOVA_v; idx++) {
                    __m256i bshT[4], bshP[4];
                    rct_gf16_bshuf_exp(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&T12[(idx * SNOVA_o + j1) * SNOVA_l2])), bshT);
                    rct_gf16_bshuf_exp(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&rw.cells[(size_t)idx * SNOVA_l2])), bshP);
                    for (int k1 = 0; k1 < SNOVA_o; k1++) {
                        acc[k1] = _mm256_xor_si256(acc[k1], rct_gf16_mm4_bs(bshT,
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&F12[((size_t)(i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]))));
                        acc[k1] = _mm256_xor_si256(acc[k1], rct_gf16_mm4_bs(bshP,
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&T12[(idx * SNOVA_o + k1) * SNOVA_l2]))));
                    }
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&P22[((size_t)(i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2],
                                     gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc[k1])));
            }
        }
#else
        while (snova_pgen_next_row(&pg, &rw)) {
            const int i1 = rw.mi, j1 = rw.ni;
            if (rw.block == 0) {
                gf_t *frow = &F12[((size_t)(i1 * SNOVA_v + j1) * SNOVA_o) * SNOVA_l2];
#if RCT_F5_A4 || RCT_F5_M4
                rct_f5_fold_row_bsec(frow, rw.cells, T12);
#else
                for (int idx = 0; idx < SNOVA_v; idx++)
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        RCT_MATMUL_ADD_BSEC(&frow[k1 * SNOVA_l2],
                                            &rw.cells[(size_t)idx * SNOVA_l2],
                                            &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
#endif
            } else if (rw.block == 1) {
                gf_t *frow = &F12[((size_t)(i1 * SNOVA_v + j1) * SNOVA_o) * SNOVA_l2];
                for (int t = 0; t < SNOVA_o * SNOVA_l2; t++) gf_set_add(&frow[t], rw.cells[t]);
            } else {
                gf_t *prow = &P22[((size_t)(i1 * SNOVA_o + j1) * SNOVA_o) * SNOVA_l2];
#if RCT_F5_A4 || RCT_F5_M4
                rct_f5_fold_row_bsec(prow, rw.cells, T12);
#else
                for (int idx = 0; idx < SNOVA_v; idx++)
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        RCT_MATMUL_ADD_BSEC(&prow[k1 * SNOVA_l2],
                                            &rw.cells[(size_t)idx * SNOVA_l2],
                                            &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
#endif
            }
        }
#if RCT_F5_A4 || RCT_F5_M4
        rct_f5_fold_P22_pass1(P22, T12, F12);
#else
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int idx = 0; idx < SNOVA_v; idx++)
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    for (int j1 = 0; j1 < SNOVA_o; j1++)
                        RCT_MATMUL_ADD_ASEC(&P22[((size_t)(i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2],
                                            &T12[(idx * SNOVA_o + j1) * SNOVA_l2],
                                            &F12[((size_t)(i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]);
#endif
#if SNOVA_q != 16
        for (int t = 0; t < SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2; t++)
            P22[t] = (SNOVA_q - P22[t]) % SNOVA_q;
#endif
#endif
    }
#else
#if SNOVA_WRAPPER_STACK
    _Alignas(32) gf_t P_matrix[NUM_PUB_GF];
#else
    gf_t *P_matrix = rct_pub_Pmatrix;
#endif
    RCT_PT(_pa);
    expand_public(P_matrix, seed);
    RCT_PT(_pb); RCT_PACC(0,_pa,_pb);

    gf_t *P11 = P_matrix;
    gf_t *P12 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_v * SNOVA_l2;
    gf_t *P21 = P_matrix + SNOVA_m1 * SNOVA_v * SNOVA_n * SNOVA_l2;
#if SNOVA_WRAPPER_STACK
    _Alignas(32) gf_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2 + 16];
#else
    static gf_t F12[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2 + 16];
#endif
    memset(F12, 0, sizeof(F12));

    RCT_PT(_pa);
#if RCT_Q_SIMD
    _Static_assert((uint32_t)SNOVA_v * (SNOVA_q - 1) * (SNOVA_q - 1)
                       + (SNOVA_q - 1) < 65536u,
                   "u16 accumulation overflow guard");
    _Static_assert(2u * (uint32_t)SNOVA_v * (SNOVA_q - 1) * (SNOVA_q - 1) < 65536u,
                   "u16 accumulation overflow guard");
    {
        static uint16_t F12u[SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2];
        memset(F12u, 0, sizeof(F12u));
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++)
                for (int idx = 0; idx < SNOVA_v; idx++)
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        rct_q_matmul4_add(&F12u[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
                                          &P11[((i1 * SNOVA_v + j1) * SNOVA_v + idx) * SNOVA_l2],
                                          &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++)
            F12[i1] = (gf_t)((F12u[i1] + P12[i1]) % SNOVA_q);
        SNOVA_CLEAR_OBJ(F12u);
    }
    RCT_PT(_pb); RCT_PACC(1,_pa,_pb);

    RCT_PT(_pa);
    {
        static uint16_t P22u[SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2];
        memset(P22u, 0, sizeof(P22u));
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int idx = 0; idx < SNOVA_v; idx++)
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    for (int j1 = 0; j1 < SNOVA_o; j1++)
                        rct_q_matmul4_add(&P22u[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2],
                                          &T12[(idx * SNOVA_o + j1) * SNOVA_l2],
                                          &F12[((i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]);
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_o; j1++)
                for (int idx = 0; idx < SNOVA_v; idx++)
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        rct_q_matmul4_add(&P22u[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2],
                                          &P21[((i1 * SNOVA_o + j1) * SNOVA_v + idx) * SNOVA_l2],
                                          &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
        for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2; i1++)
            P22[i1] = (gf_t)((SNOVA_q - (P22u[i1] % SNOVA_q)) % SNOVA_q);
        SNOVA_CLEAR_OBJ(P22u);
    }
    RCT_PT(_pb); RCT_PACC(2,_pa,_pb);
#elif RCT_KGM_ACTIVE
    {
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_v; j1++) {
                __m256i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm256_setzero_si256();
                for (int idx = 0; idx < SNOVA_v; idx++) {
                    __m256i bsh[4];
                    rct_gf16_bshuf_exp(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&P11[((i1 * SNOVA_v + j1) * SNOVA_v + idx) * SNOVA_l2])), bsh);
                    for (int k1 = 0; k1 < SNOVA_o; k1++)
                        acc[k1] = _mm256_xor_si256(acc[k1], rct_gf16_mm4_bs(bsh,
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&T12[(idx * SNOVA_o + k1) * SNOVA_l2]))));
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
                                     gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc[k1])));
            }
    }
    RCT_PT(_pb); RCT_PACC(1,_pa,_pb);

    for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++)
        gf_set_add(&F12[i1], P12[i1]);

    RCT_PT(_pa);
    {
        for (int i1 = 0; i1 < SNOVA_m1; i1++)
            for (int j1 = 0; j1 < SNOVA_o; j1++) {
                __m256i acc[SNOVA_o];
                for (int k1 = 0; k1 < SNOVA_o; k1++) acc[k1] = _mm256_setzero_si256();
                for (int idx = 0; idx < SNOVA_v; idx++) {
                    __m256i bshT[4], bshP[4];
                    rct_gf16_bshuf_exp(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&T12[(idx * SNOVA_o + j1) * SNOVA_l2])), bshT);
                    rct_gf16_bshuf_exp(_mm256_cvtepu8_epi16(_mm_loadu_si128(
                        (const __m128i *)&P21[((i1 * SNOVA_o + j1) * SNOVA_v + idx) * SNOVA_l2])), bshP);
                    for (int k1 = 0; k1 < SNOVA_o; k1++) {
                        acc[k1] = _mm256_xor_si256(acc[k1], rct_gf16_mm4_bs(bshT,
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&F12[((i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]))));
                        acc[k1] = _mm256_xor_si256(acc[k1], rct_gf16_mm4_bs(bshP,
                            _mm256_cvtepu8_epi16(_mm_loadu_si128(
                                (const __m128i *)&T12[(idx * SNOVA_o + k1) * SNOVA_l2]))));
                    }
                }
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    _mm_storeu_si128((__m128i *)&P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2],
                                     gf16_pack_u16_to_bytes(gf16_compress_u16x16(acc[k1])));
            }
    }
    RCT_PT(_pb); RCT_PACC(2,_pa,_pb);
#else
#if RCT_F5_A4 || RCT_F5_M4 || RCT_F5_WIDE
    rct_f5_fold_F12(F12, P11, T12);
#elif RCT_L4G_FOLD
    rct_l4g_fold_bsec(F12, P11, T12, SNOVA_v, 0);
#else
    for (int i1 = 0; i1 < SNOVA_m1; i1++)
        for (int j1 = 0; j1 < SNOVA_v; j1++)
            for (int idx = 0; idx < SNOVA_v; idx++)
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    RCT_MATMUL_ADD_BSEC(&F12[((i1 * SNOVA_v + j1) * SNOVA_o + k1) * SNOVA_l2],
                                   &P11[((i1 * SNOVA_v + j1) * SNOVA_v + idx) * SNOVA_l2],
                                   &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
#endif
    RCT_PT(_pb); RCT_PACC(1,_pa,_pb);

    for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_v * SNOVA_o * SNOVA_l2; i1++)
        gf_set_add(&F12[i1], P12[i1]);

    RCT_PT(_pa);
#if RCT_F5_A4 || RCT_F5_M4
    rct_f5_fold_P22_pass1(P22, T12, F12);
    rct_f5_fold_P22_pass2(P22, P21, T12);
#elif RCT_L4G_FOLD
    rct_l4g_fold_asec(P22, T12, F12);
    rct_l4g_fold_bsec(P22, P21, T12, SNOVA_o, 1);
#else
    for (int i1 = 0; i1 < SNOVA_m1; i1++)
        for (int idx = 0; idx < SNOVA_v; idx++)
            for (int k1 = 0; k1 < SNOVA_o; k1++)
                for (int j1 = 0; j1 < SNOVA_o; j1++)
                    RCT_MATMUL_ADD_ASEC(&P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2],
                                   &T12[(idx * SNOVA_o + j1) * SNOVA_l2],
                                   &F12[((i1 * SNOVA_v + idx) * SNOVA_o + k1) * SNOVA_l2]);

    for (int i1 = 0; i1 < SNOVA_m1; i1++)
        for (int j1 = 0; j1 < SNOVA_o; j1++)
            for (int idx = 0; idx < SNOVA_v; idx++)
                for (int k1 = 0; k1 < SNOVA_o; k1++)
                    RCT_MATMUL_ADD_BSEC(&P22[((i1 * SNOVA_o + j1) * SNOVA_o + k1) * SNOVA_l2],
                                   &P21[((i1 * SNOVA_o + j1) * SNOVA_v + idx) * SNOVA_l2],
                                   &T12[(idx * SNOVA_o + k1) * SNOVA_l2]);
#endif

#if SNOVA_q != 16
    for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_o * SNOVA_o * SNOVA_l2; i1++)
        P22[i1] = (SNOVA_q - P22[i1]) % SNOVA_q;
#endif
    RCT_PT(_pb); RCT_PACC(2,_pa,_pb);
#endif
#endif
#ifdef RCT_PROFILE
    fprintf(stderr,"[KG] expand_pub=%lu F12=%lu P22=%lu\n",(unsigned long)rct_pf[0],(unsigned long)rct_pf[1],(unsigned long)rct_pf[2]);
    rct_pf[0]=rct_pf[1]=rct_pf[2]=0;
#endif

    memcpy(pk, seed, SEED_LENGTH_PUBLIC);
    compress_pk(pk + SEED_LENGTH_PUBLIC, P22);
    memcpy(sk, seed, SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE);
#if HASH_PK
    shake256(sk + SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE, BYTES_PK_HASH, pk, BYTES_PK);
#endif
    SNOVA_CLEAR_OBJ(T12);
    SNOVA_CLEAR_OBJ(F12);
    return 0;
}
