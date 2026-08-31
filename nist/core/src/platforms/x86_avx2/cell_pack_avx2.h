/**
 * @file platforms/x86_avx2/cell_pack_avx2.h
 */
#ifndef SNOVA_PLATFORMS_X86_AVX2_CELL_PACK_AVX2_H
#define SNOVA_PLATFORMS_X86_AVX2_CELL_PACK_AVX2_H

#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include "../../snova_params.h"
#include "../../gf16_core/gf16m.h"

#ifndef CELL_PACK_V
  #define CELL_PACK_V       SNOVA_V
  #define CELL_PACK_O       SNOVA_O
  #define CELL_PACK_L       SNOVA_L
  #define CELL_PACK_LSQ     SNOVA_L2
  #define CELL_PACK_VTL_V_PAD ((((SNOVA_V) * (SNOVA_L) + 31) / 32) * 32)
  #define CELL_PACK_VTL_O_PAD ((((SNOVA_O) * (SNOVA_L) + 31) / 32) * 32)
  #define CELL_PACK_RANK_PAIR_HALF (((SNOVA_L) + 1) / 2)
  #define CELL_PACK_L_PAIR_HALF    ((SNOVA_L) / 2)
  #define CELL_PACK_L_FLOOR2       (CELL_PACK_L_PAIR_HALF * 2)
  #define CELL_PACK_L_HAS_ODD      ((SNOVA_L) & 1)
#endif

/**
 * @param AJ        Output byte-flat buffer (size: mainCol*L × mainRow_pad)
 * @param A         Input gf16_t cell array
 * @param mainCol   Number of cell rows in source
 * @param mainRow   Number of cell cols in source
 * @param mainRow_pad  Stride of AJ rows in bytes (must be ≥ mainRow * L, typically VTL_*_PAD)
 */
static inline void snova_cell_pack_jogress_LL(
    uint8_t *AJ, const gf16_t *A, int mainCol, int mainRow, int mainRow_pad)
{
#if CELL_PACK_L == 3
    for (int mi = 0; mi < mainCol; ++mi)
        for (int mj = 0; mj < mainRow; ++mj) {
            const uint8_t *cell = (const uint8_t *)(A + (mi * mainRow + mj) * 9);
            uint32_t r0 = ((uint32_t)cell[0]      ) | ((uint32_t)cell[1] <<  8)
                        | ((uint32_t)cell[2] << 16);
            uint32_t r1 = ((uint32_t)cell[3]      ) | ((uint32_t)cell[4] <<  8)
                        | ((uint32_t)cell[5] << 16);
            uint32_t r2 = ((uint32_t)cell[6]      ) | ((uint32_t)cell[7] <<  8)
                        | ((uint32_t)cell[8] << 16);
            r0 &= 0x000F0F0F; r1 &= 0x000F0F0F; r2 &= 0x000F0F0F;
            memcpy(&AJ[(mi * 3 + 0) * mainRow_pad + mj * 3], &r0, sizeof r0);
            memcpy(&AJ[(mi * 3 + 1) * mainRow_pad + mj * 3], &r1, sizeof r1);
            memcpy(&AJ[(mi * 3 + 2) * mainRow_pad + mj * 3], &r2, sizeof r2);
        }
#elif CELL_PACK_L == 5
    const __m128i mask_lo128 = _mm_set1_epi8(0x0F);
    for (int mi = 0; mi < mainCol; ++mi)
        for (int mj = 0; mj < mainRow; ++mj) {
            const uint8_t *cell = (const uint8_t *)(A + (mi * mainRow + mj) * 25);
            __m128i m_lo = _mm_loadu_si128((const __m128i *)cell);
            __m128i m_hi = _mm_loadu_si128((const __m128i *)(cell + 9));
            m_lo = _mm_and_si128(m_lo, mask_lo128);
            m_hi = _mm_and_si128(m_hi, mask_lo128);
            #define PACK_ROW5_LO(EI) do {                                         \
                __m128i mask = _mm_setr_epi8(                                     \
                    (EI) * 5,     (EI) * 5 + 1, (EI) * 5 + 2,                    \
                    (EI) * 5 + 3, (EI) * 5 + 4,                                  \
                    (char)0x80, (char)0x80, (char)0x80, (char)0x80,              \
                    (char)0x80, (char)0x80, (char)0x80, (char)0x80,              \
                    (char)0x80, (char)0x80, (char)0x80);                         \
                __m128i row = _mm_shuffle_epi8(m_lo, mask);                       \
                uint64_t v = (uint64_t)_mm_cvtsi128_si64(row);                    \
                memcpy(&AJ[(mi * 5 + (EI)) * mainRow_pad + mj * 5], &v, sizeof v);\
            } while (0)
            #define PACK_ROW5_HI(EI, BASE) do {                                   \
                __m128i mask = _mm_setr_epi8(                                     \
                    (BASE),     (BASE) + 1, (BASE) + 2,                          \
                    (BASE) + 3, (BASE) + 4,                                      \
                    (char)0x80, (char)0x80, (char)0x80, (char)0x80,              \
                    (char)0x80, (char)0x80, (char)0x80, (char)0x80,              \
                    (char)0x80, (char)0x80, (char)0x80);                         \
                __m128i row = _mm_shuffle_epi8(m_hi, mask);                       \
                uint64_t v = (uint64_t)_mm_cvtsi128_si64(row);                    \
                memcpy(&AJ[(mi * 5 + (EI)) * mainRow_pad + mj * 5], &v, sizeof v);\
            } while (0)
            PACK_ROW5_LO(0);
            PACK_ROW5_LO(1);
            PACK_ROW5_LO(2);
            PACK_ROW5_HI(3, 6);
            PACK_ROW5_HI(4, 11);
            #undef PACK_ROW5_LO
            #undef PACK_ROW5_HI
        }
#else
    for (int mi = 0; mi < mainCol; ++mi)
        for (int mj = 0; mj < mainRow; ++mj)
            for (int ei = 0; ei < CELL_PACK_L; ++ei)
                for (int ej = 0; ej < CELL_PACK_L; ++ej)
                    AJ[(mi * CELL_PACK_L + ei) * mainRow_pad + (mj * CELL_PACK_L + ej)] =
                        A[(mi * mainRow + mj) * CELL_PACK_L * CELL_PACK_L + ei * CELL_PACK_L + ej] & 0x0F;
#endif
}

static inline void snova_cell_pack_jogress_Tr_F21vo(
    uint8_t F21_vo_tr_J[CELL_PACK_V * CELL_PACK_L][CELL_PACK_VTL_O_PAD],
    const gf16m_t F21_mp[CELL_PACK_O][CELL_PACK_V])
{
#if CELL_PACK_L == 3
    for (int vj = 0; vj < CELL_PACK_V; ++vj)
        for (int ej = 0; ej < CELL_PACK_L; ++ej) {
            uint8_t *out = &F21_vo_tr_J[vj * CELL_PACK_L + ej][0];
            for (int oi = 0; oi < CELL_PACK_O; ++oi) {
                const uint8_t *cell = (const uint8_t *)F21_mp[oi][vj];
                uint32_t r = ((uint32_t)cell[0 * 3 + ej]      ) |
                             ((uint32_t)cell[1 * 3 + ej] <<  8) |
                             ((uint32_t)cell[2 * 3 + ej] << 16);
                r &= 0x000F0F0F;
                memcpy(&out[oi * CELL_PACK_L], &r, sizeof r);
            }
        }
#else
    for (int vj = 0; vj < CELL_PACK_V; ++vj)
        for (int oi = 0; oi < CELL_PACK_O; ++oi)
            for (int ej = 0; ej < CELL_PACK_L; ++ej)
                for (int ek = 0; ek < CELL_PACK_L; ++ek)
                    F21_vo_tr_J[vj * CELL_PACK_L + ej][oi * CELL_PACK_L + ek] =
                        F21_mp[oi][vj][ek * CELL_PACK_L + ej] & 0x0F;
#endif
}

static inline void snova_cell_pack_L_nibble_pair(
    uint8_t L_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF][CELL_PACK_VTL_V_PAD],
    int j,
    const uint8_t *cell)
{
#if CELL_PACK_L == 3
    uint32_t r0 = ((uint32_t)cell[0]) | ((uint32_t)cell[1] <<  8) | ((uint32_t)cell[2] << 16);
    uint32_t r1 = ((uint32_t)cell[3]) | ((uint32_t)cell[4] <<  8) | ((uint32_t)cell[5] << 16);
    uint32_t r2 = ((uint32_t)cell[6]) | ((uint32_t)cell[7] <<  8) | ((uint32_t)cell[8] << 16);
    r0 &= 0x000F0F0F; r1 &= 0x000F0F0F; r2 &= 0x000F0F0F;
    uint32_t p01 = r0 | (r1 << 4);
    memcpy(&L_J_nibble_ma[0][j * 3], &p01, sizeof p01);
    memcpy(&L_J_nibble_ma[1][j * 3], &r2,  sizeof r2);
#elif CELL_PACK_L == 5
    const __m128i mask_lo128 = _mm_set1_epi8(0x0F);
    __m128i m_lo = _mm_loadu_si128((const __m128i *)cell);
    __m128i m_hi = _mm_loadu_si128((const __m128i *)(cell + 9));
    m_lo = _mm_and_si128(m_lo, mask_lo128);
    m_hi = _mm_and_si128(m_hi, mask_lo128);
    #define L5_ROW_LO(EI) _mm_shuffle_epi8(m_lo, _mm_setr_epi8(  \
        (EI)*5, (EI)*5+1, (EI)*5+2, (EI)*5+3, (EI)*5+4,          \
        (char)0x80, (char)0x80, (char)0x80, (char)0x80,           \
        (char)0x80, (char)0x80, (char)0x80, (char)0x80,           \
        (char)0x80, (char)0x80, (char)0x80))
    #define L5_ROW_HI(EI, BASE) _mm_shuffle_epi8(m_hi, _mm_setr_epi8(  \
        (BASE), (BASE)+1, (BASE)+2, (BASE)+3, (BASE)+4,                \
        (char)0x80, (char)0x80, (char)0x80, (char)0x80,                \
        (char)0x80, (char)0x80, (char)0x80, (char)0x80,                \
        (char)0x80, (char)0x80, (char)0x80))
    __m128i r0 = L5_ROW_LO(0), r1 = L5_ROW_LO(1), r2 = L5_ROW_LO(2);
    __m128i r3 = L5_ROW_HI(3, 6), r4 = L5_ROW_HI(4, 11);
    #undef L5_ROW_LO
    #undef L5_ROW_HI
    __m128i pair01 = _mm_or_si128(r0, _mm_slli_epi16(r1, 4));
    __m128i pair23 = _mm_or_si128(r2, _mm_slli_epi16(r3, 4));
    uint64_t v01 = (uint64_t)_mm_cvtsi128_si64(pair01);
    uint64_t v23 = (uint64_t)_mm_cvtsi128_si64(pair23);
    uint64_t v4  = (uint64_t)_mm_cvtsi128_si64(r4);
    memcpy(&L_J_nibble_ma[0][j * 5], &v01, sizeof v01);
    memcpy(&L_J_nibble_ma[1][j * 5], &v23, sizeof v23);
    memcpy(&L_J_nibble_ma[2][j * 5], &v4,  sizeof v4);
#else
    for (int r_pair = 0; r_pair < CELL_PACK_L_PAIR_HALF; ++r_pair) {
        int r_lo = 2 * r_pair, r_hi = 2 * r_pair + 1;
        for (int ej = 0; ej < CELL_PACK_L; ++ej) {
            uint8_t lo = cell[r_lo * CELL_PACK_L + ej] & 0x0F;
            uint8_t hi = cell[r_hi * CELL_PACK_L + ej] & 0x0F;
            L_J_nibble_ma[r_pair][j * CELL_PACK_L + ej] = (uint8_t)(lo | (hi << 4));
        }
    }
#if CELL_PACK_L_HAS_ODD
    for (int ej = 0; ej < CELL_PACK_L; ++ej)
        L_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF - 1][j * CELL_PACK_L + ej] =
            cell[CELL_PACK_L_FLOOR2 * CELL_PACK_L + ej] & 0x0F;
#endif
#endif
}

static inline void snova_cell_pack_R_col(
    uint8_t R_tr_J_ma[CELL_PACK_L][CELL_PACK_VTL_V_PAD],
    int k,
    const uint8_t *cell)
{
#if CELL_PACK_L == 3
    uint32_t c0 = ((uint32_t)cell[0]) | ((uint32_t)cell[3] <<  8) | ((uint32_t)cell[6] << 16);
    uint32_t c1 = ((uint32_t)cell[1]) | ((uint32_t)cell[4] <<  8) | ((uint32_t)cell[7] << 16);
    uint32_t c2 = ((uint32_t)cell[2]) | ((uint32_t)cell[5] <<  8) | ((uint32_t)cell[8] << 16);
    c0 &= 0x000F0F0F; c1 &= 0x000F0F0F; c2 &= 0x000F0F0F;
    memcpy(&R_tr_J_ma[0][k * 3], &c0, sizeof c0);
    memcpy(&R_tr_J_ma[1][k * 3], &c1, sizeof c1);
    memcpy(&R_tr_J_ma[2][k * 3], &c2, sizeof c2);
#elif CELL_PACK_L == 5
    const __m128i mask_lo128 = _mm_set1_epi8(0x0F);
    __m128i m_lo = _mm_loadu_si128((const __m128i *)cell);
    __m128i m_hi = _mm_loadu_si128((const __m128i *)(cell + 9));
    m_lo = _mm_and_si128(m_lo, mask_lo128);
    m_hi = _mm_and_si128(m_hi, mask_lo128);
    #define R5_COL(C) do {                                              \
        __m128i mask_lo_c = _mm_setr_epi8(                              \
            (C), (C) + 5, (C) + 10,                                     \
            (char)0x80, (char)0x80,                                     \
            (char)0x80, (char)0x80, (char)0x80, (char)0x80,             \
            (char)0x80, (char)0x80, (char)0x80, (char)0x80,             \
            (char)0x80, (char)0x80, (char)0x80);                        \
        __m128i mask_hi_c = _mm_setr_epi8(                              \
            (char)0x80, (char)0x80, (char)0x80,                         \
            (C) + 6, (C) + 11,                                          \
            (char)0x80, (char)0x80, (char)0x80, (char)0x80,             \
            (char)0x80, (char)0x80, (char)0x80, (char)0x80,             \
            (char)0x80, (char)0x80, (char)0x80);                        \
        __m128i col_lo = _mm_shuffle_epi8(m_lo, mask_lo_c);              \
        __m128i col_hi = _mm_shuffle_epi8(m_hi, mask_hi_c);              \
        __m128i col    = _mm_or_si128(col_lo, col_hi);                   \
        uint64_t cv = (uint64_t)_mm_cvtsi128_si64(col);                  \
        memcpy(&R_tr_J_ma[(C)][k * 5], &cv, sizeof cv);                  \
    } while (0)
    R5_COL(0); R5_COL(1); R5_COL(2); R5_COL(3); R5_COL(4);
    #undef R5_COL
#else
    for (int ek = 0; ek < CELL_PACK_L; ++ek)
        for (int c = 0; c < CELL_PACK_L; ++c)
            R_tr_J_ma[c][k * CELL_PACK_L + ek] = cell[ek * CELL_PACK_L + c] & 0x0F;
#endif
}

static inline void snova_cell_pack_LR_nibble_pair(
    uint8_t L_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF][CELL_PACK_VTL_V_PAD],
    uint8_t R_tr_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF][CELL_PACK_VTL_V_PAD],
    int j,
    const uint8_t *L_cell,
    const uint8_t *R_cell)
{
#if CELL_PACK_L == 3
    uint32_t L_r0 = ((uint32_t)L_cell[0]) | ((uint32_t)L_cell[1] <<  8) | ((uint32_t)L_cell[2] << 16);
    uint32_t L_r1 = ((uint32_t)L_cell[3]) | ((uint32_t)L_cell[4] <<  8) | ((uint32_t)L_cell[5] << 16);
    uint32_t L_r2 = ((uint32_t)L_cell[6]) | ((uint32_t)L_cell[7] <<  8) | ((uint32_t)L_cell[8] << 16);
    L_r0 &= 0x000F0F0F; L_r1 &= 0x000F0F0F; L_r2 &= 0x000F0F0F;
    uint32_t R_c0 = ((uint32_t)R_cell[0]) | ((uint32_t)R_cell[3] <<  8) | ((uint32_t)R_cell[6] << 16);
    uint32_t R_c1 = ((uint32_t)R_cell[1]) | ((uint32_t)R_cell[4] <<  8) | ((uint32_t)R_cell[7] << 16);
    uint32_t R_c2 = ((uint32_t)R_cell[2]) | ((uint32_t)R_cell[5] <<  8) | ((uint32_t)R_cell[8] << 16);
    R_c0 &= 0x000F0F0F; R_c1 &= 0x000F0F0F; R_c2 &= 0x000F0F0F;
    uint32_t L_p01 = L_r0 | (L_r1 << 4);
    uint32_t R_p01 = R_c0 | (R_c1 << 4);
    memcpy(&L_J_nibble_ma[0][j * 3],     &L_p01, sizeof L_p01);
    memcpy(&R_tr_J_nibble_ma[0][j * 3],  &R_p01, sizeof R_p01);
    memcpy(&L_J_nibble_ma[1][j * 3],     &L_r2,  sizeof L_r2);
    memcpy(&R_tr_J_nibble_ma[1][j * 3],  &R_c2,  sizeof R_c2);
#elif CELL_PACK_L == 5
    const __m128i mask_lo128 = _mm_set1_epi8(0x0F);
    __m128i Lm_lo = _mm_loadu_si128((const __m128i *)L_cell);
    __m128i Lm_hi = _mm_loadu_si128((const __m128i *)(L_cell + 9));
    __m128i Rm_lo = _mm_loadu_si128((const __m128i *)R_cell);
    __m128i Rm_hi = _mm_loadu_si128((const __m128i *)(R_cell + 9));
    Lm_lo = _mm_and_si128(Lm_lo, mask_lo128); Lm_hi = _mm_and_si128(Lm_hi, mask_lo128);
    Rm_lo = _mm_and_si128(Rm_lo, mask_lo128); Rm_hi = _mm_and_si128(Rm_hi, mask_lo128);
    #define L5C_ROW_LO(EI) _mm_shuffle_epi8(Lm_lo, _mm_setr_epi8(  \
        (EI)*5, (EI)*5+1, (EI)*5+2, (EI)*5+3, (EI)*5+4,            \
        (char)0x80, (char)0x80, (char)0x80, (char)0x80,            \
        (char)0x80, (char)0x80, (char)0x80, (char)0x80,            \
        (char)0x80, (char)0x80, (char)0x80))
    #define L5C_ROW_HI(EI, BASE) _mm_shuffle_epi8(Lm_hi, _mm_setr_epi8(  \
        (BASE), (BASE)+1, (BASE)+2, (BASE)+3, (BASE)+4,                  \
        (char)0x80, (char)0x80, (char)0x80, (char)0x80,                  \
        (char)0x80, (char)0x80, (char)0x80, (char)0x80,                  \
        (char)0x80, (char)0x80, (char)0x80))
    __m128i L_r0 = L5C_ROW_LO(0), L_r1 = L5C_ROW_LO(1), L_r2 = L5C_ROW_LO(2);
    __m128i L_r3 = L5C_ROW_HI(3, 6), L_r4 = L5C_ROW_HI(4, 11);
    #undef L5C_ROW_LO
    #undef L5C_ROW_HI
    #define R5C_COL(C) _mm_or_si128(                                     \
        _mm_shuffle_epi8(Rm_lo, _mm_setr_epi8(                            \
            (C), (C)+5, (C)+10,                                          \
            (char)0x80, (char)0x80,                                      \
            (char)0x80, (char)0x80, (char)0x80, (char)0x80,              \
            (char)0x80, (char)0x80, (char)0x80, (char)0x80,              \
            (char)0x80, (char)0x80, (char)0x80)),                        \
        _mm_shuffle_epi8(Rm_hi, _mm_setr_epi8(                            \
            (char)0x80, (char)0x80, (char)0x80,                          \
            (C)+6, (C)+11,                                               \
            (char)0x80, (char)0x80, (char)0x80, (char)0x80,              \
            (char)0x80, (char)0x80, (char)0x80, (char)0x80,              \
            (char)0x80, (char)0x80, (char)0x80)))
    __m128i R_c0 = R5C_COL(0), R_c1 = R5C_COL(1), R_c2 = R5C_COL(2);
    __m128i R_c3 = R5C_COL(3), R_c4 = R5C_COL(4);
    #undef R5C_COL
    __m128i L_p01 = _mm_or_si128(L_r0, _mm_slli_epi16(L_r1, 4));
    __m128i L_p23 = _mm_or_si128(L_r2, _mm_slli_epi16(L_r3, 4));
    __m128i R_p01 = _mm_or_si128(R_c0, _mm_slli_epi16(R_c1, 4));
    __m128i R_p23 = _mm_or_si128(R_c2, _mm_slli_epi16(R_c3, 4));
    uint64_t Lv01 = (uint64_t)_mm_cvtsi128_si64(L_p01);
    uint64_t Lv23 = (uint64_t)_mm_cvtsi128_si64(L_p23);
    uint64_t Lv4  = (uint64_t)_mm_cvtsi128_si64(L_r4);
    uint64_t Rv01 = (uint64_t)_mm_cvtsi128_si64(R_p01);
    uint64_t Rv23 = (uint64_t)_mm_cvtsi128_si64(R_p23);
    uint64_t Rv4  = (uint64_t)_mm_cvtsi128_si64(R_c4);
    memcpy(&L_J_nibble_ma[0][j * 5],    &Lv01, sizeof Lv01);
    memcpy(&L_J_nibble_ma[1][j * 5],    &Lv23, sizeof Lv23);
    memcpy(&L_J_nibble_ma[2][j * 5],    &Lv4,  sizeof Lv4);
    memcpy(&R_tr_J_nibble_ma[0][j * 5], &Rv01, sizeof Rv01);
    memcpy(&R_tr_J_nibble_ma[1][j * 5], &Rv23, sizeof Rv23);
    memcpy(&R_tr_J_nibble_ma[2][j * 5], &Rv4,  sizeof Rv4);
#else
    for (int r_pair = 0; r_pair < CELL_PACK_L_PAIR_HALF; ++r_pair) {
        int r_lo = 2 * r_pair, r_hi = 2 * r_pair + 1;
        for (int ej = 0; ej < CELL_PACK_L; ++ej) {
            L_J_nibble_ma[r_pair][j * CELL_PACK_L + ej] = (uint8_t)(
                (L_cell[r_lo * CELL_PACK_L + ej] & 0x0F) |
                ((L_cell[r_hi * CELL_PACK_L + ej] & 0x0F) << 4));
            R_tr_J_nibble_ma[r_pair][j * CELL_PACK_L + ej] = (uint8_t)(
                (R_cell[ej * CELL_PACK_L + r_lo] & 0x0F) |
                ((R_cell[ej * CELL_PACK_L + r_hi] & 0x0F) << 4));
        }
    }
#if CELL_PACK_L_HAS_ODD
    for (int ej = 0; ej < CELL_PACK_L; ++ej) {
        L_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF - 1][j * CELL_PACK_L + ej] =
            L_cell[CELL_PACK_L_FLOOR2 * CELL_PACK_L + ej] & 0x0F;
        R_tr_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF - 1][j * CELL_PACK_L + ej] =
            R_cell[ej * CELL_PACK_L + CELL_PACK_L_FLOOR2] & 0x0F;
    }
#endif
#endif
}

#endif
