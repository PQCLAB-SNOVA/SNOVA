/**
 * @file platforms/ref/cell_pack_ref.h
 */
#ifndef SNOVA_PLATFORMS_REF_CELL_PACK_REF_H
#define SNOVA_PLATFORMS_REF_CELL_PACK_REF_H

#include <stdint.h>
#include <string.h>
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

static inline void snova_cell_pack_jogress_LL(
    uint8_t *AJ, const gf16_t *A, int mainCol, int mainRow, int mainRow_pad)
{
    for (int mi = 0; mi < mainCol; ++mi)
        for (int mj = 0; mj < mainRow; ++mj)
            for (int ei = 0; ei < CELL_PACK_L; ++ei)
                for (int ej = 0; ej < CELL_PACK_L; ++ej)
                    AJ[(mi * CELL_PACK_L + ei) * mainRow_pad + (mj * CELL_PACK_L + ej)] =
                        A[(mi * mainRow + mj) * CELL_PACK_L * CELL_PACK_L + ei * CELL_PACK_L + ej] & 0x0F;
}

static inline void snova_cell_pack_jogress_Tr_F21vo(
    uint8_t F21_vo_tr_J[CELL_PACK_V * CELL_PACK_L][CELL_PACK_VTL_O_PAD],
    const gf16m_t F21_mp[CELL_PACK_O][CELL_PACK_V])
{
    for (int vj = 0; vj < CELL_PACK_V; ++vj)
        for (int oi = 0; oi < CELL_PACK_O; ++oi)
            for (int ej = 0; ej < CELL_PACK_L; ++ej)
                for (int ek = 0; ek < CELL_PACK_L; ++ek)
                    F21_vo_tr_J[vj * CELL_PACK_L + ej][oi * CELL_PACK_L + ek] =
                        F21_mp[oi][vj][ek * CELL_PACK_L + ej] & 0x0F;
}

static inline void snova_cell_pack_L_nibble_pair(
    uint8_t L_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF][CELL_PACK_VTL_V_PAD],
    int j,
    const uint8_t *cell)
{
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
}

static inline void snova_cell_pack_R_col(
    uint8_t R_tr_J_ma[CELL_PACK_L][CELL_PACK_VTL_V_PAD],
    int k,
    const uint8_t *cell)
{
    for (int ek = 0; ek < CELL_PACK_L; ++ek)
        for (int c = 0; c < CELL_PACK_L; ++c)
            R_tr_J_ma[c][k * CELL_PACK_L + ek] = cell[ek * CELL_PACK_L + c] & 0x0F;
}

static inline void snova_cell_pack_LR_nibble_pair(
    uint8_t L_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF][CELL_PACK_VTL_V_PAD],
    uint8_t R_tr_J_nibble_ma[CELL_PACK_RANK_PAIR_HALF][CELL_PACK_VTL_V_PAD],
    int j,
    const uint8_t *L_cell,
    const uint8_t *R_cell)
{
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
}

#endif
