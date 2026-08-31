/**
 * @file gf16_spec.h
 */
#ifndef GF16_SPEC_H
#define GF16_SPEC_H

#define GF16_REDUCTION_POLY 0x13u

#define GF16_GENERATOR 0x2u

#define GF16_ORDER 16
#define GF16_MULORDER 15

#define GF16_S_ENTRY(i, j) ((unsigned char)(8 - ((i) + (j))))

#endif
