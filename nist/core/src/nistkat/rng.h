/**
 * @file nistkat/rng.h
 */
#ifndef SNOVA_NISTKAT_RNG_H
#define SNOVA_NISTKAT_RNG_H

#include "../snova_core/drbg.h"

typedef struct {
    unsigned char buffer[16];
    int buffer_pos;
    unsigned long length_remaining;
    unsigned char key[32];
    unsigned char ctr[16];
} AES_XOF_struct;

#define RNG_SUCCESS 0
#define RNG_BAD_MAXLEN -1
#define RNG_BAD_OUTBUF -2
#define RNG_BAD_REQ_LEN -3

#endif
