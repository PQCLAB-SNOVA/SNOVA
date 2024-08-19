/**
 * Constant time functions and miscellaneous.
 * Put in a separate file to prevent compiler over-optimisation.
 *
 * See "Compiler-introduced timing leak in Kyber reference implementation"
 * https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/hqbtIGFKIpU
 */

#ifndef CTFUN_H
#define CTFUN_H

#include <stdint.h>
#include <stddef.h>

int ct_is_negative(int val);
uint32_t ct_xgf16_is_not_zero(uint32_t val);
uint32_t ct_gf16_is_not_zero(uint8_t val);

void snova_set_zero(void *ptr, size_t size);
#define SNOVA_CLEAR(x) snova_set_zero(x, sizeof(x));

#define SNOVA_CLEAR_BYTE(x, byte) snova_set_zero(x, byte);

#endif
