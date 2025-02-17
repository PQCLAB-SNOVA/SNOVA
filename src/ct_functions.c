#include "ct_functions.h"
#include <string.h>

int ct_is_negative(int val)
{
    return ((val >> 31) & 1);
}

// Constant time version of: (val != 0)
uint32_t ct_gf16_is_not_zero(uint8_t val)
{
    return (val | (val >> 1) | (val >> 2) | (val >> 3)) & 1;
}

uint32_t ct_xgf16_is_not_zero(uint32_t val)
{
    return (val | (val >> 3) | (val >> 6) | (val >> 9)) & 1;
}

void snova_set_zero(void *ptr, size_t size)
{
    memset(ptr, 0, size);
}
