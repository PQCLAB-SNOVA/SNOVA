/**
 * @file secure_clear.h
 */
#ifndef SNOVA_SECURE_CLEAR_H
#define SNOVA_SECURE_CLEAR_H

#include <stddef.h>
#include <string.h>

#if defined(__GNUC__) || defined(__clang__)
static inline void snova_secure_clear(void *p, size_t n) {
    if (!p || !n) return;
    memset(p, 0, n);
    __asm__ __volatile__("" : : "r"(p) : "memory");
}
#elif defined(__STDC_LIB_EXT1__)
static inline void snova_secure_clear(void *p, size_t n) {
    if (p && n) memset_s(p, n, 0, n);
}
#elif defined(_WIN32)
#include <windows.h>
static inline void snova_secure_clear(void *p, size_t n) {
    if (p && n) SecureZeroMemory(p, n);
}
#else
static inline void snova_secure_clear(void *p, size_t n) {
    volatile unsigned char *vp = (volatile unsigned char *)p;
    while (n--) *vp++ = 0;
}
#endif

#define SNOVA_CLEAR(p, n) snova_secure_clear((p), (n))
#define SNOVA_CLEAR_OBJ(x) snova_secure_clear(&(x), sizeof(x))

#endif
