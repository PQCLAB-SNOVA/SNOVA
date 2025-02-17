#ifndef CYCLES_H
#define CYCLES_H

#include <stdint.h>
#include <time.h>

double get_cpu_f(void);

#ifdef __ARM_ARCH
// cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_cur_freq
// sudo cpufreq-set -d 2.4GHz
static inline uint64_t get_cycles(void)
{
    struct timespec time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time);
    return (int64_t)((time.tv_sec * 1e9 + time.tv_nsec) * 2.4);
}
#else
static inline uint64_t get_cycles() {
    uint32_t lo, hi;
    uint64_t o;
    __asm__ __volatile__("rdtscp" : "=a"(lo), "=d"(hi) : : "%ecx");
    o = hi;
    o <<= 32;
    return (o | lo);
}
#endif

static inline uint64_t cpucycles_overhead(void) {
    uint64_t t0, t1, overhead = -1LL;
    unsigned int i;

    for (i = 0; i < 100000; i++) {
        t0 = get_cycles();
        __asm__ volatile("");
        t1 = get_cycles();
        if (t1 - t0 < overhead) overhead = t1 - t0;
    }

    return overhead;
}

#endif
