#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "cycles.h"

void analysis(uint64_t *t, uint32_t len, const uint64_t CPU_CLOCK);

static int cmp_uint64(const void *a, const void *b) {
    if (*(uint64_t *)a < *(uint64_t *)b) return -1;
    if (*(uint64_t *)a > *(uint64_t *)b) return 1;
    return 0;
}

static uint64_t median(uint64_t *l, size_t len) {
    qsort(l, len, sizeof(uint64_t), cmp_uint64);

    if (len % 2)
        return l[len / 2];
    else
        return (l[len / 2 - 1] + l[len / 2]) / 2;
}

static uint64_t average(uint64_t *t, size_t len) {
    uint64_t acc = 0;
    for (size_t i = 0; i < len; i++) {
        acc += t[i];
    }
    return acc / len;
}

void analysis(uint64_t *t, uint32_t len, const uint64_t CPU_CLOCK) {
    uint64_t *td = malloc(sizeof(uint64_t) * len);

    for (uint32_t i = 0; i < len; ++i) {
        td[i] = t[i * 2 + 1] - t[i * 2];
    }

    uint64_t md = median(td, len);
    uint64_t avg = average(td, len);
    printf("median cycle:\t%lu cycles/ticks\n", md);
    printf("average cycle:\t%lu cycles/ticks\n", avg);
    printf("median:\t\t%lf ticks/sec\n", (double)CPU_CLOCK / (double)md);
    printf("average:\t%lf ticks/sec\n", (double)CPU_CLOCK / (double)avg);

    free(td);
}

#endif