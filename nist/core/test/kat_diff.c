/**
 * @file kat_diff.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../src/snova.h"
#include "../src/snova_core/drbg.h"

static long parse_hex_line(const char *line, const char *key, uint8_t *out) {
    size_t kl = strlen(key);
    if (strncmp(line, key, kl) != 0) return -1;
    const char *p = line + kl;
    while (*p == ' ' || *p == '=') p++;
    long n = 0;
    while (p[0] && p[1] && p[0] != '\n' && p[0] != '\r') {
        int hi, lo;
        char c = p[0];
        if (c >= '0' && c <= '9') hi = c - '0';
        else if (c >= 'A' && c <= 'F') hi = c - 'A' + 10;
        else if (c >= 'a' && c <= 'f') hi = c - 'a' + 10;
        else break;
        c = p[1];
        if (c >= '0' && c <= '9') lo = c - '0';
        else if (c >= 'A' && c <= 'F') lo = c - 'A' + 10;
        else if (c >= 'a' && c <= 'f') lo = c - 'a' + 10;
        else break;
        out[n++] = (uint8_t)((hi << 4) | lo);
        p += 2;
    }
    return n;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: kat_diff <PQCsignKAT_*.rsp> [count]\n");
        return 2;
    }
    const char *rsp = argv[1];
    int want = (argc > 2) ? atoi(argv[2]) : 3;
    FILE *f = fopen(rsp, "r");
    if (!f) { fprintf(stderr, "kat_diff: cannot open %s\n", rsp); return 1; }

    static uint8_t seed[64], msg[1 << 16], pk_e[1 << 20], sk_e[1 << 20],
        sm_e[1 << 20];
    static uint8_t pk[1 << 20], sk[1 << 20], sm[1 << 20], mo[1 << 16];
    long msgn = 0, pkn = 0, skn = 0, smn = 0;
    long smlen_e = 0;
    static char line[1 << 22];
    int count = -1, fails = 0, done = 0;

    while (done < want && fgets(line, sizeof line, f)) {
        long n;
        if (strncmp(line, "count = ", 8) == 0) count = atoi(line + 8);
        else if ((n = parse_hex_line(line, "seed", seed)) >= 0) (void)n;
        else if (strncmp(line, "mlen = ", 7) == 0) {   }
        else if ((n = parse_hex_line(line, "msg", msg)) >= 0) msgn = n;
        else if ((n = parse_hex_line(line, "pk", pk_e)) >= 0) pkn = n;
        else if ((n = parse_hex_line(line, "sk", sk_e)) >= 0) skn = n;
        else if (strncmp(line, "smlen = ", 8) == 0) smlen_e = atol(line + 8);
        else if ((n = parse_hex_line(line, "sm", sm_e)) >= 0) {
            smn = n;
            randombytes_init(seed, NULL, 256);
            size_t smlen = 0, mlen = 0;
            snova_crypto_sign_keypair(pk, sk);
            snova_crypto_sign(sm, &smlen, msg, (size_t)msgn, sk);
            int open_ok =
                snova_crypto_sign_open(mo, &mlen, sm, smlen, pk) == 0 &&
                mlen == (size_t)msgn && memcmp(mo, msg, msgn) == 0;

            int ok_pk = (pkn > 0 && memcmp(pk, pk_e, pkn) == 0);
            int ok_sk = (skn > 0 && memcmp(sk, sk_e, skn) == 0);
            int ok_sm = ((long)smlen == smlen_e && smn == smlen_e &&
                         memcmp(sm, sm_e, smn) == 0);
            printf("count=%d  pk:%s  sk:%s  sm:%s  open:%s\n", count,
                   ok_pk ? "OK" : "DIFF", ok_sk ? "OK" : "DIFF",
                   ok_sm ? "OK" : "DIFF", open_ok ? "OK" : "DIFF");
            if (!ok_pk) {
                long d = 0;
                while (d < pkn && pk[d] == pk_e[d]) d++;
                printf("  pk first diff @%ld got %02X exp %02X (pkn=%ld)\n", d,
                       pk[d], pk_e[d], pkn);
            }
            if (ok_pk && !ok_sm) {
                long d = 0;
                while (d < smn && sm[d] == sm_e[d]) d++;
                printf("  sm first diff @%ld got %02X exp %02X (smlen=%zu exp %ld)\n",
                       d, sm[d], sm_e[d], smlen, smlen_e);
            }
            if (!ok_pk || !ok_sk || !ok_sm || !open_ok) fails++;
            done++;
        }
    }
    fclose(f);
    if (fails == 0) { printf("== KAT %d/%d MATCH ==\n", done, done); return 0; }
    printf("== KAT %d FAIL(S) of %d ==\n", fails, done);
    return 1;
}
