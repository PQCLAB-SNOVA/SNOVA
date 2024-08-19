#include <stdio.h>

#define NO_MT4B
#include "../snova.h"
#include "../util/util.h"
#include "analysis.h"
#include "cycles.h"

#ifndef TEST_SPEED_N
#define TEST_SPEED_N 2048
#endif

// #define CPU_CLOCK 2294684000

double get_cpu_f(void) {
    FILE* fp;
    char buffer[1024];
    double cpu_frequency;

    fp = fopen("/proc/cpuinfo", "r");
    while (fgets(buffer, sizeof(buffer), fp)) {
        if (sscanf(buffer, "cpu MHz : %lf", &cpu_frequency) == 1) {
            break;
        }
    }

    printf("CPU Frequency: %.2f MHz\n", cpu_frequency);

    fclose(fp);
    return cpu_frequency * 1000000;
}

int main() {
    snova_init();
    double CPU_CLOCK = get_cpu_f();
    uint8_t array_digest[64];
    uint8_t array_signature1[bytes_signature + bytes_salt];
    uint8_t array_signature2[bytes_signature + bytes_salt];
    uint8_t seed[seed_length];
    uint8_t* pt_private_key_seed;
    uint8_t* pt_public_key_seed;
    uint8_t pk[bytes_pk], ssk[seed_length], esk[bytes_sk];
    uint8_t array_salt[bytes_salt];
    uint8_t entropy_input[48];

    uint64_t t0[TEST_SPEED_N * 2] = {0};
    uint64_t t1[TEST_SPEED_N * 2] = {0};
    uint64_t t2[TEST_SPEED_N * 2] = {0};
    uint64_t t3[TEST_SPEED_N * 2] = {0};
    uint64_t t4[TEST_SPEED_N * 2] = {0};
    int r = 0;

    printf("------------------------------------------------\n");
    printf("SNOVA params:\t\t\n");
    printf("(V, O, L) =\t\t(%d, %d, %d)\n", v_SNOVA, o_SNOVA, l_SNOVA);
#if PK_EXPAND_SHAKE
    printf("PK_EXPAND =\t\tPK_EXPAND_SHAKE\n");
#else
    printf("PK_EXPAND =\t\tAES\n");
#endif
    printf("PK size =\t\t%d\n", bytes_pk);
    printf("Sign size =\t\t%d\n", bytes_signature + bytes_salt);
    printf("CRYPTO_ALGNAME =\t%s\n", CRYPTO_ALGNAME);
    printf("------------------------------------------------\n");
    printf("SNOVA TEST SPEED N=\t%d\n", TEST_SPEED_N);
    printf("SNOVA TEST SPEED start...\n");
    printf("================================================\n");
    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }
    randombytes_init(entropy_input, NULL, 256);
    randombytes(seed, seed_length);

    pt_public_key_seed = seed;
    pt_private_key_seed = seed + seed_length_public;

    for (int i = 0; i < TEST_SPEED_N; i++) {
		randombytes(array_salt, bytes_salt);

        t0[i * 2] = get_cycles();
        generate_keys_ssk(pk, ssk, pt_public_key_seed, pt_private_key_seed);
        t0[i * 2 + 1] = get_cycles();

        t1[i * 2] = get_cycles();
        generate_keys_esk(pk, esk, pt_public_key_seed, pt_private_key_seed);
        t1[i * 2 + 1] = get_cycles();

        t2[i * 2] = get_cycles();
        sign_digest_ssk(array_signature2, array_digest, 64, array_salt, ssk);
        t2[i * 2 + 1] = get_cycles();

        t3[i * 2] = get_cycles();
        sign_digest_esk(array_signature1, array_digest, 64, array_salt, esk);
        t3[i * 2 + 1] = get_cycles();

        t4[i * 2] = get_cycles();
        r += verify_signture(array_digest, 64, array_signature2, pk);
        t4[i * 2 + 1] = get_cycles();
    }

    printf("generate_keys_ssk: \n");
    analysis(t0, TEST_SPEED_N, CPU_CLOCK);
    printf("------------------------------------------------\n");
    printf("generate_keys_esk: \n");
    analysis(t1, TEST_SPEED_N, CPU_CLOCK);
    printf("------------------------------------------------\n");
    printf("sign_digest_ssk: \n");
    analysis(t2, TEST_SPEED_N, CPU_CLOCK);
    printf("------------------------------------------------\n");
    printf("sign_digest_esk: \n");
    analysis(t3, TEST_SPEED_N, CPU_CLOCK);
    printf("------------------------------------------------\n");
    printf("verify_signture: \n");
    analysis(t4, TEST_SPEED_N, CPU_CLOCK);
    printf("================================================\n");
    printf("verify_signture fail count: %d / %d\n", -r, TEST_SPEED_N);
    printf("Total cycle: %lu, Time: %lfs\n", (t4[TEST_SPEED_N * 2 - 1] - t0[0]),
           (double)(t4[TEST_SPEED_N * 2 - 1] - t0[0]) / (double)CPU_CLOCK);
    printf("------------------------------------------------\n");
    return 0;
}
