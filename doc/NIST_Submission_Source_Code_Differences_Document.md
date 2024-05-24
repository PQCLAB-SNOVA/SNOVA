NIST Submission Source Code Differences Document
=======

Optimized Implementation
-------
We have abandoned the previously submitted optimized implementation in favor of a new approach, which includes added AVX optimization. All optimized code can be found in this folder [path](https://github.com/PQCLAB-SNOVA/SNOVA/tree/main/snova_plasma)."

And optimization modes are configured through this [file](https://github.com/PQCLAB-SNOVA/SNOVA/blob/main/snova_plasma/snova_plasms_option.h).


OpenSSL clean
-------
Remove OpenSSL, replace it with [AES](https://github.com/PQCLAB-SNOVA/SNOVA/tree/main/aes) and [shake](https://github.com/PQCLAB-SNOVA/SNOVA/tree/main/shake) source code, making SNOVA a separate source code.

Constant-time
-------
Added some Constant-time methods to mitigate side-channel attacks.

KAT
-------
Except for the parameters (v, o, l): (37, 17, 2), (56, 25, 2), (75, 33, 2), which are new and differ from the KAT submitted to NIST, the remaining parameters' KATs are the same as those submitted to NIST.

The reason for the modifications to the new parameters: [https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/m11kg20sTyU/m/cLkGIDaiBAAJ?utm_medium=email&utm_source=footer](https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/m11kg20sTyU/m/cLkGIDaiBAAJ?utm_medium=email&utm_source=footer)

Variants
-------
An experimental variant using SHAKE for the public key expansion have been added. Shake can be configured via the PK_EXPAND_SHAKE setting in the makefile. (Note: This will generate different KATs.)

source code: snova_kernel.h -> line:114
```c
/**
 * pk expand from seed
 */
#if PK_EXPAND_SHAKE
void pk_expand(const uint8_t* pt_public_key_seed, uint8_t* out_pk) {
    uint64_t vexof_array[(bytes_prng_public + 7) / 8];
    vexof(pt_public_key_seed, 16, vexof_array, 8 * ((bytes_prng_public + 7) / 8));
    memcpy(out_pk, vexof_array, bytes_prng_public);
}
#else
void pk_expand(const uint8_t* pt_public_key_seed, uint8_t* out_pk) {
    hash_aes128(pt_public_key_seed, out_pk);
}
#endif
```