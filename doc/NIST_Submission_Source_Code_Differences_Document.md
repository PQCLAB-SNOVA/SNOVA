NIST Submission Source Code Differences Document
=======

Optimized Implementation
-------
We have abandoned the previously submitted optimized implementation in favor of a new approach, which includes added AVX optimization. All optimized code can be found in this folder [path](https://github.com/pqclab-zero/SNOVA/tree/main/snova_plasma)."

And optimization modes are configured through this [file](https://github.com/pqclab-zero/SNOVA/blob/main/snova_plasma/snova_plasms_option.h).


OpenSSL clean
-------
Remove OpenSSL, replace it with [AES](https://github.com/pqclab-zero/SNOVA/tree/main/aes) and [shake](https://github.com/pqclab-zero/SNOVA/tree/main/shake) source code, making SNOVA a separate source code.

Constant-time
-------
Added some Constant-time methods to mitigate side-channel attacks.

KAT
-------
Except for the parameters (v, o, l): (37, 17, 2), (56, 25, 2), (75, 33, 2), which are new and differ from the KAT submitted to NIST, the remaining parameters' KATs are the same as those submitted to NIST.

The reason for the modifications to the new parameters: [https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/m11kg20sTyU/m/cLkGIDaiBAAJ?utm_medium=email&utm_source=footer](https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/m11kg20sTyU/m/cLkGIDaiBAAJ?utm_medium=email&utm_source=footer)

variants
-------
New experimental variants Vexof have been added, allowing different methods for expanding the public key. Shake can be configured via the PK_EXPAND_VEXOF setting in the makefile. (Note: This will generate different KATs.)

source code: snova_kernel.h -> line:114
```c
/**
 * pk expand from seed
 */
#if PK_EXPAND_VEXOF
#include "vexof/vexof.h"
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