NIST Submission Source Code Differences Document
=======

Optimized Implementation
-------
We have abandoned the previously submitted optimized implementation in favor of a new approach, which includes added AVX optimization. All optimized code can be found in this folder [path](https://github.com/PQCLAB-SNOVA/SNOVA/tree/main/snova_plasma)."

And optimization modes are configured through this [file](https://github.com/PQCLAB-SNOVA/SNOVA/blob/main/snova_plasma/snova_plasma_option.h).


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

SHAKE
-------
A variant that uses SHAKE (Secure Hash Algorithm KECCAK) for the public key expansion has been added as an alternative to AES-CTR. We denote the XOF of this variant as SNOVA_SHAKE. SNOVA_SHAKE can be configured via the PK_EXPAND_SHAKE setting in the makefile. (Note: This will generate different KATs.)

The algorithm of SNOVA_SHAKE can be described by: extract the bytes from SHAKE128 in 168 bytes blocks where a block index is appended to the seed. The block size 168 follows from the rate of SHAKE128. SNOVA_SHAKE is specified by:

Let `SHAKE128(seed, n)` denote the n-th byte of the `SHAKE128` XOF (eXtendable Output Function) when instantiated with `seed` as input, and similarly `SNOVA_SHAKE(seed, n)`. Then for all required bytes:
```
SNOVA_SHAKE(seed, n) = SHAKE128(seed || floor(n / 168), (n mod 168))
```
where `floor(n / 168)` is the 8 bytes little-endian representation of `n / 168` rounded to below, `||` represents concatenation of bytes and `mod` is the integer modulus operation.
