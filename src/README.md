SNOVA
=======
This directory contains the official constant-time implementation of the SNOVA signature scheme.

Five implementations are provided:
1. The reference implementation in `snova_ref.c`,
2. An optimized implementation for odd prime $q$ in `snova_opt_q.c`. This version uses plain-C. For additional optimization on older compilers, it will use AVX2 instructions if available,
3. An optimized implementation for $q=16$ in `snova_opt_16.c`. This version uses plain-C. For additional optimization it will use AVX2 instructions if available,
4. Another optimized version for $q=16$ and $l=r$ in `snova_avx2_16.c`. This version uses AVX2 or ARM NEON instructions. For $l \neq 4$ this version is substantially faster than `snova_opt_16.c`.


Building
-------

Building SNOVA requires a C compiler and `make`. There are no other dependencies.
The SNOVA parameters are set in `snova_params.h`. The SNOVA parameters can also be changed by the command line parameters of the `make` command, e.g.
```
make clean all P="-D SNOVA_v=37 -D SNOVA_o=8 -D SNOVA_q=19 -D SNOVA_l=4"
```
An example command line build for $q=16$, using AES-CTR for the public key expansion, is
```
make clean all P="-D SNOVA_v=24 -D SNOVA_o=5 -D SNOVA_q=16 -D SNOVA_l=4 -D AESCTR"
```

Available optimization options are:
1. Use `make OPT=REF` to build the reference implementation.
2. Use `make OPT=OPT` for the optimized version.
3. Use `make OPT=AVX2` to use `snova_avx2_16.c` for $q=16$ and $l=2$. For other parameter sets `OPT=AVX2` is identical to `OPT=OPT`.


Symmetric Primitives
-------

The distribution comes with implementations of AES and SHAKE. It is also possible to use the OpenSSL library (version 3.3 or higher), which may be faster on non-AVX2 platforms.
To use the OpenSSL library build as
```
make clean all P="-D USE_OPENSSL" LIBS=-lcrypto
```
