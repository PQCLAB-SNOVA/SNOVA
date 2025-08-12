SNOVA
=======
This directory contains a constant-time implementation of the SNOVA signature scheme in C language which allows the field order $q$ to be an odd prime in the range 7...31. In addition, this directory contains a new optimized and constant-time implementation for $q=16, l=4$.


Building
-------

Building SNOVA requires a C compiler, `make` and the OpenSSL library. The SNOVA parameters are set in `snova_params.h`. The SNOVA parameters can also be changed by the command line parameters of the `make` command, e.g.
```
make clean all P="-D SNOVA_v=37 -D SNOVA_o=8 -D SNOVA_q=19 -D SNOVA_l=4"
```
An example command line build for $q=16$, using AES-CTR for the public key expansion, is
```
make clean all P="-D SNOVA_v=24 -D SNOVA_o=5 -D SNOVA_q=16 -D SNOVA_l=4 -D AESCTR"
```

Available optimization options are:
1. Use `make OPT=REF` to build the reference implementation.
2. Use `make OPT=OPT` (default) for the plain-C optimized version.


Recommended parameters
-------

Our recommended parameter sets for odd prime $q$ all feature a matrix rank $l=4$. The following are the preliminary recommended parameters:

| SL |        Name      |  V |  O |   q |  L |  sk size |  pk size | sign size  |
|----|------------------|----|----|-----|----|----------|----------|------------|
|  1 |  SNOVA_24_5_23_4 | 24 |  5 |  23 |  4 |      48  |      616 |        282 |
|  3 |  SNOVA_37_8_19_4 | 37 |  8 |  19 |  4 |      48  |     2269 |        400 |
|  5 | SNOVA_60_10_23_4 | 60 | 10 |  23 |  4 |      48  |     4702 |        656 |


Performance
-------

While the optimized versions uses only C statements, the compiler will actually vectorize the code. We found that the level of vectorization that the compiler produces depends significantly on the compiler used, and also the version of the compiler used. The best performance was obtained using gcc version 15.1.1 20250729 on Arch Linux.

We have observed the following cycle counts for SNOVA_24_5_23_4:

| Compiler | version |   Genkey  |    Sign   | Verify  |
|----------|---------|-----------|-----------|---------|
|   gcc    |  15.1.1 |   455,475 |   775,198 | 363,258 |
|   gcc    |  13.3.0 |   462,056 |   954,194 | 852,167 |
|   gcc    |  11.4.1 |   508,161 | 1,445,794 | 769,554 |
|  clang   |  20.1.8 |   952,002 | 1,955,005 | 559,704 |

These cycle counts were collected on an Intel(R) Core(TM) Ultra 7 155H (Meteor Lake) laptop running Arch Linux and using gcc version 15.1.1 20250729. We report the median over 2048 benchmark runs.

The performance on older compilers can be made to be close to that of gcc 15.1.1 by using explicit AVX2 vectorization instructions. We have chosen not to create a version with explicit AVX2 vectorization at this point in time. This would be more appropriate when preparing SNOVA for use in production.
