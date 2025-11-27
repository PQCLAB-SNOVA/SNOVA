SNOVA
=======
This directory contains the official constant-time implementation of the SNOVA signature scheme.


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
2. Use `make OPT=OPT` for the plain-C optimized version.
3. Use `make OPT=AVX2` for a version using AVX2 vectorization instructions, if available.

The distribution comes with implementations of AES and SHAKE. It is also possible to use the OpenSSL library (version 3.3 or higher), which may be faster on non-AVX2 platforms.
To use the OpenSSL library build as
```
make clean all P="-D USE_OPENSSL" LIBS=-lcrypto
```

Recommended parameters
-------

Our recommended parameter sets for odd prime $q$ all have a matrix rank $l=4$. The following are the recommended parameters:

| SL |        Name      |  v |  o |   q |  l |  sk size |  pk size |  sign size |
|----|------------------|----|----|-----|----|----------|----------|------------|
|  1 |  SNOVA_24_5_23_4 | 24 |  5 |  23 |  4 |      48  |      616 |        282 |
|  1 |  SNOVA_24_5_16_4 | 24 |  5 |  16 |  4 |      48  |     1016 |        248 |
|  1 | SNOVA_43_17_16_2 | 43 | 17 |  16 |  2 |      48  |     9842 |        136 |
|  3 |  SNOVA_37_8_19_4 | 37 |  8 |  19 |  4 |      48  |     2269 |        400 |
|  3 |  SNOVA_37_8_16_4 | 37 |  8 |  16 |  4 |      48  |     4112 |        376 |
|  3 | SNOVA_69_25_16_2 | 69 | 25 |  16 |  2 |      48  |    31266 |        204 |
|  5 | SNOVA_60_10_23_4 | 60 | 10 |  23 |  4 |      48  |     4702 |        656 |
|  5 | SNOVA_60_10_16_4 | 60 | 10 |  16 |  4 |      48  |     8016 |        576 |
|  5 | SNOVA_99_25_16_2 | 99 | 25 |  16 |  2 |      48  |    71890 |        280 |


Performance
-------

The performance of SNOVA depends significantly on the version of the compiler used. The best performance was obtained using gcc version 15.2.1 20250813 on Arch Linux.
We have observed the following cycle counts for the vectorized `OPT=AVX2` version of SNOVA_24_5_23_4:

| Compiler | version |   Genkey  |     Sign   |  Verify  |
|----------|---------|-----------|------------|----------|
|   gcc    |  15.2.1 |   404,634 |    692,838 | 317,819 |
|   gcc    |  13.3.0 |   405,211 |    900,276 | 322,433 |
|   gcc    |  11.4.1 |   448,117 |  1,013,698 | 340,407 |
|  clang   |  20.1.8 |   893,024 |    970,647 | 377,861 |

These cycle counts were collected on an Intel(R) Core(TM) Ultra 7 155H (Meteor Lake) laptop running Arch Linux. We report the median over 2048 benchmark runs.
