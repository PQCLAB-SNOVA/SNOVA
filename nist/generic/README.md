SNOVA
=======
This directory contains the official constant-time optimized implementation of the SNOVA signature scheme.

The following implementations are provided:
1. The reference implementation in `snova_ref.c`,
2. Optimized implementations in `snova_opt_*.c`. These versions use GFNI if available or else AVX2. If AVX2 is not available, plain-C is used.


Building
-------

Building SNOVA requires a C compiler and `make`. There are no other dependencies.
The SNOVA parameters are set in `snova_params.h`. The SNOVA parameters can also be changed by the command line parameters of the `make` command, e.g.
```
make clean all P="-D SNOVA_v=22 -D SNOVA_o=4 -D SNOVA_q=23 -D SNOVA_l=4"
```
An example command line build for $q=16$, using AES-CTR for the public key expansion, is
```
make clean all P="-D SNOVA_v=27 -D SNOVA_o=5 -D SNOVA_q=16 -D SNOVA_l=4 -D AESCTR"
```

Available optimization options are:
1. Use `make OPT=REF` to build the reference implementation.
2. Use `make OPT=OPT` (default) for the optimized version using explicit AVX2 or GFNI instructions if available. The file `snova_opt.c` will include the appropriate actual implementation.


Symmetric Primitives
-------

The distribution comes with implementations of AES and SHAKE. It is also possible to use the OpenSSL library (version 3.3 or higher), which may be faster on non-AVX2 platforms.
To use the OpenSSL library build as
```
make clean all P="-D USE_OPENSSL" LIBS=-lcrypto
```


Compatibility with Round 2 SNOVA
-------

While the SNOVA parameter space has been expanded and the recommended parameters have changed, the algorithm itself has not changed between rounds 2 and 3. The Round 2 KAT files can still be created by using appropriate parameters. For example
```
make clean kat P="-D SNOVA_v=24 -D SNOVA_o=5 -D SNOVA_q=16 -D SNOVA_l=4 -D SNOVA_r=4 -D SNOVA_m1=5 -D SNOVA_alpha=20 -D FIXED_ABQ=0 -D HASH_PK=0 -D ROUND2_T12=1 -D SNOVA_NAME=SNOVA_24_5_4_SHAKE"
```
This will create a response file `PQCsignKAT_SNOVA_24_5_4_SHAKE.rsp`. To check for changes against the Round 2 version:
```
diff PQCsignKAT_SNOVA_24_5_4_SHAKE.rsp $SNOVA_KAT/Round2/PQCsignKAT_SNOVA_24_5_4_SHAKE_SSK.rsp
```
