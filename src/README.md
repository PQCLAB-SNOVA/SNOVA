SNOVA
=======
This directory contains the official constant-time optimized implementation of the SNOVA signature scheme.


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
1. Use `make OPT=REF` to build the reference implementation in `snova_ref.c`,
2. Use `make OPT=OPT` (default) for the optimized version.
3. Use `make OPT=OPT2` for version specifically optimized for $l=2$. On x86 this is the fastest version when $l=2$.
4. Use `make OPT=AVX2` for a faster optimized version for $l=4,5$ that uses explicit AVX2 or GFNI instructions.
5. Use `make OPT=MEM` for a plain-C version that uses substantially less memory. On x86 this version is about three to five times slower than the `AVX2` and `OPT2` versions. We expect that further reductions in both compute time and memory usage are possible.


Symmetric Primitives
-------

The distribution comes with implementations of AES and SHAKE. It is also possible to use the AES implementation in the OpenSSL library, which may be faster on non-AVX2 platforms.
To use the OpenSSL library for AES, build as
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
