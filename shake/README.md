# Vectorized SHAKE XOF

This folder contains an implementation of SHAKE128 CTR-XOF in snova_shake.c.

The implementation can be tested using the provided test, e.g.
```
gcc -march=native -DTEST_SNOVA_SHAKE snova_shake.c ../shake/KeccakHash.c ../shake/SimpleFIPS202.c ../shake/KeccakP-1600-opt64.c ../shake/KeccakSponge.c ../shake/KeccakP-1600-times4-SIMD256.c -o test
./test
```

# XKCP source code

The distribution includes a copy of XKCP as source code. It was retrieved using:

```
git clone https://github.com/XKCP/XKCP
cd XKCP
git checkout 908262a95140e21598ca0a703e45d770a168e8ec
git submodule update --init
```

Then add the following to doc/HOWTO-customize.build
```
    <target name="FIPS202-timesx" inherits="FIPS202 K1600-plain-64bits-ua K1600x4-AVX2-ua K1600x8-AVX512-ua"/>
```
This will create a new target with the modules required by vexof. To extract the sources:
```
make FIPS202-timesx.pack
cd ..
tar -xf XKCP/bin/FIPS202-timesx.tar.gz
```

See LICENSE_XKCP for the XKCP license conditions.

----
This version used the following commit:

commit 908262a95140e21598ca0a703e45d770a168e8ec (HEAD -> master, origin/master, origin/HEAD)
Author: Hadi El-Yakhni <hadiyakhni33@gmail.com>
Date:   Sat Dec 23 17:20:52 2023 +0100

    rename chunk to output in the Xoodyak examples
