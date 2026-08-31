# Vectorized SHAKE XOF

This folder contains implementations of SHAKE128 CTR-XOF in snova_shake_ref.c and  snova_shake_opt.c.
The optimised implementation can be tested against the reference using the genkat function of SNOVA
for optimisation levels OPTIMISATION=0 versus OPTIMISATION=1 or OPTIMISATION=2.

# XKCP source code

The distribution includes a copy of XKCP as source code. It was retrieved using:

```
git clone https://github.com/XKCP/XKCP
cd XKCP
git checkout ade40f8e46298ba3cab7af2f2e5ec739f6160407
git submodule update --init
```

Then add the following to doc/HOWTO-customize.build
```
    <target name="snova_shake" inherits="FIPS202 K1600-plain-64bits-ua K1600x4-AVX2-ua K1600x8-AVX512-ua"/>
```
This will create a new target with the modules required by SNOVA. To extract the sources:
```
make snova_shake.pack
cd ..
tar -xf XKCP/bin/snova_shake.tar.gz
```

See LICENSE_XKCP for the XKCP license conditions.

----
This version used the following commit:

commit ade40f8e46298ba3cab7af2f2e5ec739f6160407 (grafted, HEAD -> master, origin/master, origin/HEAD)
Author: Ryad Benadjila <ryadbenadjila@gmail.com>
Date:   Wed Jul 3 08:26:09 2024 +0200

    Fix Xoodoo API discrepancy between declaration and definition (should fix issue #147).