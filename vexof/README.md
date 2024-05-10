# VeXOF

Implementation of SHAKE128 CTR-XOF

Test using
```
gcc -march=skylake -DDEBUG test.c vexof.c ../shake/KeccakHash.c ../shake/SimpleFIPS202.c ../shake/KeccakP-1600-opt64.c ../shake/KeccakSponge.c ../shake/KeccakP-1600-times4-SIMD256.c -o test
./test
```
