SNOVA
=======
This directory contains a tool to create the recommended SNOVA Round 3 instances, as well as the short signature alternatives.

Use `make` to create `ref`, `opt`, `avx2`, and `mem` directories.

In one of those directories use e.g.
```
make clean kat speed
make speed
make digest
```

Building SNOVA requires a C compiler and `make`.


# KAT digests

This directory contains a file `KATs` with 24 byte SHAKE256 digests. It is the output of `make digest` after building.

The official SNOVA Round 3 KAT files can be found online in the repository https://github.com/PQCLAB-SNOVA/SNOVA_KAT.
