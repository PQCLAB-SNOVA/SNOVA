SNOVA
=======
This directory contains a tool to create the recommended SNOVA 2.2 instances.

Use `make` to create `ref`, `opt`, and `avx2` directories.

In one of those directories use e.g.
```
make clean kat speed
make speed
make digest
```

Building SNOVA requires a C compiler and `make`.


# KAT digests

This directory contains KAT digests for SNOVA 2.2. It is the output of `make digest` after building.

The official KAT files can be found in the repository https://github.com/PQCLAB-SNOVA/SNOVA_KAT. The KAT files for q=16, l=4 have not been changed since Round 2.
