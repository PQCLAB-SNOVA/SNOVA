SNOVA
=======
This repository contains the latest official Reference, Optimized, and AVX2 implementations in C of the SNOVA signature scheme.

There is only a single version of the source code of SNOVA for q=16. The parameter sets and the optimization level are controlled in the makefile or by the command line parameters of the `make` command  in the `./src` subdirectory. The makefile in the top level directory will create subdirectories corresponding to the various SNOVA parameter sets.

The `oddqsrc` directory contains a version of SNOVA for odd primes in the range q=7...31. This directory also contains an improved and constant-time version for q=16, l=4. See the README.md in the `oddqsrc` directory.

## Build instructions

Building requires `gcc` and `make`. There are no other dependencies.

To create the various SNOVA instances in separate directories enter `make` in this folder.  This will also generate the known-answer tests, the KAT files. Enter `make clean` to delete all generated folders and files contained therein. For more options and details on the SNOVA parameters see `./src/readme.md`.

## Options

The `platform` option sets the optimization level. Supported values are `ref` (default), `opt`, and `avx2`. 
For example: `make platform=avx2`
