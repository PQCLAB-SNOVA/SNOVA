SNOVA
=======
This repository contains the latest official reference & optimisation & AVX2 implementation of the SNOVA signature scheme in C language.
The offical repository is at https://github.com/PQCLAB-SNOVA/SNOVA

There is only a single version of the source code of SNOVA. The optimization level and the parameter sets are controlled in the makefile or by the command line parameters of the `make` command.

The main parts of Optimized (including AVX2) are in the `snova_plasma` folder.

Frequently Used Table of Contents
-------
- [SNOVA parameter](#snova-parameter)
- [Test programs](#test-programs)
- [Set the parameters](#set-the-parameters)
- [Optimization Settings](#optimization-settings)

SNOVA parameter
-------
| SL |         Name  |  V |  O |  L | sk size (esk) | sk size (ssk) |    pk size   | sign size  |
|----| --------------|----|----|----|---------------|---------------|--------------|------------|
|  1 | SNOVA_37_17_2 | 37 | 17 |  2 |         91440 |            48 |         9842 |        124 |
|  1 |  SNOVA_25_8_3 | 25 |  8 |  3 |         39576 |            48 |         2320 |      164.5 |
|  1 |  SNOVA_24_5_4 | 24 |  5 |  4 |         36848 |            48 |         1016 |        248 |
|  3 | SNOVA_56_25_2 | 56 | 25 |  2 |        300848 |            48 |        31266 |        178 |
|  3 | SNOVA_49_11_3 | 49 | 11 |  3 |        177060 |            48 |       6005.5 |        286 |
|  3 |  SNOVA_37_8_4 | 37 |  8 |  4 |        133040 |            48 |         4112 |        376 |
|  3 |  SNOVA_24_5_5 | 24 |  5 |  5 |         60048 |            48 |       1578.5 |      378.5 |
|  5 | SNOVA_75_33_2 | 75 | 33 |  2 |        704532 |            48 |        71890 |        232 |
|  5 | SNOVA_66_15_3 | 66 | 15 |  3 |        435423 |            48 |      15203.5 |      380.5 |
|  5 | SNOVA_60_10_4 | 60 | 10 |  4 |        395248 |            48 |         8016 |        576 |
|  5 |  SNOVA_29_6_5 | 29 |  6 |  5 |        100398 |            48 |         2716 |      453.5 |

Build instructions
-------
These implementations include multiple test programs and an easy-to-compile Makefile.

Prerequisites for compiling this repository
-------
None.
This repository contains the symmetric primitive libraries as source code. See the respective files for credits and copyright notices.

Test programs
-------
The steps to compile and run the SNOVA signature scheme's test program on Linux are as follows:

### Test keypair, sign and verify:
```bash
make clean
make test
```
### Test nist api:
```bash
make clean
make test_api
```
### Generate KAT:
```bash
make clean
make PQCgenKAT
```
### Test speed:
```bash
make clean
make test_speed
```
```bash
make clean
make test_speed  TEST_SPEED_N=4096
```
Tips. TEST_SPEED_N specifies the number of test iterations.


Set the parameters
-------
Only need to modify the SNOVA_V, SNOVA_O, SNOVA_L, SK_IS_SEED, PK_EXPAND_SHAKE in the Makefile file.

Example: (makefile line 37)
```bash
SNOVA_V ?= 24
SNOVA_O ?= 5
SNOVA_L ?= 4

# 0: sk = ssk; 1: sk = esk
SK_IS_SEED ?= 0

# 0: disable; 1: enable
PK_EXPAND_SHAKE ?= 0
```

Tip: The following parameters can all be input during MAKE execution, such as
```bash
make test_speed SNOVA_V=24 SNOVA_O=5 SNOVA_L=4 SK_IS_SEED=1 PK_EXPAND_SHAKE=1 TEST_SPEED_N=4096
```

### The following are the latest reference parameters.

SL 1: 
| SNOVA_V | SNOVA_O | SNOVA_L |
|---------|---------|---------|
|      37 |      17 |       2 |
|      25 |       8 |       3 |
|      24 |       5 |       4 |

SL 3: 
| SNOVA_V | SNOVA_O | SNOVA_L |
|---------|---------|---------|
|      56 |      25 |       2 |
|      49 |      11 |       3 |
|      37 |       8 |       4 |
|      24 |       5 |       5 |

SL 5: 
| SNOVA_V | SNOVA_O | SNOVA_L |
|---------|---------|---------|
|      75 |      33 |       2 |
|      66 |      15 |       3 |
|      60 |      10 |       4 |
|      29 |       6 |       5 |


### Example: SL 5 (60, 10, 4) 
```bash
make test_speed SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 TEST_SPEED_N=512
```
```bash
make PQCgenKAT SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 SK_IS_SEED=1
```
```bash
make test_api SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 SK_IS_SEED=0
```

Optimization Source Code
-------
All optimized code can be found in this folder "snova_plasma", You can refer to line 6 in "snova.c" to see how to switch optimization methods.

Optimization Settings
-------
To configure optimization settings, you only need to adjust the OPTIMISATION parameter in the makefile or use the "OPTIMISATION=" command. A detailed description can also be found within the makefile.

If no configuration is provided, the most optimal method will be used automatically.

### Tips.
```bash
OPTIMISATION = 0  # Using Reference
OPTIMISATION = 1  # Using General Optimization
OPTIMISATION = 2  # Using AVX2 Optimization
```

### Example: (Using Reference)
```bash
make test_api SNOVA_V=24 SNOVA_O=5 SNOVA_L=4 SK_IS_SEED=1 OPTIMISATION=0
```
### Example: (Using General Optimization)
```bash
make test_api SNOVA_V=24 SNOVA_O=5 SNOVA_L=4 SK_IS_SEED=1 OPTIMISATION=1
```
Example: (Using AVX2 Optimization)
```bash
make test_api SNOVA_V=24 SNOVA_O=5 SNOVA_L=4 SK_IS_SEED=1 OPTIMISATION=2
```

SHAKE
-------
A variant that uses SHAKE128 for the public key expansion has been added as an alternative to AES-CTR. SHAKE can be configured via the PK_EXPAND_SHAKE setting in the makefile. (Note: This will generate different KATs.)

### Example: (Using PK_EXPAND_SHAKE)
```bash
make test_api SNOVA_V=24 SNOVA_O=5 SNOVA_L=4 SK_IS_SEED=1 PK_EXPAND_SHAKE=1
```

Team Members
-------
Thank you to our team members from SNOVA:

- Lih-Chung Wang
- Chun-Yen Chou
- Jintai Ding
- Yen-Liang Kuan
- Jan Adriaan Leegwater
- Ming-Siou Li
- Bo-Shu Tseng
- Po-En Tseng
- Chia-Chun Wang
