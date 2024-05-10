SNOVA
=======
This repository contains the latest official reference & AVX implementation of the SNOVA signature scheme in C language.

Please refer to this [document](https://github.com/pqclab-zero/SNOVA/blob/main/doc/NIST_Submission_Source_Code_Differences_Document.md) for the main differences between the previously submitted code to NIST and the current version.


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

1. test keypair, sign and verify
```bash
make clean
make test
```
2. test nist api
```bash
make clean
make test_api
```
3. test speed
```bash
make clean
make test_speed
```
```bash
make clean
make test_speed  TEST_SPEED_N=4096
```

Set the parameters
-------
Only need to modify the SNOVA_V, SNOVA_O, SNOVA_L, SK_IS_SEED in the Makefile file.

Example: (makefile line 37)
```make
SNOVA_V ?= 24
SNOVA_O ?= 5
SNOVA_L ?= 4

SK_IS_SEED ?= 0 # 0: sk = ssk; 1: sk = esk
```

Genarate KAT
-------
```bash
make clean
make PQCgenKAT
```

Tip: The following parameters can all be input during MAKE execution, such as
-------
```bash
make test_speed SNOVA_V=24 SNOVA_O=5 SNOVA_L=4 TEST_SPEED_N=4096
```

The following are the latest reference parameters.
-------

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

SL 5: 
| SNOVA_V | SNOVA_O | SNOVA_L |
|---------|---------|---------|
|      75 |      33 |       2 |
|      66 |      15 |       3 |
|      60 |      10 |       4 |


Example:
SL 5 (60, 10, 4)
```
make test_speed SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 TEST_SPEED_N=512
```
```
make PQCgenKAT SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 SK_IS_SEED=1
```
```
make test_api SNOVA_V=60 SNOVA_O=10 SNOVA_L=4 SK_IS_SEED=1
```


## Team Members

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

