ARCH=native

# Autodetect AVX
AVX2 := $(findstring AVX2, $(shell gcc -march=$(ARCH) -dM -E - < /dev/null))
AVX512 := $(findstring AVX512, $(shell gcc -march=$(ARCH) -dM -E - < /dev/null))

# Autodetect AES support on ARM
ARM_FEATURE_AES := $(findstring ARM_FEATURE_AES, $(shell gcc -march=$(ARCH) -dM -E - < /dev/null))

CC = gcc
# CC = clang

ifeq ($(ARM_FEATURE_AES), ARM_FEATURE_AES)
CFLAGS = -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -ftree-vectorize -O3 -march=armv8-a+crypto -mtune=$(ARCH)
else
CFLAGS = -std=c99 -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -ftree-vectorize -O3 -march=$(ARCH) -mtune=$(ARCH)
endif

BUILD_OUT_PATH = ./build/

OLIST = $(BUILD_OUT_PATH)rng.o $(BUILD_OUT_PATH)snova.o

SYMMETRICLIB = shake/KeccakHash.c shake/SimpleFIPS202.c shake/KeccakP-1600-opt64.c shake/KeccakSponge.c shake/snova_shake.c
SYMMETRICLIB += aes/aes_c.c aes/snova_aes.c
ifeq ($(AVX512), AVX512)
SYMMETRICLIB += shake/KeccakP-1600-times4-SIMD512.c shake/KeccakP-1600-times8-SIMD512.c aes/aes128_ni.c aes/aes256_ni.c
else
ifeq ($(AVX2), AVX2)
SYMMETRICLIB += shake/KeccakP-1600-times4-SIMD256.c aes/aes128_ni.c aes/aes256_ni.c
else
ifeq ($(ARM_FEATURE_AES), ARM_FEATURE_AES)
SYMMETRICLIB += aes/aes128_armv8.c aes/aes256_armv8.c
endif
endif
endif
SYMMETRICLIBO = $(SYMMETRICLIB:.c=.o)
STATICLIB = build/libsnovasym.a

LIBS = -L$(BUILD_OUT_PATH) -lsnovasym

# snova params
SNOVA_V ?= 24
SNOVA_O ?= 5
SNOVA_L ?= 4

# 0: sk = esk; 1: sk = ssk
SK_IS_SEED ?= 0

# 0: disable; 1: enable
PK_EXPAND_SHAKE ?= 0

# 0: Reference, 1:Optimised, 2: AVX2
ifeq ($(AVX2), AVX2)
OPTIMISATION ?= 2
else
OPTIMISATION ?= 1
endif

ifeq ($(PK_EXPAND_SHAKE), 1)
    PK_EXPAND = _SHAKE
else
    PK_EXPAND = 
endif

CRYPTO_ALGNAME = SNOVA_$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L)$(PK_EXPAND)
SNOVA_PARAMS = -D v_SNOVA=$(SNOVA_V) -D o_SNOVA=$(SNOVA_O) -D l_SNOVA=$(SNOVA_L) -D sk_is_seed=$(SK_IS_SEED) -D CRYPTO_ALGNAME=\"$(CRYPTO_ALGNAME)\" -D PK_EXPAND_SHAKE=$(PK_EXPAND_SHAKE)
SNOVA_PARAMS += -D OPTIMISATION=$(OPTIMISATION)

TEST_SPEED_N ?= 2048

info:
	@echo "SNOVA MAKE ..."
	@echo "CRYPTO_ALGNAME: $(CRYPTO_ALGNAME)"
	@echo "PK_EXPAND:      $(PK_EXPAND)"
	@echo "OPTIMISATION:   $(OPTIMISATION)"
	@echo "==================================="

clean:
	rm -f ./build/*.o ./*/*.o ./build/*.a *.a 

clean_all: 
	rm -f ./build/*.o ./*/*.o ./*/*.a *.a *.log *.req *.rsp

build/libsnovasym.a: $(SYMMETRICLIBO)
	ar rcs build/libsnovasym.a $(SYMMETRICLIBO)

build/rng.o: 
	$(CC) $(CFLAGS) -c -o ./build/rng.o ./rng.c

build/snova.o: build/rng.o
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c -o ./build/snova.o ./snova.c

build/sign.o: build/snova.o
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c -o ./build/sign.o ./sign.c

test: info build/rng.o build/snova.o $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) $(OLIST) ./test/test.c -o test.a ${LIBS}
	./test.a > test_$(CRYPTO_ALGNAME).log

test_api: info build/rng.o build/snova.o build/sign.o $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) $(OLIST) ./test/test_api.c ./build/sign.o -o test_api.a ${LIBS}
	./test_api.a > test_api_$(CRYPTO_ALGNAME).log

test_speed: info build/rng.o build/snova.o $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -D TEST_SPEED_N=$(TEST_SPEED_N) $(OLIST) ./test_speed/test_speed.c -o test_speed.a ${LIBS}
	./test_speed.a

PQCgenKAT: info build/sign.o $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) $(OLIST) ./build/sign.o ./PQCgenKAT_sign.c -o ./PQCgenKAT.a ${LIBS}
	./PQCgenKAT.a

