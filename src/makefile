-include snova_config.mk

ARCH=native

# Autodetect AVX
AVX2 := $(findstring AVX2, $(shell gcc -march=$(ARCH) -dM -E - < /dev/null))
AVX512 := $(findstring AVX512, $(shell gcc -march=$(ARCH) -dM -E - < /dev/null))

# Autodetect ARM
ARM_FEATURE_AES := $(findstring ARM_FEATURE_AES, $(shell gcc -march=$(ARCH) -dM -E - < /dev/null))
ARM_NEON := $(findstring ARM_NEON, $(shell gcc -march=$(ARCH) -dM -E - < /dev/null))

# snova params
SNOVA_V ?= 24
SNOVA_O ?= 5
SNOVA_L ?= 4

# 0: sk = esk; 1: sk = ssk
SK_IS_SEED ?= 0
ifeq ($(SK_IS_SEED), 1)
    SK_TYPE = SSK
else
    SK_TYPE = ESK
endif

# 0: disable; 1: enable
PK_EXPAND_SHAKE ?= 0

# 0: Reference, 1:Optimised, 2:Vectorized AVX2/NEON
ifeq ($(AVX2), AVX2)
OPTIMISATION ?= 2
else
ifeq ($(ARM_NEON), ARM_NEON)
OPTIMISATION ?= 2
else
OPTIMISATION ?= 1
endif
endif

ifeq ($(PK_EXPAND_SHAKE), 1)
    PK_EXPAND = _SHAKE
else
    PK_EXPAND = 
endif

CC = gcc
AR = ar

ifeq ($(ARM_FEATURE_AES), ARM_FEATURE_AES)
CFLAGS = -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -ftree-vectorize -O3 -march=armv8-a+crypto -mtune=$(ARCH)
else
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -ftree-vectorize -O3 -march=$(ARCH) -mtune=$(ARCH)
endif
CFLAGS += -I. -I./snova_plasma -I./snova_plasma/plasma_general
# CFLAGS += -g -fstack-usage

BUILD_OUT_PATH = build/
BUILD_LIBO_OUT_PATH = build_libo

OLIST = $(BUILD_OUT_PATH)snova.o $(BUILD_OUT_PATH)snova_common.o $(BUILD_OUT_PATH)snova_aes.o
ifeq ($(OPTIMISATION), 0)
OLIST += $(BUILD_OUT_PATH)snova_shake_ref.o 
else
OLIST += $(BUILD_OUT_PATH)snova_shake.o 
endif

SYMMETRICLIB = shake/KeccakHash.c shake/KeccakP-1600-opt64.c shake/KeccakSponge.c
SYMMETRICLIB += aes/aes_c.c aes/oqs_secure_free.c
ifeq ($(AVX512), AVX512)
SYMMETRICLIB += shake/KeccakP-1600-times4-SIMD256.c shake/KeccakP-1600-times8-SIMD512.c
SYMMETRICLIB += aes/aes128_ni.c aes/aes256_ni.c
else
ifeq ($(AVX2), AVX2)
SYMMETRICLIB += shake/KeccakP-1600-times4-SIMD256.c
SYMMETRICLIB += aes/aes128_ni.c aes/aes256_ni.c
else
ifeq ($(ARM_FEATURE_AES), ARM_FEATURE_AES)
SYMMETRICLIB += aes/aes128_armv8.c aes/aes256_armv8.c
endif
endif
endif
# SYMMETRICLIBO = $(SYMMETRICLIB:.c=.o)
SYMMETRICLIBO = $(patsubst %.c, $(BUILD_LIBO_OUT_PATH)/%.o, $(notdir $(SYMMETRICLIB)))
STATICLIB = $(BUILD_LIBO_OUT_PATH)/libsnovasym.a
LIBS = -L$(BUILD_LIBO_OUT_PATH) -lsnovasym

CRYPTO_ALGNAME = SNOVA_$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L)$(PK_EXPAND)
SNOVA_PARAMS = -D v_SNOVA=$(SNOVA_V) -D o_SNOVA=$(SNOVA_O) -D l_SNOVA=$(SNOVA_L) -D sk_is_seed=$(SK_IS_SEED) -D PK_EXPAND_SHAKE=$(PK_EXPAND_SHAKE)
SNOVA_PARAMS += -D OPTIMISATION=$(OPTIMISATION)

TEST_SPEED_N ?= 2048
WASM = 0

all: clean_without_libsnovasym test_speed

init:
	@mkdir -p $(BUILD_OUT_PATH)

init_libo:
	@mkdir -p $(BUILD_LIBO_OUT_PATH)

info:
	@echo "==================================="
	@echo "SNOVA MAKE ..."
	@echo "CRYPTO_ALGNAME: $(CRYPTO_ALGNAME)"
	@echo "SK_IS_SEED:     $(SK_IS_SEED)"
	@echo "PK_EXPAND:      $(PK_EXPAND)"
	@echo "OPTIMISATION:   $(OPTIMISATION)"
	@echo "AVX2:           $(AVX2)"
	@echo "AVX512:         $(AVX512)"
	@echo "==================================="

wasm:
	$(eval CC := emcc)
	$(eval AR := emar)
	$(eval CFLAGS := $(filter-out -march%, $(CFLAGS)))
	$(eval CFLAGS := $(filter-out -mtune%, $(CFLAGS)))
	$(eval CFLAGS += -msimd128)
	$(eval OPTIMISATION := 1)
	$(eval WASM := 1)
	$(eval WASM_OUT_PATH := ./)
	$(eval WASM_OUT_EXT := wasm)
	@echo "==================================="
	@echo "SETTING WASM ..."
	@echo "WASM:           $(WASM)"
	@echo "WASM_OUT_PATH:  $(WASM_OUT_PATH)"
	@echo "WASM_OUT_EXT:   $(WASM_OUT_EXT)"
	@echo "==================================="

clean: clean_kat
	rm -f ./*/*.o ./*/*.a *.a *.wasm *.js
	rm -rf ./build/ 2>/dev/null
	rm -rf ./build_libo/* 2>/dev/null

clean_all: clean
	rm -f *.log

clean_without_libsnovasym:
	rm -f *.wasm *.js
	rm -rf ./build/* 2>/dev/null

clean_kat: 
	rm -f *.req *.rsp

$(BUILD_LIBO_OUT_PATH)/%.o: shake/%.c | init_libo
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_LIBO_OUT_PATH)/%.o: aes/%.c | init_libo
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_LIBO_OUT_PATH)/libsnovasym.a: $(SYMMETRICLIBO)
	$(AR) rcs $@ $^

build/rng.o:
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c nistkat/rng.c -o build/rng.o

build/snova_common.o:
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c snova_common.c -o build/snova_common.o

build/snova_shake_ref.o:
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c snova_shake_ref.c -o build/snova_shake_ref.o

build/snova_shake.o:
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c snova_plasma/snova_shake.c -o build/snova_shake.o

build/snova_aes.o:
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c snova_aes.c -o build/snova_aes.o

build/snova.o:
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c snova.c -o build/snova.o

build/sign.o: build/rng.o build/snova.o 
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -c sign.c -o build/sign.o

test: init info build/rng.o $(OLIST) $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) $(OLIST) build/rng.o test/test.c -o test.a $(LIBS)
	./test.a > test_$(CRYPTO_ALGNAME)_OPT_$(OPTIMISATION).log

test_api: init info build/sign.o $(OLIST) $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) $(OLIST) build/rng.o test/test_api.c build/sign.o -o test_api.a $(LIBS)
	./test_api.a > test_api_$(CRYPTO_ALGNAME)_$(SK_TYPE)_OPT_$(OPTIMISATION).log

test_speed: init info build/rng.o $(OLIST) $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -D TEST_SPEED_N=$(TEST_SPEED_N) -D CRYPTO_ALGNAME=\"$(CRYPTO_ALGNAME)\" $(OLIST) build/rng.o test_speed/test_speed.c -o test_speed.a $(LIBS)
	./test_speed.a

PQCgenKAT: init info build/sign.o $(OLIST) $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -D CRYPTO_ALGNAME=\"$(CRYPTO_ALGNAME)\" $(OLIST) build/rng.o build/sign.o nistkat/PQCgenKAT_sign.c -o PQCgenKAT.a $(LIBS)
	./PQCgenKAT.a

build_wasmapi: init wasm info $(OLIST) $(STATICLIB)
	$(CC) $(CFLAGS) $(SNOVA_PARAMS) -D CRYPTO_ALGNAME=\"$(CRYPTO_ALGNAME)\" $(OLIST) wasmapi.c -o $(WASM_OUT_PATH)$(CRYPTO_ALGNAME).$(WASM_OUT_EXT) -s WASM=1 -s EXPORTED_FUNCTIONS='["_malloc", "_free"]' -s SAFE_HEAP=1 -s TOTAL_STACK=17825792 -s ALLOW_MEMORY_GROWTH=1 $(LIBS) --no-entry -s INITIAL_MEMORY=20971520
