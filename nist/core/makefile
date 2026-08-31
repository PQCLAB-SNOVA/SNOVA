# snova - SNOVA Round 3 reference implementation
#
# Build flow:
# ALWAYS pass the complete 7-tuple v/o/q/l/r/m1/m2 - not "add SNOVA_R only when rectangular".
#   3 of the 9 default sets OVERRIDE m1 (29_3/27_4/38_4); omitting it silently builds a DIFFERENT
#   scheme. m2 cannot BE overridden; passing it only feeds a _Static_assert cross-check.
#   make kat-verify    SNOVA_V=27 SNOVA_O=5 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=4 SNOVA_M1=5 SNOVA_M2=80
#   make bench-details SNOVA_V=29 SNOVA_O=3 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=8 SNOVA_M1=5 SNOVA_M2=96 BENCH_N=512
#   AES=1 -> AES XOF, else SHAKE (AES is the knob; PK_EXPAND_SHAKE is derived and cannot be set).
#
# Default params: 27_5_16_4_4_5_80 (full ordered name), SHAKE. sk is always seed-form.

CC      ?= gcc
ARCH    ?= x86_avx2
CFLAGS  ?= -std=c11 -O3 -Wall -ftree-vectorize
-include snova_config.mk
SNOVA_V ?= 27
SNOVA_O ?= 5
SNOVA_L ?= 4
SNOVA_Q ?= 16
SNOVA_R ?= $(SNOVA_L)
PK_EXPAND_SHAKE ?= 1
KAT_N    ?= 10
BENCH_N  ?= 512
KAT_RSP  ?= ref_kat/$(KAT_RSP_NAME)
SNOVA_NAME ?= SNOVA

AES    ?= 0
ifeq ($(AES),1)
override PK_EXPAND_SHAKE := 0
R3_XOF_SUFFIX := _AES
else
override PK_EXPAND_SHAKE := 1
R3_XOF_SUFFIX :=
endif

GFNI ?= 1
ifeq ($(GFNI),0)
NOGFNI_KERNEL_DEFAULT := 1
else
NOGFNI_KERNEL_DEFAULT := 0
endif
VNNI ?= $(GFNI)
RCT_PSHUFB ?= 0
JOG ?= 0
SNOVA_M1 ?=
RCT_QRP16    ?= $(NOGFNI_KERNEL_DEFAULT)
RCT_MULLO    ?= $(NOGFNI_KERNEL_DEFAULT)
ifeq ($(ARCH), x86_avx2)
ifeq ($(GFNI),1)
ARCH_FLAGS := -DSNOVA_ARCH_X86_AVX2=1 -mavx2 -maes -mpclmul -mgfni
else
ARCH_FLAGS := -DSNOVA_ARCH_X86_AVX2=1 -mavx2 -maes -mpclmul
endif
ifeq ($(VNNI),1)
ARCH_FLAGS += -mavxvnni
endif
PLATFORM_SRC := \
    src/platforms/x86_avx2/xgf16_avx.c \
    src/platforms/x86_avx2/aes_ni.c
else ifeq ($(ARCH), portable_opt)
ARCH_FLAGS := -DSNOVA_ARCH_PORTABLE_OPT=1
PLATFORM_SRC := \
    src/platforms/portable_opt/xgf16_opt.c \
    src/platforms/ref/aes_ref.c
else ifeq ($(ARCH), ref)
ARCH_FLAGS := -DSNOVA_ARCH_REF=1
PLATFORM_SRC := \
    src/platforms/ref/xgf16_ref.c \
    src/platforms/ref/aes_ref.c
else
$(error Unknown ARCH=$(ARCH). Use ARCH=x86_avx2 (default), portable_opt, or ref)
endif

SNOVA_FLAGS := \
    -DSNOVA_V=$(SNOVA_V) -DSNOVA_O=$(SNOVA_O) -DSNOVA_L=$(SNOVA_L) \
    -DSNOVA_Q=$(SNOVA_Q) -DSNOVA_R=$(SNOVA_R) \
    -DSNOVA_PK_EXPAND_SHAKE=$(PK_EXPAND_SHAKE)
ifeq ($(RCT_PSHUFB),1)
SNOVA_FLAGS += -DRCT_FORCE_PSHUFB
endif
ifeq ($(JOG),1)
SNOVA_FLAGS += -DRCT_FORCE_JOG=1
endif
ifneq ($(SNOVA_M1),)
SNOVA_FLAGS += -DSNOVA_M1=$(SNOVA_M1)
endif
ifneq ($(SNOVA_M2),)
SNOVA_FLAGS += -DSNOVA_M2_ASSERT=$(SNOVA_M2)
endif
SNOVA_PARAM_FLAGS := -DSNOVA_V=$(SNOVA_V) -DSNOVA_O=$(SNOVA_O) -DSNOVA_Q=$(SNOVA_Q) \
    -DSNOVA_L=$(SNOVA_L) -DSNOVA_R=$(SNOVA_R) \
    $(if $(SNOVA_M1),-DSNOVA_M1=$(SNOVA_M1),) $(if $(SNOVA_M2),-DSNOVA_M2_ASSERT=$(SNOVA_M2),)
VERIFY_STREAM ?= 0
ifeq ($(VERIFY_STREAM),1)
SNOVA_FLAGS += -DSNOVA_VERIFY_STREAM=1
endif
SIGN_STREAM ?= 0
ifeq ($(SIGN_STREAM),1)
SNOVA_FLAGS += -DSNOVA_SIGN_STREAM=1
endif
KEYGEN_STREAM ?= 0
ifeq ($(KEYGEN_STREAM),1)
SNOVA_FLAGS += -DSNOVA_KEYGEN_STREAM=1
endif
PKX_PGEN ?= 0
ifeq ($(PKX_PGEN),1)
SNOVA_FLAGS += -DSNOVA_PKX_PGEN=1
endif
RCT_TILE4 ?=
ifeq ($(RCT_TILE4),1)
SNOVA_FLAGS += -DRCT_TILE4=1
endif
ifeq ($(RCT_TILE4),0)
SNOVA_FLAGS += -DRCT_TILE4=0
endif
RCT_TILE4_PREBAKE ?= 0
ifeq ($(RCT_TILE4_PREBAKE),1)
SNOVA_FLAGS += -DRCT_TILE4_PREBAKE=1
endif
RCT_PROFILE ?= 0
ifeq ($(RCT_PROFILE),1)
SNOVA_FLAGS += -DRCT_PROFILE
$(warning RCT_PROFILE=1: bench numbers include instrumentation overhead - not a baseline)
endif
RCT_FOLD_SEG ?= 1
ifeq ($(RCT_FOLD_SEG),0)
SNOVA_FLAGS += -DRCT_FOLD_SEG=0
endif
RCT_FOLD_SEG_M4 ?= 1
ifeq ($(RCT_FOLD_SEG_M4),0)
SNOVA_FLAGS += -DRCT_FOLD_SEG_M4=0
endif
RCT_SIGN_P2_RSWAP ?= 1
ifeq ($(RCT_SIGN_P2_RSWAP),0)
SNOVA_FLAGS += -DRCT_SIGN_P2_RSWAP=0
endif
RCT_SIGN_P2_B1LANE ?= 1
ifeq ($(RCT_SIGN_P2_B1LANE),0)
SNOVA_FLAGS += -DRCT_SIGN_P2_B1LANE=0
endif
WRAPPER_STACK ?= 1
SNOVA_FLAGS += -DSNOVA_WRAPPER_STACK=$(WRAPPER_STACK)
ifeq ($(RCT_QRP16),1)
SNOVA_FLAGS += -DRCT_QRP16=1
endif
ifeq ($(RCT_MULLO),1)
SNOVA_FLAGS += -DRCT_MULLO=1
endif
ifeq ($(SNOVA_Q),16)
ifneq ($(SNOVA_R),$(SNOVA_L))
RCT_COEF_MULLO ?= 1
endif
endif
RCT_COEF_MULLO ?= $(NOGFNI_KERNEL_DEFAULT)
ifeq ($(RCT_COEF_MULLO),1)
SNOVA_FLAGS += -DRCT_COEF_MULLO=1
endif
RCT_PKX_JOGSIMD ?= 1
ifeq ($(RCT_PKX_JOGSIMD),0)
SNOVA_FLAGS += -DRCT_PKX_JOGSIMD=0
endif
SKX_SLIM ?= 0
ifeq ($(SKX_SLIM),1)
SNOVA_FLAGS += -DRCT_SKX_SLIM=1
endif
ifdef RCT_FOLD_WIDE
SNOVA_FLAGS += -DRCT_FOLD_WIDE=$(RCT_FOLD_WIDE)
endif
CM_SQ_L4 := 0
ifeq ($(SNOVA_Q),16)
ifeq ($(filter-out $(SNOVA_L),$(SNOVA_R)),)
ifeq ($(SNOVA_L),4)
CM_SQ_L4 := 1
endif
endif
endif
SQ_CM_DEFAULT := $(if $(filter 1,$(NOGFNI_KERNEL_DEFAULT) $(CM_SQ_L4)),1,0)
RCT_SQ_CM_LEFT ?= $(SQ_CM_DEFAULT)
ifeq ($(RCT_SQ_CM_LEFT),1)
SNOVA_FLAGS += -DRCT_SQ_CM_LEFT=1
endif
RCT_SQ_CM_TIGHT ?= $(SQ_CM_DEFAULT)
ifeq ($(RCT_SQ_CM_TIGHT),1)
SNOVA_FLAGS += -DRCT_SQ_CM_TIGHT=1
endif
RCT_SQ_CM_S2 ?= $(NOGFNI_KERNEL_DEFAULT)
ifeq ($(RCT_SQ_CM_S2),1)
SNOVA_FLAGS += -DRCT_SQ_CM_S2=1
endif
RCT_SQ_RIGHT_DUAL ?= $(NOGFNI_KERNEL_DEFAULT)
ifeq ($(RCT_SQ_RIGHT_DUAL),1)
SNOVA_FLAGS += -DRCT_SQ_RIGHT_DUAL=1
endif
RCT_SQ_CM_S3 ?= $(SQ_CM_DEFAULT)
ifeq ($(RCT_SQ_CM_S3),1)
SNOVA_FLAGS += -DRCT_SQ_CM_S3=1
endif
RCT_T12_MULLO ?= $(NOGFNI_KERNEL_DEFAULT)
ifeq ($(RCT_T12_MULLO),1)
SNOVA_FLAGS += -DRCT_T12_MULLO=1
endif
RCT_KG_MULLO ?= $(NOGFNI_KERNEL_DEFAULT)
ifeq ($(RCT_KG_MULLO),1)
SNOVA_FLAGS += -DRCT_KG_MULLO=1
endif
RCT_OQ_WHIPVEC ?= 1
ifeq ($(RCT_OQ_WHIPVEC),1)
SNOVA_FLAGS += -DRCT_OQ_WHIPVEC=1
endif
RCT_OQ_DEVFLOW ?= 1
ifeq ($(RCT_OQ_DEVFLOW),1)
SNOVA_FLAGS += -DRCT_OQ_DEVFLOW=1
endif

INCLUDES := -Isrc -Isrc/third_party/shake_xkcp

SNOVA_CORE_SRC := \
    src/snova_core/snova.c \
    src/snova_core/drbg.c \
    src/primitives/snova_aes.c

GF16_CORE_SRC := \
    src/gf16_core/gf16.c \
    src/gf16_core/gf16m.c \
    src/gf16_core/gf16m_mul.c

ifeq ($(PK_EXPAND_SHAKE), 1)
SHAKE_SRC := \
    src/third_party/shake_xkcp/KeccakHash.c \
    src/third_party/shake_xkcp/SimpleFIPS202.c \
    src/third_party/shake_xkcp/KeccakP-1600-opt64.c \
    src/third_party/shake_xkcp/KeccakSponge.c
ifeq ($(ARCH), x86_avx2)
SHAKE_SRC += \
    src/third_party/shake_xkcp/snova_shake_opt.c \
    src/third_party/shake_xkcp/KeccakP-1600-times4-SIMD256.c
else
SHAKE_SRC += src/third_party/shake_xkcp/snova_shake_ref.c
endif
else
SHAKE_SRC :=
endif

LEGACY_KECCAK_SRC := src/third_party/keccak_legacy/keccak.c

BUILD_DIR := build/$(ARCH)_$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)_$(SNOVA_R)_$(PK_EXPAND_SHAKE)$(if $(SNOVA_M1),_m$(SNOVA_M1),)

ALL_SRC := $(SNOVA_CORE_SRC) $(GF16_CORE_SRC) $(PLATFORM_SRC) $(SHAKE_SRC) $(LEGACY_KECCAK_SRC)
COMPILE_FLAGS := $(CFLAGS) $(ARCH_FLAGS) $(SNOVA_FLAGS) $(INCLUDES)

ALL_HDR := $(wildcard src/*.h src/*/*.h src/*/*/*.h)
FLAGSTAMP := $(BUILD_DIR)/.flags
$(FLAGSTAMP): FORCE | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo 'CC=$(CC) COMPILE_FLAGS=$(COMPILE_FLAGS)' | cmp -s - $@ || \
	 echo 'CC=$(CC) COMPILE_FLAGS=$(COMPILE_FLAGS)' > $@
FORCE:
.PHONY: FORCE
$(BUILD_DIR)/kat_diff $(BUILD_DIR)/bench_cycles $(BUILD_DIR)/bench_details \
$(BUILD_DIR)/memtest \
$(BUILD_DIR)/PQCgenKAT_sign \
$(BUILD_DIR)/fuzz_verify $(BUILD_DIR)/kg_stream_test $(BUILD_DIR)/pkx_pgen_test \
$(BUILD_DIR)/sign_stream_test $(BUILD_DIR)/stream_test \
$(BUILD_DIR)/wasm_selftest: \
	$(ALL_HDR) $(FLAGSTAMP)


.PHONY: all kat-verify bench bench-details memtest static-size PQCgenKAT wasm build_wasmapi clean_all clean_without_libsnovasym fuzz clean help

all: kat-verify


kat-verify: $(BUILD_DIR)/PQCgenKAT_sign
	@echo "[kat-verify] $(KAT_RSP_NAME): regenerate .rsp + digest vs ref_kat/KAT_DIGESTS.sha256"
	@cd $(BUILD_DIR) && ./PQCgenKAT_sign >/dev/null
	@exp=$$(awk -v f="$(KAT_RSP_NAME)" '$$2==f {print $$1}' ref_kat/KAT_DIGESTS.sha256); \
	 if [ -z "$$exp" ]; then \
	   echo "== KAT-VERIFY FAIL: no digest for $(KAT_RSP_NAME) in ref_kat/KAT_DIGESTS.sha256 =="; \
	   exit 1; \
	 fi; \
	 act=$$(sha256sum $(BUILD_DIR)/$(KAT_RSP_NAME) | awk '{print $$1}'); \
	 if [ "$$exp" = "$$act" ]; then \
	   echo "== KAT-VERIFY $(KAT_RSP_NAME) DIGEST MATCH =="; \
	 else \
	   echo "== KAT-VERIFY $(KAT_RSP_NAME) DIGEST MISMATCH =="; \
	   echo "   expected $$exp"; \
	   echo "   actual   $$act"; \
	   exit 1; \
	 fi

# Single-param cycles bench (E2E full chain)
bench: $(BUILD_DIR)/bench_cycles
	@$(BUILD_DIR)/bench_cycles $(BENCH_N)

bench-details: $(BUILD_DIR)/bench_details
	@$(BUILD_DIR)/bench_details $(BENCH_N)

# Single-param NIST KAT generator (produces PQCsignKAT_*.req/.rsp)
PQCgenKAT: $(BUILD_DIR)/PQCgenKAT_sign
	@echo "Generating PQCgenKAT for $(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L) (shake=$(PK_EXPAND_SHAKE))"
	@cd $(BUILD_DIR) && ./PQCgenKAT_sign

$(BUILD_DIR)/kat_diff: $(ALL_SRC) test/kat_diff.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L)/shake=$(PK_EXPAND_SHAKE)] build kat_diff"
	@$(CC) $(COMPILE_FLAGS) $(ALL_SRC) test/kat_diff.c -o $@

$(BUILD_DIR)/bench_cycles: $(ALL_SRC) test/bench_cycles.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L)/shake=$(PK_EXPAND_SHAKE)] build bench_cycles"
	@$(CC) $(COMPILE_FLAGS) $(ALL_SRC) test/bench_cycles.c -o $@ -DCRYPTO_ALGNAME='"$(SNOVA_NAME)"'

$(BUILD_DIR)/bench_details: $(ALL_SRC) test/bench_details.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L)/shake=$(PK_EXPAND_SHAKE)] build bench_details"
	@$(CC) $(COMPILE_FLAGS) $(ALL_SRC) test/bench_details.c -o $@

MEM_API  ?= sign
MEM_REPS ?= 1
$(BUILD_DIR)/memtest: $(ALL_SRC) test/memtest.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L)/shake=$(PK_EXPAND_SHAKE)] build memtest"
	@$(CC) $(COMPILE_FLAGS) $(ALL_SRC) test/memtest.c -o $@
memtest: $(BUILD_DIR)/memtest
	@$(BUILD_DIR)/memtest $(MEM_API) $(MEM_REPS)
static-size: $(BUILD_DIR)/memtest
	@size $(BUILD_DIR)/memtest

R3_M1 := $(if $(SNOVA_M1),$(SNOVA_M1),$(shell expr $(SNOVA_O) '*' $(SNOVA_R) / $(SNOVA_L)))
R3_M2 := $(shell expr $(SNOVA_O) '*' $(SNOVA_L) '*' $(SNOVA_R))
PQCgenKAT_ALG := SNOVA_$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)_$(SNOVA_R)_$(R3_M1)_$(R3_M2)$(R3_XOF_SUFFIX)
KAT_RSP_NAME := PQCsignKAT_$(PQCgenKAT_ALG).rsp

$(BUILD_DIR)/PQCgenKAT_sign: $(ALL_SRC) src/nistkat/PQCgenKAT_sign.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L)/shake=$(PK_EXPAND_SHAKE)] build PQCgenKAT_sign ($(PQCgenKAT_ALG))"
	@$(CC) $(COMPILE_FLAGS) -DCRYPTO_ALGNAME='"$(PQCgenKAT_ALG)"' \
	    $(ALL_SRC) src/nistkat/PQCgenKAT_sign.c -o $@

nistkat: $(BUILD_DIR)/distPQCgenKAT_sign
	cd $(BUILD_DIR) && ./distPQCgenKAT_sign

$(BUILD_DIR)/distPQCgenKAT_sign: $(ALL_SRC) src/nistkat/PQCgenKAT_sign.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_L)/shake=$(PK_EXPAND_SHAKE)] build PQCgenKAT_sign"
	@$(CC) $(COMPILE_FLAGS) -DCRYPTO_ALGNAME='"$(SNOVA_NAME)"' \
	    $(ALL_SRC) src/nistkat/PQCgenKAT_sign.c -o $@

$(BUILD_DIR):
	@mkdir -p $@

FUZZ_SAN     ?= -fsanitize=address,undefined -fno-omit-frame-pointer -g
FUZZ_ARGS    ?= 200000
FUZZ_CLANG   := $(shell command -v clang 2>/dev/null)
FUZZ_LF_OK   := $(shell printf 'int LLVMFuzzerTestOneInput(const unsigned char*d,unsigned long s){(void)d;(void)s;return 0;}\n' > /tmp/.snova_lfprobe.c 2>/dev/null && clang -fsanitize=fuzzer,address,undefined /tmp/.snova_lfprobe.c -o /tmp/.snova_lfprobe.out >/dev/null 2>&1 && echo 1 || echo 0)

$(BUILD_DIR)/fuzz_verify: $(ALL_SRC) test/fuzz_verify.c | $(BUILD_DIR)
	@mkdir -p $(@D)
ifeq ($(FUZZ_LF_OK),1)
	@echo "[fuzz] clang libFuzzer AVAILABLE -> coverage-guided build (-fsanitize=fuzzer,address,undefined)"
	@clang $(COMPILE_FLAGS) -DFUZZ_LIBFUZZER=1 -fsanitize=fuzzer,address,undefined -g \
	    $(ALL_SRC) test/fuzz_verify.c -o $@ $(LIBS)
else
	@echo "[fuzz] clang/libFuzzer UNAVAILABLE -> gcc ASan+UBSan deterministic corpus-replay fallback"
	@$(CC) $(COMPILE_FLAGS) $(FUZZ_SAN) $(ALL_SRC) test/fuzz_verify.c -o $@ $(LIBS)
endif

fuzz: $(BUILD_DIR)/fuzz_verify
ifeq ($(FUZZ_LF_OK),1)
	@echo "[fuzz/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)_$(SNOVA_R)] libFuzzer campaign"
	@ASAN_OPTIONS=abort_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
	    $(BUILD_DIR)/fuzz_verify $(FUZZ_LIBFUZZER_ARGS)
else
	@echo "[fuzz/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)_$(SNOVA_R)] fallback ASan+UBSan campaign"
	@UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 $(BUILD_DIR)/fuzz_verify $(FUZZ_ARGS)
endif

.PHONY: stream-test
$(BUILD_DIR)/stream_test: $(ALL_SRC) test/stream_verify_test.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)] build stream_test (VERIFY_STREAM=$(VERIFY_STREAM))"
	@$(CC) $(COMPILE_FLAGS) $(ALL_SRC) test/stream_verify_test.c -o $@
stream-test: $(BUILD_DIR)/stream_test
	@$(BUILD_DIR)/stream_test

.PHONY: sign-stream-test
$(BUILD_DIR)/sign_stream_test: $(ALL_SRC) test/sign_stream_test.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)_$(SNOVA_R)] build sign_stream_test (SIGN_STREAM=$(SIGN_STREAM))"
	@$(CC) $(COMPILE_FLAGS) $(ALL_SRC) test/sign_stream_test.c -o $@
sign-stream-test: $(BUILD_DIR)/sign_stream_test
	@$(BUILD_DIR)/sign_stream_test

.PHONY: kg-stream-test pkx-pgen-test
$(BUILD_DIR)/kg_stream_test: $(ALL_SRC) test/kg_stream_test.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)_$(SNOVA_R)] build kg_stream_test (KEYGEN_STREAM=$(KEYGEN_STREAM))"
	@$(CC) $(COMPILE_FLAGS) $(ALL_SRC) test/kg_stream_test.c -o $@
kg-stream-test: $(BUILD_DIR)/kg_stream_test
	@$(BUILD_DIR)/kg_stream_test

$(BUILD_DIR)/pkx_pgen_test: $(ALL_SRC) test/pkx_pgen_test.c | $(BUILD_DIR)
	@mkdir -p $(@D)
	@echo "[snova/$(ARCH)/$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)_$(SNOVA_R)] build pkx_pgen_test (PKX_PGEN=$(PKX_PGEN))"
	@$(CC) $(COMPILE_FLAGS) $(ALL_SRC) test/pkx_pgen_test.c -o $@
pkx-pgen-test: $(BUILD_DIR)/pkx_pgen_test
	@$(BUILD_DIR)/pkx_pgen_test

WASM_ARCH ?= wasm_simd128
ifeq ($(WASM_ARCH), wasm_simd128)
WASM_ARCH_FLAGS := -DSNOVA_ARCH_WASM_SIMD128=1
WASM_PLATFORM := \
    src/platforms/wasm_simd128/xgf16_wasm.c \
    src/platforms/ref/aes_ref.c
else ifeq ($(WASM_ARCH), portable_opt)
WASM_ARCH_FLAGS := -DSNOVA_ARCH_PORTABLE_OPT=1
WASM_PLATFORM := \
    src/platforms/portable_opt/xgf16_opt.c \
    src/platforms/ref/aes_ref.c
else
WASM_ARCH_FLAGS := -DSNOVA_ARCH_REF=1
WASM_PLATFORM := \
    src/platforms/ref/xgf16_ref.c \
    src/platforms/ref/aes_ref.c
endif
WASM_BUILD := build/wasm_$(SNOVA_V)_$(SNOVA_O)_$(SNOVA_Q)_$(SNOVA_L)_$(SNOVA_R)_$(PK_EXPAND_SHAKE)$(if $(SNOVA_M1),_m$(SNOVA_M1),)
WASM_ALG := $(PQCgenKAT_ALG)$(if $(filter 1,$(PK_EXPAND_SHAKE)),_SHAKE,)
WASM_OUT := $(WASM_BUILD)/$(WASM_ALG).wasm

ifeq ($(PK_EXPAND_SHAKE), 1)
WASM_SHAKE_SRC := \
    src/third_party/shake_xkcp/KeccakHash.c \
    src/third_party/shake_xkcp/SimpleFIPS202.c \
    src/third_party/shake_xkcp/KeccakP-1600-opt64.c \
    src/third_party/shake_xkcp/KeccakSponge.c \
    src/third_party/shake_xkcp/snova_shake_ref.c
else
WASM_SHAKE_SRC :=
endif

WASM_SRC := \
    $(SNOVA_CORE_SRC) \
    $(GF16_CORE_SRC) \
    $(WASM_PLATFORM) \
    $(WASM_SHAKE_SRC) \
    $(LEGACY_KECCAK_SRC) \
    src/wasmapi.c

EMCC := $(shell command -v emcc 2>/dev/null)

$(WASM_BUILD):
	@mkdir -p $@

wasm: | $(WASM_BUILD)
ifeq ($(EMCC),)
	@echo "[snova] emcc not found; falling back to gcc compile-test of wasmapi.c"
	@$(CC) -O2 -std=c11 $(WASM_ARCH_FLAGS) -Isrc -Isrc/third_party/shake_xkcp \
	    $(SNOVA_PARAM_FLAGS) -DSNOVA_PK_EXPAND_SHAKE=$(PK_EXPAND_SHAKE) \
	    -c src/wasmapi.c -o $(WASM_BUILD)/wasmapi.o
	@echo "  $(WASM_BUILD)/wasmapi.o ok (compile-test only - install emscripten for real wasm)"
else
	@echo "[snova] emcc found at $(EMCC); building wasm (WASM_ARCH=$(WASM_ARCH), $(WASM_ALG))"
	@$(EMCC) -O3 -msimd128 $(WASM_ARCH_FLAGS) -Isrc -Isrc/third_party/shake_xkcp \
	    $(SNOVA_PARAM_FLAGS) -DSNOVA_PK_EXPAND_SHAKE=$(PK_EXPAND_SHAKE) \
	    -DCRYPTO_ALGNAME='"$(WASM_ALG)"' \
	    -s WASM=1 \
	    -s EXPORTED_FUNCTIONS='["_malloc","_free"]' \
	    -s ALLOW_MEMORY_GROWTH=1 \
	    -s TOTAL_STACK=16777216 \
	    --no-entry \
	    $(WASM_SRC) \
	    -o $(WASM_OUT)
	@echo "  Wasm output: $(WASM_OUT) (size: $$(stat -c%s $(WASM_OUT)) bytes)"
endif

WASM_OUT_PATH ?= $(WASM_BUILD)/
WASM_OUT_EXT  ?= wasm
WASM_APINAME  := $(PQCgenKAT_ALG)
build_wasmapi: | $(WASM_BUILD)
	@if [ -z "$(EMCC)" ]; then echo "[snova] build_wasmapi needs emscripten (emcc)"; exit 1; fi
	@mkdir -p $(WASM_OUT_PATH)
	@echo "[snova] build_wasmapi $(WASM_APINAME) (WASM_ARCH=$(WASM_ARCH)) -> $(WASM_OUT_PATH)"
	@$(EMCC) -O3 -msimd128 $(WASM_ARCH_FLAGS) -Isrc -Isrc/third_party/shake_xkcp \
	    $(SNOVA_PARAM_FLAGS) -DSNOVA_PK_EXPAND_SHAKE=$(PK_EXPAND_SHAKE) \
	    -DCRYPTO_ALGNAME='"$(WASM_APINAME)"' \
	    -s WASM=1 -s EXPORTED_FUNCTIONS='["_malloc","_free"]' \
	    -s ALLOW_MEMORY_GROWTH=1 -s TOTAL_STACK=16777216 --no-entry \
	    $(WASM_SRC) -o $(WASM_OUT_PATH)$(WASM_APINAME).$(WASM_OUT_EXT)
	@echo "  -> $(WASM_OUT_PATH)$(WASM_APINAME).$(WASM_OUT_EXT) ($$(stat -c%s $(WASM_OUT_PATH)$(WASM_APINAME).$(WASM_OUT_EXT)) bytes)"

clean_all: clean
clean_without_libsnovasym: clean

build:
	@mkdir -p $@

clean:
	@rm -rf build/
	@echo "[snova] cleaned build/"

help:
	@echo "snova - SNOVA Round 3 reference implementation"
	@echo ""
	@echo "Single-param targets:"
	@echo "  Parameters are ALWAYS the complete 7-tuple (see snova_config.src):"
	@echo "    SNOVA_V= SNOVA_O= SNOVA_Q= SNOVA_L= SNOVA_R= SNOVA_M1= SNOVA_M2="
	@echo "    XOF knob is AES=0|1 (SHAKE is the default); PK_EXPAND_SHAKE is derived and ignored."
	@echo ""
	@echo "  make kat-verify SNOVA_V=27 SNOVA_O=5 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=4 SNOVA_M1=5 SNOVA_M2=80"
	@echo "                 (default target; self-contained: regen .rsp + sha256 vs ref_kat/KAT_DIGESTS.sha256)"
	@echo "  make kat       ... KAT_RSP=<path> KAT_N=<n>"
	@echo "                 (byte-diff vs an EXTERNAL oracle .rsp; this tree ships none)"
	@echo "  make bench     SNOVA_V=... ... BENCH_N=<n>"
	@echo "  make PQCgenKAT SNOVA_V=... ... (gen .req/.rsp in build dir)"
	@echo ""
	@echo "WASM (emscripten):"
	@echo "  make wasm SNOVA_V=27 SNOVA_O=5 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=4 SNOVA_M1=5 SNOVA_M2=80"
	@echo "                                        (build .wasm at build/wasm_*/)"
	@echo "  make build_wasmapi SNOVA_V=27 SNOVA_O=5 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=4 SNOVA_M1=5 SNOVA_M2=80 WASM_OUT_PATH=<dir>/"
	@echo ""
	@echo "Other:"
	@echo "  make clean     - remove build/"
	@echo "  make help"
