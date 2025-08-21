ifeq ($(ARCH),)
ARCH = native
endif

ifeq ($(OPT),)
# Default to OPT
OPT = OPT
endif

CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -O3 -march=$(ARCH) -mtune=$(ARCH)
# CFLAGS += -g -fstack-usage

# Select appropriate source file
ifeq ($(OPT), REF)
OBJS = ./build/snova_ref.o ./build/symmetric_ref.o
else ifeq ($(OPT), OPT)
OBJS = ./build/snova_opt.o ./build/symmetric.o
else
$(error OPT=$(OPT) is not supported)
endif

OBJS += ./build/rng.o ./build/sign.o
HDRS = api.h rng.h symmetric.h snova.h snova_params.h
LIBS = -lcrypto


all: build kat speed

clean:
	rm -rf build *.req *.rsp

build:
	mkdir -p build

./build/%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $(P) -o $@ -c $<

./build/PQCgenKAT_sign: $(OBJS) ./build/PQCgenKAT_sign.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

kat: build ./build/PQCgenKAT_sign
	./build/PQCgenKAT_sign

./build/speed: $(OBJS) ./build/speed.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

speed: build ./build/speed
	@./build/speed
