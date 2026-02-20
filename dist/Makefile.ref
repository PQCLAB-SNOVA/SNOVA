CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -O3
OBJS = ./build/snova_ref.o ./build/symmetric_ref.o ./build/rng.o ./build/sign.o ./build/aes.o
HDRS = api.h rng.h symmetric.h snova.h snova_params.h
# LIBS = -lcrypto


all: build kat

clean:
	rm -rf build *.req *.rsp *.log

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
