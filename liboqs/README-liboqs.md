SNOVA-LIBOQS
=======
This is an experimental integration of SNOVA into the liboqs library. It works but it will need to be updated if it is to comply with the liboqs coding conventions.

Building:

Download SNOVA with the liboqs integration
```
git clone https://github.com/vacuas/SNOVA-OQS SNOVA
export SNOVA_DIR=$(pwd)/SNOVA
```

Download liboqs
```
git clone https://github.com/open-quantum-safe/liboqs
cd liboqs
git checkout 6337a8424deae09aa835ddd920faff83a8f0e1d7
```

Add SNOVA to the algorithms known to liboqs:
```
cp $SNOVA_DIR/liboqs/copy_from_upstream_snova.yml scripts/copy_from_upstream/copy_from_upstream.yml
cp $SNOVA_DIR/liboqs/snova.yml docs/algorithms/sig
```

Create some files required for copy_from_upstream to work without errors (this appears to be a bug in copy_from_upstream):
```
mkdir -p src/sig/snova/
touch src/sig/snova/CMakeLists.txt src/sig/snova/sig_snova.h
touch src/sig/snova/sig_snova_SNOVA_24_5_4.c src/sig/snova/sig_snova_SNOVA_24_5_4_SHAKE.c
touch src/sig/snova/sig_snova_SNOVA_24_5_4_esk.c src/sig/snova/sig_snova_SNOVA_24_5_4_SHAKE_esk.c
touch src/sig/snova/sig_snova_SNOVA_37_17_2.c src/sig/snova/sig_snova_SNOVA_25_8_3.c
touch src/sig/snova/sig_snova_SNOVA_56_25_2.c  src/sig/snova/sig_snova_SNOVA_49_11_3.c
touch src/sig/snova/sig_snova_SNOVA_37_8_4.c src/sig/snova/sig_snova_SNOVA_24_5_5.c 
touch src/sig/snova/sig_snova_SNOVA_66_15_3.c src/sig/snova/sig_snova_SNOVA_60_10_4.c
touch src/sig/snova/sig_snova_SNOVA_29_6_5.c docs/algorithms/sig/snova.md
```

Generate the source using copy_from_upstream
```
export LIBOQS_DIR=$(pwd) ; cd scripts/copy_from_upstream ; python3 ./copy_from_upstream.py -k copy ; cd ../..
```

Build and optionally install liboqs:
```
rm -r build/ ; mkdir build/ && cd build/ ; cmake -GNinja .. ; ninja
ninja install
```

## OQS provider

Download and go to the root:
```
git clone https://github.com/open-quantum-safe/oqs-provider
cd oqs-provider
git checkout 1f190e87557b017530418ee3cc8b7d2fde3f740b
```

Add snova to the supported algorithms
```
export LIBOQS_SRC_DIR=$(pwd)/../liboqs
cp $SNOVA_DIR/liboqs/generate.yml oqs-template/
python3 oqs-template/generate.py
```

Build using the fullbuild script
```
export liboqs_DIR=$(pwd)/../liboqs.install
./scripts/fullbuild.sh
```

### OpenSSL OQS provider testing

Obtain info
```
openssl list -provider-path _build/lib/ -provider oqsprovider -providers -verbose
openssl list -signature-algorithms -provider-path _build/lib/ -provider oqsprovider
```

Create private and public keys
```
openssl genpkey -algorithm snova2454 -out private.pem  -provider-path _build/lib/ -provider oqsprovider
openssl pkey -in private.pem -pubout -out public.pem -provider-path _build/lib/ -provider oqsprovider
```

Sign README.md
```
openssl dgst -sha3-256 -sign private.pem -provider-path _build/lib/ -provider oqsprovider -out README.md.signature README.md
```

Verify using
```
openssl dgst -sha3-256 -verify public.pem -provider-path _build/lib/ -provider oqsprovider -signature README.md.signature README.md
```

### Certificates

Create keypair and view
```
openssl genpkey -algorithm snova2454 -out private.pem  -provider-path _build/lib/ -provider oqsprovider
openssl pkey -in private.pem -pubout -out public.pem -provider-path _build/lib/ -provider oqsprovider
```

Create sign request and view
```
openssl req -key private.pem -new -out domain.csr -provider-path _build/lib/ -provider oqsprovider
openssl req -in domain.csr -text -noout -provider-path _build/lib/ -provider oqsprovider
```

Self-sign and view
```
openssl x509 -signkey private.pem -in domain.csr -req -days 365 -out domain.crt -provider-path _build/lib/ -provider oqsprovider
openssl x509 -in domain.crt -text -noout -provider-path _build/lib/ -provider oqsprovider
```
