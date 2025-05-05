SNOVA-LIBOQS
=======
This directory contains the files used by the integration of SNOVA into the liboqs library.

Download liboqs
```
git clone https://github.com/open-quantum-safe/liboqs
cd liboqs
git checkout d91aad0431176567c4ceafc1b4fed74b7645a1f1
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
