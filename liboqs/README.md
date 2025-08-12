SNOVA-LIBOQS
=======
SNOVA has been merged into liboqs and oqs-provider. See https://github.com/open-quantum-safe/liboqs and https://github.com/open-quantum-safe/oqs-provider

This directory contains the glue code required to use SNOVA in liboqs.

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
