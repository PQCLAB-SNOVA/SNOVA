SNOVA-LIBOQS
=======
SNOVA has been merged into liboqs and oqs-provider. See https://github.com/open-quantum-safe/liboqs and https://github.com/open-quantum-safe/oqs-provider

This directory contains the glue code required to use SNOVA in liboqs.

### Support limitations

This project is not commercially supported. The software is maintained on a best-effort basis and this may change at any time. We will attempt to respond when able but offer no timing commitment whatsoever. 
Support contact: SNOVA team, pqclaborg@gmail.com

### OpenSSL OQS provider testing

Obtain info
```
openssl list -provider-path _build/lib/ -provider oqsprovider -providers -verbose
openssl list -signature-algorithms -provider-path _build/lib/ -provider oqsprovider
```

Create private and public keys
```
openssl genpkey -algorithm snova1b -out private.pem  -provider-path _build/lib/ -provider oqsprovider -provider default
openssl pkey -in private.pem -pubout -out public.pem -provider-path _build/lib/ -provider oqsprovider -provider default
```

Sign README.md
```
openssl dgst -sha3-256 -sign private.pem -provider-path _build/lib/ -provider oqsprovider -out sign.README.md README.md
```

Verify using
```
openssl dgst -sha3-256 -verify public.pem -provider-path _build/lib/ -provider oqsprovider -signature sign.README.md README.md
```

### Certificates

Create keypair and view
```
openssl genpkey -algorithm snova1s -out private.pem  -provider-path _build/lib/ -provider oqsprovider -provider default
openssl pkey -in private.pem -pubout -out public.pem -provider-path _build/lib/ -provider oqsprovider -provider default
```

Create sign request and view
```
openssl req -key private.pem -new -out domain.csr -provider-path _build/lib/ -provider oqsprovider -provider default
openssl req -in domain.csr -text -noout -provider-path _build/lib/ -provider oqsprovider -provider default
```

Self-sign and view
```
openssl x509 -signkey private.pem -in domain.csr -req -days 365 -out domain.crt -provider-path _build/lib/ -provider oqsprovider -provider default
openssl x509 -in domain.crt -text -noout -provider-path _build/lib/ -provider oqsprovider -provider default
```

Verify the self-signed certificate:
```
openssl verify -CAfile domain.crt -provider-path _build/lib/ -provider oqsprovider -provider default domain.crt
```
