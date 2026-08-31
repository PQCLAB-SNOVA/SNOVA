
#if !SNOVA_WRAPPER_STACK
static rct_pk_t rct_nist_pkx;
#endif

int snova_crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
    uint8_t seed[SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE];
    randombytes(seed, SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE);
    SNOVA_CT_POISON(seed + SEED_LENGTH_PUBLIC, SEED_LENGTH_PRIVATE);
    int rc = rct_genkeys(pk, sk, seed);
    SNOVA_CLEAR_OBJ(seed);
    return rc;
}

int snova_crypto_sign(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen, const uint8_t *sk) {
    uint8_t digest[BYTES_DIGEST], salt[BYTES_SALT];
    uint8_t sigbuf[CRYPTO_BYTES_R];
    shake256(digest, BYTES_DIGEST, m, mlen);
    randombytes(salt, BYTES_SALT);
    SNOVA_CT_POISON(salt, BYTES_SALT);
    int rc = rct_sign(sk, sigbuf, digest, BYTES_DIGEST, salt);
    if (rc) return rc;
    memmove(sm + CRYPTO_BYTES_R, m, mlen);
    memcpy(sm, sigbuf, CRYPTO_BYTES_R);
    *smlen = mlen + CRYPTO_BYTES_R;
    return 0;
}

int snova_crypto_sign_open(uint8_t *m, size_t *mlen, const uint8_t *sm, size_t smlen, const uint8_t *pk) {
    if (smlen < CRYPTO_BYTES_R) return -1;
    uint8_t digest[BYTES_DIGEST];
#if SNOVA_WRAPPER_STACK
    _Alignas(64) rct_pk_t pkx_l;
    rct_pk_t *pkx = &pkx_l;
#else
    rct_pk_t *pkx = &rct_nist_pkx;
#endif
    if (rct_pk_expand(pkx, pk)) return -1;
    shake256(digest, BYTES_DIGEST, sm + CRYPTO_BYTES_R, smlen - CRYPTO_BYTES_R);
    if (rct_verify(pkx, sm, digest, BYTES_DIGEST)) return -1;
    memmove(m, sm + CRYPTO_BYTES_R, smlen - CRYPTO_BYTES_R);
    *mlen = smlen - CRYPTO_BYTES_R;
    return 0;
}

int snova_keygen(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
    return rct_genkeys(pk, sk, seed);
}

size_t snova_skx_bytes(void) { return sizeof(rct_skx_t); }

void snova_sk_expand_skx(uint8_t *skx, const uint8_t *sk) {
    rct_sk_expand(sk, (rct_skx_t *)skx);
}

int snova_sign_digest_skx(uint8_t *sm, const uint8_t *digest, size_t dlen,
                          const uint8_t *salt, const uint8_t *skx) {
    return rct_sign_expanded((rct_skx_t *)(uintptr_t)skx, sm, digest, dlen, salt);
}

void snova_skx_clear(uint8_t *skx) {
    SNOVA_CLEAR(skx, sizeof(rct_skx_t));
}

size_t snova_pkx_bytes(void) { return sizeof(rct_pk_t); }

int snova_pk_expand(uint8_t *pkx, const uint8_t *pk) {
    return rct_pk_expand((rct_pk_t *)pkx, pk);
}

int snova_verify_digest_pkx(const uint8_t *digest, size_t dlen,
                            const uint8_t *sm, const uint8_t *pkx) {
    return rct_verify((const rct_pk_t *)pkx, sm, digest, dlen);
}
