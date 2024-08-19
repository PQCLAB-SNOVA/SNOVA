#include "snova.h"

#include "gf16_matrix_inline.h"

#if OPTIMISATION != 0
    #define snova_shake snova_shake_opt
    #include "snova_kernel.h"
    #include "snova_plasma/snova_plasma_option.h"
#else
    #define snova_shake snova_shake_ref
    #include "snova_kernel.h"
    #define gen_F gen_F_ref
    #define gen_P22 gen_P22_ref
    #define sign_digest_core sign_digest_core_ref
    #define verify_core verify_signture_ref
#endif


/**
 * SNOVA init
 */
void snova_init(void) {
    static int first_time = 1;
    if (first_time) {
        first_time = 0;
        init_gf16_tables();
        gen_S_array();

#if OPTIMISATION != 0
        snova_plasma_init();
#endif    
    }
}

/**
 * generate snova key elements.
 * @param key_elems - pointer to output snova key elements.
 * @param pk_seed - pointer to input public key seed.
 * @param sk_seed - pointer to input private key elements.
 */
void generate_keys_core(snova_key_elems* key_elems, const uint8_t* pk_seed, const uint8_t* sk_seed) {
    gen_seeds_and_T12(key_elems->T12, sk_seed);
    memcpy(key_elems->pk.pt_public_key_seed, pk_seed, seed_length_public);
    gen_A_B_Q_P(&(key_elems->map1), pk_seed);
    gen_F(&(key_elems->map2), &(key_elems->map1), key_elems->T12);
    gen_P22(key_elems->pk.P22, key_elems->T12, key_elems->map1.P21, key_elems->map2.F12);
}

/**
 * Generates public and private key. where private key is the seed of private
 * key.
 * @param pk - pointer to output public key.
 * @param ssk - pointer to output private key.
 * @param pkseed - pointer to input public key seed.
 * @param skseed - pointer to input private key seed.
 */
void generate_keys_ssk(uint8_t* pk, uint8_t* ssk, const uint8_t* pkseed, const uint8_t* skseed) {
    snova_init();
    snova_key_elems key_elems;
    generate_keys_core(&key_elems, pkseed, skseed);
    pk_pack(pk, &key_elems);
    memcpy(ssk, pkseed, seed_length_public);
    memcpy(ssk + seed_length_public, skseed, seed_length_private);

    // Clear Secret!
    SNOVA_CLEAR_BYTE(&key_elems, sizeof(key_elems));
}

/**
 * Generates public and private key. where private key is the expanded version.
 * @param pk - pointer to output public key.
 * @param esk - pointer to output private key. (expanded)
 * @param pkseed - pointer to input public key seed.
 * @param skseed - pointer to input private key seed.
 */
void generate_keys_esk(uint8_t* pk, uint8_t* esk, const uint8_t* pkseed, const uint8_t* skseed) {
    snova_init();
    snova_key_elems key_elems;
    generate_keys_core(&key_elems, pkseed, skseed);
    pk_pack(pk, &key_elems);
    sk_pack(esk, &key_elems, skseed);
    
    // Clear Secret!
    SNOVA_CLEAR_BYTE(&key_elems, sizeof(key_elems));
}

/**
 * Generates public key. where private key is the seed of private
 * @param pk - pointer to output public key.
 * @param ssk - pointer to output private key.
 */
void generate_pk_with_ssk(uint8_t* pk, const uint8_t* ssk) {
    snova_init();
    const uint8_t* pkseed = ssk;
    const uint8_t* skseed = ssk + seed_length_public;
    snova_key_elems key_elems;

    generate_keys_core(&key_elems, pkseed, skseed);
    pk_pack(pk, &key_elems);

    // Clear Secret!
    SNOVA_CLEAR_BYTE(&key_elems, sizeof(key_elems));
}

/**
 * Generates public key. where private key is the expanded version.
 * @param pk - pointer to output public key.
 * @param @param esk - pointer to output private key. (expanded)
 */
void generate_pk_with_esk(uint8_t* pk, const uint8_t* esk) {
    snova_init();
    snova_key_elems key_elems;
    sk_gf16 sk_upk;
    
    sk_unpack(&sk_upk, esk);
    const uint8_t* pkseed = sk_upk.pt_public_key_seed;

    memcpy(key_elems.pk.pt_public_key_seed, pkseed, seed_length_public);
    gen_A_B_Q_P(&(key_elems.map1), pkseed);
    gen_F(&(key_elems.map2), &(key_elems.map1), sk_upk.T12);
    gen_P22(key_elems.pk.P22, sk_upk.T12, key_elems.map1.P21, key_elems.map2.F12);
    pk_pack(pk, &key_elems);
    
    // Clear Secret!
    SNOVA_CLEAR_BYTE(&key_elems, sizeof(key_elems));
}



/**
 * Compute the signature using ssk (private key seed). some preparatory work
 * before using sign_digest_core()
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param ssk - pointer to input private key (seed).
 */
void sign_digest_ssk(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, uint8_t* array_salt,
                     const uint8_t* ssk) {
    snova_init();
    snova_key_elems key_elems;
    const uint8_t* pk_seed = ssk;
    const uint8_t* sk_seed = ssk + seed_length_public;

    gen_seeds_and_T12(key_elems.T12, sk_seed);
    gen_A_B_Q_P(&(key_elems.map1), pk_seed);
    gen_F(&(key_elems.map2), &(key_elems.map1), key_elems.T12);
    sign_digest_core(pt_signature, digest, bytes_digest, array_salt, key_elems.map1.Aalpha, key_elems.map1.Balpha,
                     key_elems.map1.Qalpha1, key_elems.map1.Qalpha2, key_elems.T12, key_elems.map2.F11, key_elems.map2.F12,
                     key_elems.map2.F21, pk_seed, sk_seed);

    // Clear Secret!
    SNOVA_CLEAR_BYTE(&key_elems, sizeof(key_elems));
}

/**
 * Compute the signature using esk (). some preparatory work before using
 * sign_digest_core()
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param esk - pointer to input private key (expanded).
 */
void sign_digest_esk(uint8_t* pt_signature, const uint8_t* digest, uint64_t bytes_digest, uint8_t* array_salt,
                     const uint8_t* esk) {
    snova_init();
    sk_gf16 sk_upk;
    sk_unpack(&sk_upk, esk);
    sign_digest_core(pt_signature, digest, bytes_digest, array_salt, sk_upk.Aalpha, sk_upk.Balpha, sk_upk.Qalpha1,
                     sk_upk.Qalpha2, sk_upk.T12, sk_upk.F11, sk_upk.F12, sk_upk.F21, sk_upk.pt_public_key_seed,
                     sk_upk.pt_private_key_seed);

    // Clear Secret!
    SNOVA_CLEAR_BYTE(&sk_upk, sizeof(sk_upk));
}

/**
 * Verifies signature.
 * @param pt_digest - pointer to input digest.
 * @param pt_signature - pointer to output signature.
 * @param pk - pointer to output public key.
 * @returns - 0 if signature could be verified correctly and -1 otherwise
 */
int verify_signture(const uint8_t* pt_digest, uint64_t bytes_digest, const uint8_t* pt_signature, const uint8_t* pk) {
    snova_init();
    return verify_core(pt_digest, bytes_digest, pt_signature, pk);
}
