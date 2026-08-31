/**
 * @file primitives/snova_aes.c
 */
#include "snova_aes.h"
#include "platform.h"

#if defined(SNOVA_ARCH_X86_AVX2)
  #include "../platforms/x86_avx2/aes_ni.h"
#elif defined(SNOVA_ARCH_WASM_SIMD128)
  #include "../platforms/wasm_simd128/aes_wasm.h"
#elif defined(SNOVA_ARCH_PORTABLE_OPT)
  #include "../platforms/portable_opt/aes_opt.h"
#elif defined(SNOVA_ARCH_REF)
  #include "../platforms/ref/aes_ref.h"
#endif

void AES_256_ECB(const unsigned char *key, const uint8_t *input, unsigned char *output)
{
    uint8_t rk[240];
    aes256_keyexp(key, rk);
    aes256_encrypt(rk, input, output);
}

int AES_128_CTR(unsigned char *output, size_t outputByteLen,
                const unsigned char *input, size_t inputByteLen)
{
    (void)inputByteLen;
    aes128_ctr_zero(output, outputByteLen, input);
    return 0;
}
