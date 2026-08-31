/**
 * @file snova.c
 */
#include "../snova.h"
#include "../gf16_core/gf16.h"
#include "../gf16_core/gf16m.h"
#include "drbg.h"
#include "../primitives/sym_aes.h"
#include "../primitives/sym_shake.h"

#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#include "ct_poison.h"
#include "secure_clear.h"

#if defined(SNOVA_DUDECT_RETRYCOUNT) && SNOVA_DUDECT_RETRYCOUNT
#include <stdint.h>
volatile uint32_t snova_dudect_sign_attempts = 0;
#define SNOVA_DUDECT_RETRY_RESET() (snova_dudect_sign_attempts = 0u)
#define SNOVA_DUDECT_RETRY_TICK()  (snova_dudect_sign_attempts++)
#else
#define SNOVA_DUDECT_RETRY_RESET() ((void)0)
#define SNOVA_DUDECT_RETRY_TICK()  ((void)0)
#endif

#include "snova_rect.h"
