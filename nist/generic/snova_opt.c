#include "snova.h"

#if SNOVA_q == 16
#if SNOVA_l == 2
#include "snova_opt_16_2.c"
#elif SNOVA_l == 5
#include "snova_opt_16_5.c"
#else
#include "snova_opt_16.c"
#endif
#elif defined(SNOVA_r) && (SNOVA_r != SNOVA_l)
#include "snova_opt_q_r.c"
#else
#include "snova_opt_q_s.c"
#endif
