#include "snova.h"

#if SNOVA_q == 16
#if SNOVA_l == 4
#include "snova_opt_16.c"
#else
#include "snova_avx2_16.c"
#endif
#else
#include "snova_opt_q.c"
#endif
