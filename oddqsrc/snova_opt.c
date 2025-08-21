#include "snova.h"

#if SNOVA_q == 16
#include "snova_opt_16.c"
#else
#include "snova_opt_q.c"
#endif
