#ifndef MC04TYPES_H
#define MC04TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

/* MISRA C 2004 compliant numeric typedef */
#ifndef USE_MPCSTDINT
#include <stdint.h>
#else
#include "mpcstdint.h"
#endif

typedef float float32_t;
typedef double float64_t;
typedef long double float128_t;

#ifdef __cplusplus
}
#endif

#endif
