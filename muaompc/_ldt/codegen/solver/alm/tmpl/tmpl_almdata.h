#ifndef {PREFIX}{DATANAME}ALMDATA_H
#define {PREFIX}{DATANAME}ALMDATA_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}alm.h"

extern void {prefix}{dataname}_alm_setup_solver(
                struct {prefix}_alm *alm,
                struct {prefix}_cvp_prb *prb);

#ifdef __cplusplus
}}
#endif

#endif
