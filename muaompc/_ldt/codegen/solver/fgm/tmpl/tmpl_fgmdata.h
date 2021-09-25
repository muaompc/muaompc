#ifndef {PREFIX}{DATANAME}FGMDATA_H
#define {PREFIX}{DATANAME}FGMDATA_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}fgm.h"

extern void {prefix}{dataname}_fgm_setup_solver(
                struct {prefix}_fgm *fgm,
                struct {prefix}_cvp_prb *prb);

#ifdef __cplusplus
}}
#endif

#endif
