#ifndef {PREFIX}{DATANAME}PBMDATA_H
#define {PREFIX}{DATANAME}PBMDATA_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}pbm.h"
#include "{prefix}cvp.h"
#include "{prefix}mtxops.h"

extern void {prefix}{dataname}_pbm_setup_solver(
                struct {prefix}_pbm *pbm,
                struct {prefix}_cvp_prb *prb);

#ifdef __cplusplus
}}
#endif

#endif /* {PREFIX}{DATANAME}PBMDATA_H */