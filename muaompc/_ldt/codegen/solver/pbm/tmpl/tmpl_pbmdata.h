#ifndef {PREFIX}{DATANAME}PBMDATA_H
#define {PREFIX}{DATANAME}PBMDATA_H

#include "{prefix}pbm.h"
#include "{prefix}cvp.h"
#include "{prefix}mtxops.h"

extern void {prefix}{dataname}_pbm_setup_solver(
                struct {prefix}_pbm *pbm,
                struct {prefix}_cvp_prb *prb);

#endif /* {PREFIX}{DATANAME}PBMDATA_H */