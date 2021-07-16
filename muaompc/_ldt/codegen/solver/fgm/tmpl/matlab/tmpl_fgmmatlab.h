#ifndef {PREFIX}_FGMMATLAB_H
#define {PREFIX}_FGMMATLAB_H

#include "mex.h"
#include "{prefix}fgm.h"

extern mxArray* {prefix}_fgm_matlab_create_conf(const struct {prefix}_fgm *fgm);
extern void {prefix}_fgm_matlab_apply_conf(struct {prefix}_fgm *fgm,
                const mxArray* conf);

#endif
