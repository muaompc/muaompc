#ifndef {PREFIX}_ALMMATLAB_H
#define {PREFIX}_ALMMATLAB_H

#include "mex.h"
#include "{prefix}alm.h"

extern mxArray* {prefix}_alm_matlab_create_conf(const struct {prefix}_alm *alm);
extern void {prefix}_alm_matlab_apply_conf(struct {prefix}_alm *alm,
                const mxArray* conf);

#endif
