#ifndef {PREFIX}_CVPMATLAB_H
#define {PREFIX}_CVPMATLAB_H

#include "mex.h"
#include "{prefix}cvp.h"

extern mxArray* {prefix}_cvp_matlab_create_prb(const struct {prefix}_cvp *cvp);
extern mxArray* {prefix}_cvp_matlab_create_parameters(const struct {prefix}_cvp *cvp);
extern void {prefix}_cvp_matlab_apply_parameters(struct {prefix}_cvp *cvp, const mxArray* p);

#endif
