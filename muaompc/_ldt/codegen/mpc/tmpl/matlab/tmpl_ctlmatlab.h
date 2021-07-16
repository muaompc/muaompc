#ifndef {PREFIX}CTLMATLAB_H
#define {PREFIX}CTLMATLAB_H

#include "mex.h"
#include "{prefix}ctl.h"

extern mxArray* {prefix}_matlab_get_ctl_data(const struct {prefix}_ctl *data);
extern void {prefix}_matlab_apply_ctl(struct {prefix}_ctl *data, const mxArray* ctl);
extern struct {prefix}_ctl* {prefix}_matlab_get_data_from_ctl(const mxArray* ctl);

#endif /* {PREFIX}CTLMATLAB_H */
