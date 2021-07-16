#include "{prefix}matlab.h"
#include "{prefix}cvpdynmem.h"
#include "{prefix}cvpmatlab.h"

mxArray* {prefix}_cvp_matlab_create_prb(const struct {prefix}_cvp *cvp) {{
    mxArray* prb = mxCreateStructMatrix(1, 1, 0, 0);

{create_prb}
    return prb;
}}

mxArray* {prefix}_cvp_matlab_create_parameters(const struct {prefix}_cvp *cvp) {{
    mxArray* parameters = mxCreateStructMatrix(1, 1, 0, 0);

{create_parameters}
    return parameters;
}}

void {prefix}_cvp_matlab_apply_parameters(struct {prefix}_cvp *cvp, const mxArray* parameters) {{

{apply_parameters}
}}
