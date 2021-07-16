#include "{prefix}matlab.h"
#include "{prefix}ctl.h"
#include "{prefix}ctlmatlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{{
    const mxArray* ctl;
    struct {prefix}_ctl *data;

    if (nrhs != 1) {{
        mexErrMsgTxt("Function needs one input variable.");
        return;
    }}

    if (!mxIsStruct(prhs[0])) {{
        mexErrMsgTxt("Input argument has to be a mpc structure.");
        return;
    }}

    if (nlhs != 1) {{
        mexErrMsgTxt("Function accepts only one output variable.");
        return;
    }}

    ctl = prhs[0];
    data = {prefix}_matlab_get_data_from_ctl(ctl);
    {prefix}_matlab_apply_ctl(data, ctl);

    {prefix}_ctl_form_problem(data);

    plhs[0] = {prefix}_matlab_get_ctl_data(data);

    return;
}}
