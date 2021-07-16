#include "{prefix}matlab.h"
#include "{prefix}ctlmatlab.h"
#include "{prefix}ctldynmem.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{{
    {prefix}_dynmem_error_t ret;
    struct {prefix}_ctl *data;
    char *fname;

    if (nrhs != 1) {{
        mexErrMsgTxt("Function accepts only one input arguments.");
        return;
    }}

    if (nlhs != 1) {{
        mexErrMsgTxt("Function accepts only one output variable.");
        return;
    }}

    data = {prefix}_ctl_allocate_ctl();
    if (NULL == data) {{
        mexErrMsgTxt("Could not allocate memory for data.");
        return;
    }}

    fname = mxArrayToString(prhs[0]);
    ret = {prefix}_ctl_setup_ctl(data, fname);
    mxFree(fname);
    if ({PREFIX}_DYNMEM_OK != ret) {{
        mexErrMsgTxt("Could not parse data.");
        return;
    }}

    plhs[0] = {prefix}_matlab_get_ctl_data(data);
    return;
}}

