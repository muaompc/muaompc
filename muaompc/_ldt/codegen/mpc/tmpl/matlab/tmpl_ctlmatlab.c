#include "{prefix}matlab.h"
#include "{prefix}{former}matlab.h"
#include "{prefix}{solver}matlab.h"
#include "{prefix}ctl.h"

static mxArray* {prefix}_matlab_create_data_handle(const struct {prefix}_ctl *data);

mxArray* {prefix}_matlab_get_ctl_data(const struct {prefix}_ctl *data)
{{
    mxArray* ctl;
    mxArray* prb;
    mxArray* parameters;
    mxArray* conf;
    mxArray* handle;

    ctl = mxCreateStructMatrix(1, 1, 0, 0);

    parameters = {prefix}_{former}_matlab_create_parameters(data->former);
    mxAddField(ctl, "parameters");
    mxSetField(ctl, 0, "parameters", parameters);

    conf = {prefix}_{solver}_matlab_create_conf(data->solver);
    mxAddField(ctl, "conf");
    mxSetField(ctl, 0, "conf", conf);

    prb = {prefix}_cvp_matlab_create_prb(data->former);
    mxAddField(ctl, "prb");
    mxSetField(ctl, 0, "prb", prb);

    {prefix}_matlab_add_matrix(ctl, "u_opt", data->solver->u_opt,
        data->solver->optvar_seqlen, 1);

    handle = {prefix}_matlab_create_data_handle(data);
    mxAddField(ctl, "handle");
    mxSetField(ctl, 0, "handle", handle);
    return ctl;
}}

void {prefix}_matlab_apply_ctl(struct {prefix}_ctl *data, const mxArray* ctl)
{{

    mxArray* parameters;
    mxArray* conf;
    parameters = mxGetField(ctl, 0, "parameters");
    conf = mxGetField(ctl, 0, "conf");
    {prefix}_{former}_matlab_apply_parameters(data->former, parameters);
    {prefix}_{solver}_matlab_apply_conf(data->solver, conf);
    {prefix}_matlab_read_matrix(ctl, "u_opt", data->solver->u_opt, data->solver->optvar_seqlen, 1);
}}

struct {prefix}_ctl* {prefix}_matlab_get_data_from_ctl(const mxArray* ctl)
{{
    mxArray* handle = mxGetField(ctl, 0, "handle");
    struct {prefix}_ctl *data = (struct {prefix}_ctl *)((uint64_t *)mxGetData(handle))[0];

    return data;
}}

static mxArray* {prefix}_matlab_create_data_handle(const struct {prefix}_ctl *data)
{{
    mxArray* handle  = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    ((uint64_t *)mxGetData(handle))[0] = (uint64_t)data;
    return handle;
}}

