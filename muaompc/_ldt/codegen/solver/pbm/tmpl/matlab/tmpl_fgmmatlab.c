#include "{prefix}matlab.h"
#include "{prefix}fgmdynmem.h"
#include "{prefix}fgmmatlab.h"

mxArray* {prefix}_fgm_matlab_create_conf(const struct {prefix}_fgm *fgm)
{{
    mxArray* conf = mxCreateStructMatrix(1, 1, 0, 0);
    {prefix}_matlab_add_scalar(conf, "in_iter", (real_t)fgm->conf->in_iter);
    {prefix}_matlab_add_scalar(conf, "warm_start", (real_t)fgm->conf->warm_start);

    return conf;
}}

void {prefix}_fgm_matlab_apply_conf(struct {prefix}_fgm *fgm, const mxArray* conf)
{{
    {prefix}_matlab_read_uint(conf, "in_iter", &(fgm->conf->in_iter));
    {prefix}_matlab_read_uint(conf, "warm_start", &(fgm->conf->warm_start));

    return;
}}
