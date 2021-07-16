#include "{prefix}matlab.h"
#include "{prefix}almdynmem.h"
#include "{prefix}almmatlab.h"

mxArray* {prefix}_alm_matlab_create_conf(const struct {prefix}_alm *alm)
{{
    mxArray* conf = mxCreateStructMatrix(1, 1, 0, 0);
    {prefix}_matlab_add_scalar(conf, "in_iter", (real_t)alm->conf->in_iter);
    {prefix}_matlab_add_scalar(conf, "ex_iter", (real_t)alm->conf->ex_iter);
    {prefix}_matlab_add_scalar(conf, "warm_start", (real_t)alm->conf->warm_start);

    return conf;
}}

void {prefix}_alm_matlab_apply_conf(struct {prefix}_alm *alm, const mxArray* conf)
{{
    {prefix}_matlab_read_uint(conf, "in_iter", &(alm->conf->in_iter));
    {prefix}_matlab_read_uint(conf, "ex_iter", &(alm->conf->ex_iter));
    {prefix}_matlab_read_uint(conf, "warm_start", &(alm->conf->warm_start));

    return;
}}
