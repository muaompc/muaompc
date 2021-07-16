#include "{prefix}cvp.h"
#include "{prefix}fgm.h"

real_t {prefix}{dataname}_nu = {nu};
uint32_t {prefix}{dataname}_optvar_veclen = {optvar_veclen};

struct {prefix}_fgm_conf {prefix}{dataname}_fgm_conf = {{1, 0}};

real_t {prefix}{dataname}_u_opt[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_u_ini[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_fgm_tmp1_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_fgm_tmp2_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_fgm_tmp3_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_fgm_tmp4_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_fgm_tmp5_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_fgm_tmp6_optvar_seqlen[] = {zeros_optvar_seqlen};

void {prefix}{dataname}_fgm_setup_solver(
                struct {prefix}_fgm *fgm,
                struct {prefix}_cvp_prb *prb)
{{
    fgm->goL = prb->g->data;
    fgm->HoL = prb->H->data;
    fgm->u_lb = prb->u_lb->data;
    fgm->u_ub = prb->u_ub->data;
    fgm->optvar_seqlen = prb->u_lb->rows;
    fgm->sizeof_optvar_seqlen = sizeof(real_t) * prb->u_lb->rows;

    fgm->nu = &{prefix}{dataname}_nu;
    fgm->optvar_veclen = {prefix}{dataname}_optvar_veclen;

    fgm->conf = &{prefix}{dataname}_fgm_conf;
    fgm->j_in = &(fgm->conf->in_iter);

    fgm->u_opt = {prefix}{dataname}_u_opt;
    fgm->u_ini = {prefix}{dataname}_u_ini;
    fgm->tmp1_optvar_seqlen = {prefix}{dataname}_fgm_tmp1_optvar_seqlen;
    fgm->tmp2_optvar_seqlen = {prefix}{dataname}_fgm_tmp2_optvar_seqlen;
    fgm->tmp3_optvar_seqlen = {prefix}{dataname}_fgm_tmp3_optvar_seqlen;
    fgm->tmp4_optvar_seqlen = {prefix}{dataname}_fgm_tmp4_optvar_seqlen;
    fgm->tmp5_optvar_seqlen = {prefix}{dataname}_fgm_tmp5_optvar_seqlen;
    fgm->tmp6_optvar_seqlen = {prefix}{dataname}_fgm_tmp6_optvar_seqlen;

    return;
}}
