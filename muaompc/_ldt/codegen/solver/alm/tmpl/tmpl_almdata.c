#include "{prefix}cvp.h"
#include "{prefix}fgm.h"
#include "{prefix}alm.h"
#include "{prefix}{dataname}fgmdata.h"

struct {prefix}_fgm {prefix}{dataname}_fgm;

real_t {prefix}{dataname}_mu = {mu};
real_t {prefix}{dataname}_Linv = {Linv};
uint32_t {prefix}{dataname}_lagmul_veclen = {lagmul_veclen};
struct {prefix}_alm_conf {prefix}{dataname}_alm_conf = {{1, 1, 0}};

real_t {prefix}{dataname}_l_opt[] = {zeros_lagmul_seqlen};
real_t {prefix}{dataname}_l_ini[] = {zeros_lagmul_seqlen};

real_t {prefix}{dataname}_alm_tmp1_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_alm_tmp2_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_alm_tmp3_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_alm_tmp4_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_alm_tmp5_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_alm_tmp6_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_alm_tmp7_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_alm_tmp1_lagmul_seqlen[] = {zeros_lagmul_seqlen};
real_t {prefix}{dataname}_alm_tmp2_lagmul_seqlen[] = {zeros_lagmul_seqlen};
real_t {prefix}{dataname}_alm_tmp3_lagmul_seqlen[] = {zeros_lagmul_seqlen};
real_t {prefix}{dataname}_alm_tmp4_lagmul_seqlen[] = {zeros_lagmul_seqlen};
real_t {prefix}{dataname}_alm_tmp5_lagmul_seqlen[] = {zeros_lagmul_seqlen};
real_t {prefix}{dataname}_alm_tmp6_lagmul_seqlen[] = {zeros_lagmul_seqlen};
real_t {prefix}{dataname}_alm_tmp7_lagmul_seqlen[] = {zeros_lagmul_seqlen};
real_t {prefix}{dataname}_alm_tmp1_lagmul_optvar_seqlen[] = {zeros_lagmul_optvar_seqlen};

void {prefix}{dataname}_alm_setup_solver(
                struct {prefix}_alm *alm,
                struct {prefix}_cvp_prb *prb)
{{
    {prefix}{dataname}_fgm_setup_solver(&{prefix}{dataname}_fgm, prb);
    alm->fgm = &{prefix}{dataname}_fgm;

    alm->V = prb->V->data;
    alm->v_lb = prb->v_lb->data;
    alm->v_ub = prb->v_ub->data;
    alm->optvar_seqlen = prb->u_lb->rows;
    alm->lagmul_seqlen = prb->v_lb->rows;
    alm->sizeof_optvar_seqlen = sizeof(real_t) * prb->u_lb->rows;
    alm->sizeof_lagmul_seqlen = sizeof(real_t) * prb->v_lb->rows;

    alm->mu = &{prefix}{dataname}_mu;
    alm->Linv = &{prefix}{dataname}_Linv;
    alm->lagmul_veclen = {prefix}{dataname}_lagmul_veclen;

    alm->conf = &{prefix}{dataname}_alm_conf;
    alm->i_ex = &(alm->conf->ex_iter);
    alm->fgm->j_in = &(alm->conf->in_iter);

    alm->u_ini = alm->fgm->u_ini;
    alm->u_opt = alm->fgm->u_opt;
    alm->l_ini = {prefix}{dataname}_l_ini;
    alm->l_opt = {prefix}{dataname}_l_opt;

    alm->tmp1_optvar_seqlen = {prefix}{dataname}_alm_tmp1_optvar_seqlen;
    alm->tmp2_optvar_seqlen = {prefix}{dataname}_alm_tmp2_optvar_seqlen;
    alm->tmp3_optvar_seqlen = {prefix}{dataname}_alm_tmp3_optvar_seqlen;
    alm->tmp4_optvar_seqlen = {prefix}{dataname}_alm_tmp4_optvar_seqlen;
    alm->tmp5_optvar_seqlen = {prefix}{dataname}_alm_tmp5_optvar_seqlen;
    alm->tmp6_optvar_seqlen = {prefix}{dataname}_alm_tmp6_optvar_seqlen;
    alm->tmp7_optvar_seqlen = {prefix}{dataname}_alm_tmp7_optvar_seqlen;

    alm->tmp1_lagmul_seqlen = {prefix}{dataname}_alm_tmp1_lagmul_seqlen;
    alm->tmp2_lagmul_seqlen = {prefix}{dataname}_alm_tmp2_lagmul_seqlen;
    alm->tmp3_lagmul_seqlen = {prefix}{dataname}_alm_tmp3_lagmul_seqlen;
    alm->tmp4_lagmul_seqlen = {prefix}{dataname}_alm_tmp4_lagmul_seqlen;
    alm->tmp5_lagmul_seqlen = {prefix}{dataname}_alm_tmp5_lagmul_seqlen;
    alm->tmp6_lagmul_seqlen = {prefix}{dataname}_alm_tmp6_lagmul_seqlen;
    alm->tmp7_lagmul_seqlen = {prefix}{dataname}_alm_tmp7_lagmul_seqlen;

    alm->tmp1_lagmul_optvar_seqlen = {prefix}{dataname}_alm_tmp1_lagmul_optvar_seqlen;
    return;
}}
