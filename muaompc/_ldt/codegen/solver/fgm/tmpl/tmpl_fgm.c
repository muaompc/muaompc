/** Quadratic program solver based on the fast gradient method.
 */

#include <string.h> /* sizeof */

#include "{prefix}mtxops.h"
#include "{prefix}fgm.h"

/* static functions declaration */
static void {prefix}_fgm_compute_projected_grad_step(
                const struct {prefix}_fgm *fgm,
                real_t u[], const real_t w[],
                const real_t gradoL[]);


/* external functions definition */

/* Minimize MPC quadratic program using fast gradient method */
void {prefix}_fgm_solve_problem(const struct {prefix}_fgm *fgm) {{
    real_t *u_old = fgm->tmp1_optvar_seqlen;
    real_t *w = fgm->tmp2_optvar_seqlen;
    real_t *gradoL = fgm->tmp3_optvar_seqlen;
    uint32_t j;

    memcpy(fgm->u_opt, fgm->u_ini, fgm->sizeof_optvar_seqlen);
    memcpy(w, fgm->u_ini, fgm->sizeof_optvar_seqlen);
    memcpy(u_old, fgm->u_ini, fgm->sizeof_optvar_seqlen);

    for (j = 0; j < *(fgm->j_in); j++) {{
        {prefix}_fgm_compute_grad_over_L(fgm, gradoL, w);
        {prefix}_fgm_minimize_qp_iteration(fgm, fgm->u_opt, u_old, w, gradoL);
    }}

    if (fgm->conf->warm_start) {{
        {prefix}_fgm_warm_start(fgm);
    }}

    return;
}}

void {prefix}_fgm_warm_start(const struct {prefix}_fgm *fgm) {{
    {prefix}_mtx_shift_sequence(fgm->u_ini, fgm->u_opt, fgm->optvar_veclen,
            fgm->optvar_seqlen);
    return;
}}

/* Execute one iteration of the minimization algorithm */
void {prefix}_fgm_minimize_qp_iteration(const struct {prefix}_fgm *fgm,
                real_t u[], real_t u_old[],
                real_t w[], const real_t gradoL[]) {{
    real_t *du = fgm->tmp5_optvar_seqlen;
    real_t *nu_du = fgm->tmp6_optvar_seqlen;

    {prefix}_fgm_compute_projected_grad_step(fgm, u, w, gradoL);
    {prefix}_mtx_substract(du, u, u_old, fgm->optvar_seqlen, 1);
    {prefix}_mtx_scale(nu_du, du, *(fgm->nu), fgm->optvar_seqlen, 1);
    {prefix}_mtx_add(w, u, nu_du, fgm->optvar_seqlen, 1);
    memcpy(u_old, u, fgm->sizeof_optvar_seqlen);

    return;
}}

/* Compute gradient divided by Lipschitz constant */
void {prefix}_fgm_compute_grad_over_L(const struct {prefix}_fgm *fgm, 
                real_t gradoL[], const real_t w[]) {{
    real_t *HoL_w = fgm->tmp4_optvar_seqlen;
    /* gradoL = (H/L) * w + (g/L) */
    {prefix}_mtx_multiply_mtx_vec(HoL_w, fgm->HoL, w, fgm->optvar_seqlen, fgm->optvar_seqlen);
    {prefix}_mtx_add(gradoL, HoL_w, fgm->goL, fgm->optvar_seqlen, 1);

    return;
}}

/* static functions definition */

static void {prefix}_fgm_compute_projected_grad_step(
                const struct {prefix}_fgm *fgm,
                real_t u[], const real_t w[], const real_t gradoL[]) {{
    {prefix}_mtx_substract(u, w, gradoL, fgm->optvar_seqlen, 1);
    {prefix}_mtx_saturate_vec(u, fgm->u_lb, fgm->u_ub, fgm->optvar_seqlen);

    return;
}}
