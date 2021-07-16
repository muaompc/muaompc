/** Solve state constrained (stc) Model Predictive Control (MPC) problem.
 *
 */

#include <string.h> /* sizeof */

#include "{prefix}mtxops.h"
#include "{prefix}fgm.h"
#include "{prefix}alm.h"


/* Declaration of static functions */

static void {prefix}_alm_fgm_minimize_qp(const struct {prefix}_alm *alm,
                real_t u[], const real_t u_i[], const real_t l_i[]);
static void {prefix}_alm_fgm_compute_grad_over_L(const struct {prefix}_alm *alm,
                real_t gradoL[], const real_t w[], const real_t l_i[]);
static void {prefix}_alm_fgm_compute_state_constr_violation_penalty(
                const struct {prefix}_alm *alm, real_t out[], const real_t w[],
                const real_t l_i[], const real_t v_b[]);
static void {prefix}_alm_fgm_compute_l_update(const struct {prefix}_alm *alm,
                real_t l_i1[], const real_t w[], const real_t l_i[]);
static void {prefix}_alm_fgm_compute_grad_alm(const struct {prefix}_alm *alm,
                real_t out[], const real_t w[], const real_t l_i[]);
static void {prefix}_alm_warm_start(struct {prefix}_alm *alm);


/* Definition of external functions */

/* Minimize MPC quadratic program using augmented Lagrangian method
 * together with fast gradient method */
void {prefix}_alm_solve_problem(struct {prefix}_alm *alm)
{{

    real_t *u_i = alm->tmp1_optvar_seqlen; /* iteration update of input */
    real_t *l_i = alm->tmp1_lagmul_seqlen; /* it. update of Lagrange mult.*/
    uint32_t i;

    memcpy(u_i, alm->u_ini, alm->sizeof_optvar_seqlen);
    memcpy(l_i, alm->l_ini, alm->sizeof_lagmul_seqlen);

    for (i = 0; i < *(alm->i_ex); i++) {{
        {prefix}_alm_fgm_minimize_qp(alm, alm->u_opt, u_i, l_i);
        {prefix}_alm_fgm_compute_l_update(alm, alm->l_opt, alm->u_opt, l_i);
        /* do not take the next two outside the loop. The two previous
         * functions expect different pointers for outputs u/l_upd
         * and inputs u_i/l_i respectively. */
        memcpy(u_i, alm->u_opt, alm->sizeof_optvar_seqlen);
        memcpy(l_i, alm->l_opt, alm->sizeof_lagmul_seqlen);
    }}

  if (alm->conf->warm_start) {{
    {prefix}_alm_warm_start(alm);
  }}

    return;
}}


/* Definition of static functions */

/* Minimize ALM internal iteration using a fast gradient method */
static void {prefix}_alm_fgm_minimize_qp(const struct {prefix}_alm *alm, real_t u[],
        const real_t u_i[], const real_t l_i[])
{{
    real_t *u_old = alm->tmp2_optvar_seqlen;
    real_t *w = alm->tmp3_optvar_seqlen;
    real_t *gradoL = alm->tmp4_optvar_seqlen;
    uint32_t i;

    memcpy(u, u_i, alm->sizeof_optvar_seqlen);
    memcpy(w, u, alm->sizeof_optvar_seqlen);
    memcpy(u_old, u, alm->sizeof_optvar_seqlen);

    for (i = 0; i < *(alm->fgm->j_in); i++) {{
        {prefix}_alm_fgm_compute_grad_over_L(alm, gradoL, w, l_i);
        {prefix}_fgm_minimize_qp_iteration(alm->fgm, u, u_old, w, gradoL);
    }}

    return;
}}

/* Compute gradient divided by Lipschitz constant */
static void {prefix}_alm_fgm_compute_grad_over_L(const struct {prefix}_alm *alm,
        real_t gradoL[], const real_t w[], const real_t l_i[])
{{
    real_t *gradoL_inc = alm->tmp5_optvar_seqlen; /* input constrained grad. over L*/
    real_t *grad_stc = alm->tmp6_optvar_seqlen; /* state constrained grad. */
    real_t *gradoL_stc = alm->tmp7_optvar_seqlen; /* state constrained grad. over L*/

    {prefix}_fgm_compute_grad_over_L(alm->fgm, gradoL_inc, w);
    {prefix}_alm_fgm_compute_grad_alm(alm, grad_stc, w, l_i);
    {prefix}_mtx_scale(gradoL_stc, grad_stc, *(alm->Linv), alm->optvar_seqlen, 1);
    {prefix}_mtx_add(gradoL, gradoL_inc, gradoL_stc, alm->optvar_seqlen, 1);

    return;
}}

static void {prefix}_alm_fgm_compute_state_constr_violation_penalty(
        const struct {prefix}_alm *alm, real_t out[], const real_t w[],
        const real_t l_i[], const real_t v_b[])
{{
    /* penalty = l_i + mu * (V * w - v_b) */
    real_t *V_w = alm->tmp2_lagmul_seqlen;
    real_t *diff = alm->tmp3_lagmul_seqlen;
    real_t *pen_diff = alm->tmp4_lagmul_seqlen;

    {prefix}_mtx_multiply_mtx_vec(V_w, alm->V, w, alm->lagmul_seqlen, alm->optvar_seqlen);
    {prefix}_mtx_substract(diff, V_w, v_b, alm->lagmul_seqlen, 1);
    {prefix}_mtx_scale(pen_diff, diff, *(alm->mu), alm->lagmul_seqlen, 1);
    {prefix}_mtx_add(out, l_i, pen_diff, alm->lagmul_seqlen, 1);

    return;
}}

/* Compute multiplier update */
static void {prefix}_alm_fgm_compute_l_update(const struct {prefix}_alm *alm,
        real_t l_i1[], const real_t w[], const real_t l_i[])
{{
    real_t *p_pos = alm->tmp5_lagmul_seqlen; /* positive penalty */
    real_t *p_neg = alm->tmp6_lagmul_seqlen; /* negative penalty */

    {prefix}_alm_fgm_compute_state_constr_violation_penalty(alm, p_pos, w, l_i, alm->v_ub);
    {prefix}_alm_fgm_compute_state_constr_violation_penalty(alm, p_neg, w, l_i, alm->v_lb);

    {prefix}_mtx_max_vec_zero(p_pos, alm->lagmul_seqlen);
    {prefix}_mtx_min_vec_zero(p_neg, alm->lagmul_seqlen);

    {prefix}_mtx_add(l_i1, p_pos, p_neg, alm->lagmul_seqlen, 1);

    return;
}}

static void {prefix}_alm_fgm_compute_grad_alm(const struct {prefix}_alm *alm,
        real_t out[], const real_t w[], const real_t l_i[])
{{
    real_t *l_grad = alm->tmp7_lagmul_seqlen;
    real_t *V_T = alm->tmp1_lagmul_optvar_seqlen;

    {prefix}_alm_fgm_compute_l_update(alm, l_grad, w, l_i);
    {prefix}_mtx_transpose(V_T, alm->V, alm->lagmul_seqlen, alm->optvar_seqlen);
    {prefix}_mtx_multiply_mtx_vec(out, V_T, l_grad, alm->optvar_seqlen, alm->lagmul_seqlen);

    return;
}}

/* Warm start state constrained MPC */
static void {prefix}_alm_warm_start(struct {prefix}_alm *alm) {{
    {prefix}_fgm_warm_start(alm->fgm);
    {prefix}_mtx_shift_sequence(alm->l_ini, alm->l_opt, alm->lagmul_veclen, alm->lagmul_seqlen);

    return;
}}

