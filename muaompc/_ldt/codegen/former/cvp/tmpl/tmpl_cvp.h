#ifndef {PREFIX}FORMQP_H
#define {PREFIX}FORMQP_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"

{def_flags}

/* Declarations that depend only on the structure of the problem
 * These are used to generate code independent of the data size
 */

enum {{
{par_inds}
    {PREFIX}_PAR_NUM
}};

enum {{
    {PREFIX}_H,
    {PREFIX}_U_LB,
    {PREFIX}_U_UB,
    {PREFIX}_V,
    {PREFIX}_G,
    {PREFIX}_V_LB,
    {PREFIX}_V_UB,
    {PREFIX}_PMETRIC_NUM
}};

enum {{
    {PREFIX}_SOCC_WM,
    {PREFIX}_SOCC_WVT,
    {PREFIX}_SOCC_WN,
    {PREFIX}_SOCC_WS,
    {PREFIX}_SOCC_PMETRIC_NUM
}};

struct {prefix}_term {{
    uint32_t rows;
    uint32_t cols;
    real_t *data;
}};

struct {prefix}_pmetric {{
    uint32_t *fac_num;
    uint32_t *par_id;
    struct {prefix}_term *val;
    struct {prefix}_term *aux;
    struct {prefix}_term *fac0;
    struct {prefix}_term **fac;
    struct {prefix}_term **par;
}};

struct {prefix}_cvp_socc {{
    struct {prefix}_pmetric *pmetric[{PREFIX}_SOCC_PMETRIC_NUM];
}};

struct {prefix}_cvp_prb_socc {{
    struct {prefix}_term *Wm;
    struct {prefix}_term *wn;
    struct {prefix}_term *wvT;
    struct {prefix}_term *ws;
}};

struct {prefix}_cvp_prb {{
{prb_terms}
}};

/* Second-order cone program in its most general form,
       for a given set of parameters p.
 * The parametric second-order cone program to be solved has the form:
 * minimize 0.5 * u^T * H * u + u^T * g
 * subject to u_lb(p) <= u <= u_ub(p)
 *            v_lb(p) <= V * u <= v_ub(p)
 *            u^T * Wq_i * u <= wr_i(p), i=1,...,n_q
 *            | Wm_j * u + wn_j(p) | <= wvT_j * u + ws_j(p), i=1,...,n_s
 *
 * the transpose of a matrix is denote by ^T.
 *H;  The Hessian matrix.
 *g;  The gradient vector.
 *u_lb;   The lower bound for the box constraints.
 *u_ub;   The upper bound for the box constraints.
 *V;   The mixed constraints matrix.
 *v_lb;   The lower bound for the mixed constraints.
 *v_ub;   The upper bound for the mixed constraints.
 * TODO: update info
 */

struct {prefix}_cvp {{
    struct {prefix}_term *par[{PREFIX}_PAR_NUM];
    struct {prefix}_pmetric *pmetric[{PREFIX}_PMETRIC_NUM];
    uint32_t *socc_num;
    struct {prefix}_cvp_socc **socc;
    struct {prefix}_cvp_prb *prb;
}};

struct {prefix}_cvp_parameters {{
{parameters_terms}
}};

extern void {prefix}_cvp_copy_parameters(struct {prefix}_cvp *cvp, struct {prefix}_cvp_parameters *p);
extern void {prefix}_cvp_form_problem(struct {prefix}_cvp *cvp);

#ifdef __cplusplus
}}
#endif

#endif
