#ifndef STOCMPCFORMQP_H
#define STOCMPCFORMQP_H

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"

/* Declarations that depend only on the structure of the problem
 * These are used to generate code independent of the data size
 */
enum {
    STOCMPC_X_K = 0,

    STOCMPC_PAR_NUM
};

enum {
    STOCMPC_U_LB = 0,
    STOCMPC_H = 1,
    STOCMPC_WM1 = 2,
    STOCMPC_WM0 = 3,
    STOCMPC_WVT1 = 4,
    STOCMPC_WVT0 = 5,
    STOCMPC_U_UB = 6,

    STOCMPC_CONSTANT_NUM
};

enum {
    STOCMPC_WS1 = 0,
    STOCMPC_WS0 = 1,
    STOCMPC_WN0 = 2,
    STOCMPC_WN1 = 3,
    STOCMPC_G = 4,

    STOCMPC_PMETRIC_NUM
};

struct stocmpc_term {
    uint32_t rows;
    uint32_t cols;
    real_t *data;
};

struct stocmpc_pmetric {
    uint32_t *fac_num;
    uint32_t *par_id;
    struct stocmpc_term *val;
    struct stocmpc_term *aux;
    struct stocmpc_term *fac0;
    struct stocmpc_term **fac;
    struct stocmpc_term **par;
};

struct stocmpc_cvp_prb {
    struct stocmpc_term *ws1;
    struct stocmpc_term *ws0;
    struct stocmpc_term *wn0;
    struct stocmpc_term *wn1;
    struct stocmpc_term *g;
    struct stocmpc_term *u_lb;
    struct stocmpc_term *H;
    struct stocmpc_term *Wm1;
    struct stocmpc_term *Wm0;
    struct stocmpc_term *wvT1;
    struct stocmpc_term *wvT0;
    struct stocmpc_term *u_ub;

#if 0  /* left as reference */
    struct stocmpc_term *H;  /**< The Hessian matrix. */
    struct stocmpc_term *g;  /**< The gradient vector. */
    struct stocmpc_term *u_lb;  /**< The lower bound for the box constraints. */
    struct stocmpc_term *u_ub;  /**< The upper bound for the box constraints. */
    struct stocmpc_term *V;  /**< The mixed constraints matrix. */
    struct stocmpc_term *v_lb;  /**< The lower bound for the mixed constraints. */
    struct stocmpc_term *v_ub;  /**< The upper bound for the mixed constraints. */
#endif
};  /**< STOCMPC quadratic program form for a given system state x. 
 * The quadratic program to solve has the form:
 * minimize 0.5 * u^T * HoL * u + u^T * goL
 * subject to u_lb <= u <= u_ub
 *        v_lb <= V * x <= v_ub
 *
 * the transpose of a matrix is denote by ^T.
 */

struct stocmpc_cvp {
    struct stocmpc_term *par[STOCMPC_PAR_NUM];
    struct stocmpc_term *constant[STOCMPC_CONSTANT_NUM];
    struct stocmpc_pmetric *pmetric[STOCMPC_PMETRIC_NUM];
    struct stocmpc_cvp_prb *prb;
};

struct stocmpc_cvp_parameters {
    real_t *x_k;

};

extern void stocmpc_cvp_copy_parameters(struct stocmpc_cvp *cvp, struct stocmpc_cvp_parameters *p);
extern void stocmpc_cvp_form_problem(struct stocmpc_cvp *cvp);

#endif
