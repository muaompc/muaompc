#ifndef {PREFIX}_FGM_H
#define {PREFIX}_FGM_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "mc04types.h"
#include "arithmetic.h"

/* Quadratic Programming online solver using Fast gradients.
 * Most of these functions are intended for internal use.
 */
struct {prefix}_fgm_conf {{
uint32_t in_iter;  /**< Maximum number of internal loop (FGM) iterations (j_in). */
uint32_t warm_start;  /**< If not 0, automatically warm start next algorithm iteration. */
}}; /**< Configuration parameters of the {PREFIX} algorithm. */

struct {prefix}_fgm {{
struct {prefix}_fgm_conf *conf;  /**< Algorithm configuration data. */
real_t *u_opt;  /**< Solution to the optimal control problem. */
real_t *u_ini;  /**< Initial guess for the optimal control sequence. */
real_t *goL;  /**< Gradient vector over Lipschitz for the current system state. */
uint32_t *j_in;  /**< Maximun number of internal loop (FGM) iterations .*/
real_t *HoL;  /**< Hessian matrix of QP over Lipschitz constant. */
real_t *u_lb;  /**< Lower bound constraint of the inputs for condensed QP. */
real_t *u_ub;  /**< Upper bound constraint of the inputs for condensed QP. */
real_t *nu;  /**< Fast gradient extra step constant. */
real_t *tmp1_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
real_t *tmp2_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
real_t *tmp3_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
real_t *tmp4_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
real_t *tmp5_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
real_t *tmp6_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
uint32_t optvar_veclen;  /**< The length of each vector in the optimation variable sequence. */
uint32_t optvar_seqlen;  /**< The full length of optimization variable sequence. */
uint32_t sizeof_optvar_seqlen;  /**< Number of bytes in the optimization variable sequence. */
}};  /**< Variables used by the fast gradient method. */

/* External function declarations */

extern void {prefix}_fgm_solve_problem(const struct {prefix}_fgm *fgm);

extern void {prefix}_fgm_warm_start(const struct {prefix}_fgm *fgm);

extern void {prefix}_fgm_minimize_qp_iteration(const struct {prefix}_fgm *fgm,
                real_t u[], real_t u_old[], real_t w[], const real_t gradoL[]);

extern void {prefix}_fgm_compute_grad_over_L(const struct {prefix}_fgm *fgm,
                real_t gradoL[], const real_t w[]);

extern void {prefix}_compute_gxoL(struct {prefix}_fgm *fgm, const real_t x[]);

#ifdef __cplusplus
}}
#endif

#endif /* {PREFIX}_FGM_H */
