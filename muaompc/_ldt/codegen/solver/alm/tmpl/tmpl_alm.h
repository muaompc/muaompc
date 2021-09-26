#ifndef {PREFIX}_ALM_H
#define {PREFIX}_ALM_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}fgm.h"

/* It is problem specific and defined by the automatic code generation. */

struct {prefix}_alm_conf {{
    uint32_t in_iter;  /**< Maximum number of internal loop (FGM) iterations (j_in). */
    uint32_t ex_iter;  /**< Maximum number of external loop (ALM) iterations (i_ex). */
    uint32_t warm_start;  /**< If not 0, automatically warm start next algorithm iteration. */
}}; /**< Configuration parameters of the {PREFIX} algorithm.
 * Some of them are pointed at by the respective algorithm internal variable. */


struct {prefix}_alm {{
    struct {prefix}_alm_conf *conf;  /**< Algorithm configuration data. */
    struct {prefix}_fgm *fgm;  /**< Fast gradient method data */
    real_t *u_opt;  /**< Solution to the optimal control problem. */
    real_t *u_ini;  /**< Initial guess for the optimal control sequence. */
    real_t *l_opt;  /**< Optimal Lagrange multiplier. */
    real_t *l_ini;  /**< Initial guess for the optimal multiplier. */
    real_t *V;  /**< Linear factor of 2-sided mixed constraint. */
    real_t *v_lb;  /**< Mixed constraint lower bound as function of current state. */
    real_t *v_ub;  /**< Mixed constraint upper bound as function of current state. */
    uint32_t *i_ex;  /**< Maximum number of external loop (ALM) iterations. */
    real_t *mu;  /**< Augmented Lagrange method penalty parameter. */
    real_t *Linv;  /**< Inverse of gradient Lipschitz constant (1/L). */

    real_t *tmp1_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp2_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp3_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp4_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp5_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp6_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp7_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp1_lagmul_seqlen;  /**< Temporary variable of length lagmul_seqlen. */
    real_t *tmp2_lagmul_seqlen;  /**< Temporary variable of length lagmul_seqlen. */
    real_t *tmp3_lagmul_seqlen;  /**< Temporary variable of length lagmul_seqlen. */
    real_t *tmp4_lagmul_seqlen;  /**< Temporary variable of length lagmul_seqlen. */
    real_t *tmp5_lagmul_seqlen;  /**< Temporary variable of length lagmul_seqlen. */
    real_t *tmp6_lagmul_seqlen;  /**< Temporary variable of length lagmul_seqlen. */
    real_t *tmp7_lagmul_seqlen;  /**< Temporary variable of length lagmul_seqlen. */
    real_t *tmp1_lagmul_optvar_seqlen;  /**< Temporary variable of length lagmul_seqlen * optvar_seqlen. */
    uint32_t optvar_veclen;  /**< The length of each vector in the optimation variable sequence. */
    uint32_t optvar_seqlen;  /**< The full length of optimization variable sequence. */
    uint32_t sizeof_optvar_seqlen;  /**< Number of bytes in the optimization variable sequence. */
    uint32_t lagmul_veclen;  /**< The length of each vector in the Lagrange multiplier sequence. */
    uint32_t lagmul_seqlen;  /**< The full length of Lagrange multiplier sequence. */
    uint32_t sizeof_lagmul_seqlen;  /**< Number of bytes in the Lagrange multiplier sequence. */
}};  /**< Variables used by the augmented Lagrangian method. */

/* External function declarations */

extern void {prefix}_alm_solve_problem(struct {prefix}_alm *alm);

#ifdef __cplusplus
}}
#endif

#endif /* {PREFIX}_ALM_H */
