#ifndef {PREFIX}_PBM_H
#define {PREFIX}_PBM_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "mc04types.h"
#include "arithmetic.h"

{pbm_sparse}{pbm_cond}

/* Second order cone program solver based on a primal barrier interior point
 * method.
 */ 

struct {prefix}_pbm_conf {{
    uint32_t in_iter;  /**< Maximum number of internal loop (PBM) iterations
                            (j_in). */
    uint32_t warm_start;  /**< If not 0, automatically warm start next
                               algorithm iteration. */
    real_t reg;  /**< Value of the regularization factor. */
}};  /**< Configuration parameters of the {prefix} algorithm. */

struct {prefix}_pbm_block {{  /* TODO rename block */
    uint32_t rows;
    uint32_t cols;
    real_t *data;
}};  /**< One block of a Matrix. */

struct {prefix}_pbm_qc {{
    real_t *Gamma;  /**< size (dimGamma x dimGamma) */
    uint32_t dimGamma;  /**< Dimension of square matrix Gamma */
    real_t *beta;  /**< size (dimGamma x 1) */
    real_t *alpha;  /**< size (1 x 1) */
    real_t *par;  /**< Pointer to first entry of associated part in z */
    uint32_t par_0;  /**< First entry of associated part in optim. vector z */
    /* par_l = dimGamma */
    uint32_t par_l;  /**< Length of associated part in optim. vector z */
}};  /**< Parameters for a Quadratic Constraint:
          (z')^T*Gamma*z' <= beta^T*z' + alpha
          associated part in z: z' = z[par_0 : par_0+par_l] */

struct {prefix}_pbm_socc {{
    real_t *A;  /**< size (rowsA x colsA) */
    real_t *A_T;  /**< Transposed of A: size (colsA x rowsA) */
    uint32_t rowsA;  /**< Number of rows of matrix A. */
    uint32_t colsA;  /**< Number of colums of matrix A. */
    real_t *b;  /**< size (rowsA x 1) */
    real_t *c;  /**< size (colsA x 1) */
    real_t *d;  /**< size (1 x 1) */
    real_t *AAmcc;  /**< Term: A^T*A - c*c^T */
    real_t *par;  /**< Pointer to first entry of associated part in z */
    uint32_t par_0;  /**< First entry of associated part in optim. vector z */
    /* par_l = colsA */
    uint32_t par_l;  /**< Length of associated part in optim. vector z */
}};  /**< Parameters for a Second Order Cone Constraint:
          |A*z' + b| <= c*z' + d
          associated part in z: z' = z[par_0 : par_0+par_l]*/

struct {prefix}_pbm_P_hat {{
    real_t factor;  /* for factor 2 in gradient P */
    real_t *P; /**< Unchanged matrix P for linear inequalities. */
    real_t *P_hat; /**< Extended matrix P with additionell lines for qcs and
                        soccs. */
    real_t *P2_hat; /**< P_hat calculated with argument 2*z. */
    real_t *P_hat_T; /**< Transposed of P_hat. */
    real_t *P2_hat_T; /**< Transposed of P2_hat. */
    real_t *h; /**< Unchanged vector h for linear inequalities. */
    real_t *h_hat; /**< Extended vector h with additionell lines for qcs and
                        soccs. */
    
    struct {prefix}_pbm_qc **qc; /**< Data for qcs. */
    struct {prefix}_pbm_socc **socc; /**< Data for soccs. */
    uint32_t nb_lin_constr; /**< Number of linear inequality constraints. */
    uint32_t nb_qc; /**< Number of quadratic inequality constraints. */
    uint32_t nb_socc; /**< Number of conic inequality constraints. */
}};  /**< Variables for updating matrix P_hat (depending on z). */

/* TODO Clean up variables. */
/* TODO Clean up temporary variables. */
struct {prefix}_pbm {{
    real_t *u_opt;
    struct {prefix}_pbm_conf *conf;  /* Algorithm configuration data. */
    real_t *z_opt;  /**< Solution to the optimal control problem. */
    real_t *z_ini;  /**< Initial guess for the optimal control sequence. */
    real_t *zref;  /**< Referenz value for the optimal control sequence. */
    real_t *delta_z;  /**< Search direction for the optimal control sequence. */
    real_t *v_ini;  /**< Initial guess for the dual variable of the optimal control problem. */
    real_t *v_opt;  /**< Solution to the dual variable of the optimal control problem. */
    real_t *delta_v;  /**< Search direction for the dual variable of the optimal control problem. */
    real_t *q;  /**< Vector for linear weighting of the states. */
    real_t *r;  /**< Vector for linear weighting of the inputs. */
    real_t *g;  /**< Combined vector for q und r for all time steps. */
    real_t *A;  /**< System matrix. */
    real_t *A_T;  /**< Transposed of matrix A.*/
    real_t *B;  /**< Input matrix. */
    real_t *B_T;  /**< Transposed of matrix B. */
    real_t *A_B;  /**< Term: combined matrix [A, B]. */
    real_t *A_B_T;  /**< Transposed of combined matrix [A, B]. */
    real_t *C;  /**< Combined matrix for A und B for all time steps. */
    real_t *C_T;  /**< Transposed of matrix C. */
    real_t *H;  /**< Combined Matrix for weighting matrizes for all time steps. */
    struct {prefix}_pbm_P_hat *P_of_z;  /**< Structure for updating matrix P_hat (depending on z). */
    real_t *P;  /**< Pointer to array in P_of_z (P_hat). */
    real_t *P2;  /**< Pointer to array in P_of_z (P2_hat). */
    real_t *P_T;  /**< Pointer to array in P_of_z (P_hat_T). */
    real_t *P2_T;  /**< Pointer to array in P_of_z (P2_hat_T). */
    
    real_t *Fusoft;  /**< Matrizes for soft constraints*/
    real_t *Fxsoft;
    real_t *Ffsoft;
    real_t *Psoft;
    real_t *Psoft_T;
    struct {prefix}_pbm_block **Phi_sft_blks;  /* TODO Calculate and save only blocks in Phi and Phi_soft. */
    struct {prefix}_pbm_block *tmp_Phi_sft_blk;  /**< Temporary variable of size (optvar_seqlen x nb_soft_constr) */ 
    
    real_t *b;  /**< Combined vector for disturbances. */
    real_t *h;  /**< Combined vector for inequality constraints for all time steps. */
    real_t *hsoft;  /**< Combined vector for soft constraints for all time steps. */

    real_t *u_lb;  /**< Lower bound of input box constraints from condensed problem */
    real_t *u_ub;  /**< Upper bound of input box constraints from condensed problem */
    real_t *v_lb;  /**< Lower bound of mixed box constraints from condensed problem */
    real_t *v_ub;  /**< Upper bound of mixed box constraints from condensed problem */
    
    real_t *d;  /**< Vector of length nb_of_ueq_constr. */
    real_t *dsoft;  /**< Additionell term for vector d for soft constraints. */
    /* TODO Use vector of length nb_of_ueq_constr for diagonal matrix diag_d_sq. */
    real_t *diag_d_sq;  /**< Matrix of size (nb_of_ueq_constr x nb_of_ueq_constr) */
    real_t *diag_d_soft;  /**< Additionell term for matrix diag_d_sq for soft constraints. */
    real_t *Phi;  /**< Matrix of size (optvar_seqlen x optvar_seqlen) */
    real_t *Phi_soft;  /**< Additionell term for matrix Phi for soft constraints. */
    real_t *r_p;  /**< Primal residual: size (dual_seqlen). */
    real_t *r_d;  /**< Dual residual. size (optvar_seqlen). */
    real_t *r_d_soft;  /**< Additionell term in primal residual for soft constraints. */
    
    real_t *st_size;  /**< Variable for step size in backtracking line search. */
    real_t *reg;  /**< Value of the regularization factor. */
    uint32_t *j_in;  /**< Maximum number of internal loop (PBM) iterations. */
    
    real_t *kappa;  /**< Barrier parameter. */
    real_t *roh;  /**< Parameter for soft constraint vioalation penalty */
    
    uint32_t horizon; /**< Prediction horizon. */
    uint32_t optvar_veclen;  /**< The length of each vector in the optimation variable sequence. */
    uint32_t optvar_seqlen;  /**< The full length of optimization variable sequence. */
    uint32_t state_veclen;  /**< Dimension of state variable x_k. */
    uint32_t dual_seqlen;  /**< Full length of dual variable v associated with the equalitiy constraints. */
    uint32_t control_veclen;  /**< Dimension of control variable u_k. */
    uint32_t nb_of_ueq_constr;  /**< Number of inequality constraints. */
    uint32_t nb_of_soft_constr;  /**< Number of soft constraints. */
    uint32_t rowsFusoft;  /**< Number of soft constraints in one time step. */
    uint32_t rowsFfsoft;  /**< Number of soft constraints in the last time step. */
    uint32_t sizeof_optvar_veclen;  /**< Number of bytes in the optimization variable vector. */
    uint32_t sizeof_dual_seqlen;  /**< Number of bytes in the optimization dual variable sequence. */
    uint32_t sizeof_optvar_seqlen;  /**< Number of bytes in the optimization variable sequence. */
    
    real_t *tmp1_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp2_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp6_optvar_seqlen;  /**< Temporary variable of length optvar_seqlen. */
    real_t *tmp2_dual_seqlen;  /**< Temporary variable of length dual_seqlen. */
    real_t *tmp7_dual_seqlen;  /**< Temporary variable of length dual_seqlen. */
    real_t *tmp3_state_veclen;  /**< Temporary variable of length state_veclen. */
    real_t *tmp3_mtx_optvar_nb_of_ueq;  /**< Temporary variable of length (optvar_seqlen x nb_of_ueq_constr). */
    real_t *tmp3_mtx_optvar_nb_of_soft;  /**< Temporary variable of length (optvar_seqlen x nb_of_soft_constr). */
    real_t *tmp4_mtx_optvar_optvar;  /**< Temporary variable of length (optvar_seqlen x optvar_seqlen). */
    real_t *tmp5_nb_of_constr;  /**< Temporary variable of length nb_of_ueq_constr. */
    real_t *tmp8_L_Y;  /**< Temporary variable of length (2*horizon - 1)*state_veclen*state_veclen. */
    real_t *tmp9_L_Y_T;  /**< Temporary variable of length (2*horizon - 1)*state_veclen*state_veclen. */
    real_t *tmp8_L_Phi;  /**< Temporary variable of length horizon*optvar_veclen*optvar_veclen. */
    real_t *tmp9_L_Phi_T;  /**< Temporary variable of length horizon*optvar_veclen*optvar_veclen. */
    real_t *tmp_phibl1;  /**< Temp var for matrix of size optvar_veclen x optvar_veclen */
    real_t *tmp_phibl2;  /**< Temp var for matrix of size optvar_veclen x optvar_veclen */
    real_t *tmp_phibl3;  /**< Temp var for matrix of size optvar_veclen x optvar_veclen */
    real_t *tmp10;  /**< Temp var for vectors up to size (optvar_veclen x optvar_veclen) */
    real_t *tmpYbl;  /**< Temp var for matrix of size (state_veclen x state_veclen) */
    real_t *tmpQbl;  /**< Temp var for matrix of size (state_veclen x state_veclen) */
    real_t *tmp1_res_os;  /**< Temp var for residual calculations, optvar_seqlen x 1 */
    real_t *tmp2_res_os;  /**< Temp var for residual calculations, optvar_seqlen x 1 */
    real_t *tmp3_res_ds;  /**< Temp var for residual calculations, dual_seqlen x 1 */
    real_t *eye_optvar_veclen;  /**< Identity matrix of size (optvar_veclen x optvar_veclen) */
    real_t *eye_state_veclen;  /**< Identity matrix of size (state_veclen x state_veclen) */
}};  /**< Variables used by the primal barrier interior point method. */

/* External function declarations */

extern void {prefix}_pbm_solve_problem(const struct {prefix}_pbm *pbm);

extern void {prefix}_pbm_calc_kappa(
                real_t *kappa, const struct {prefix}_pbm *pbm,
                const real_t *z);

extern void {prefix}_pbm_test_get_valid(const struct {prefix}_pbm *pbm);

extern void {prefix}_pbm_iterative_refinement(const struct {prefix}_pbm *pbm);

/* Update matrix P(z) and P(2*z) for lines of qc and socc */
extern void {prefix}_pbm_update_P(
                struct {prefix}_pbm_P_hat *P, const uint32_t optvar_seqlen,
                real_t *tmp1_optvar_seqlen, real_t *tmp2_optvar_seqlen);

extern void {prefix}_form_pbm_matrices(real_t *d, real_t *diag_d_sq,
                                       real_t *dsoft, real_t *diag_d_soft,
                                       real_t *Phi,
                                       const struct {prefix}_pbm *pbm);

extern uint32_t {prefix}_pbm_check_valid(
                const struct {prefix}_pbm *pbm, const real_t *z_check);

extern void {prefix}_pbm_get_positiv(
                const struct {prefix}_pbm *pbm, real_t *delta_to_zero);

extern void {prefix}_pbm_get_valid_lin_constr(
                const struct {prefix}_pbm *pbm, real_t *delta_to_zero);

extern void {prefix}_pbm_get_valid_trick(const struct {prefix}_pbm *pbm);

extern uint32_t {prefix}_pbm_check_positiv(
                const struct {prefix}_pbm *pbm, const real_t *z_check);

extern void residual(
                const struct {prefix}_pbm *pbm,
                const real_t *z, const real_t *v, const real_t *d,
                const real_t kappa);

#ifdef __cplusplus
}}
#endif

#endif  /* {PREFIX}_PBM_H */