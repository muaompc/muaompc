#ifndef {PREFIX}_PBMSOLVE_H
#define {PREFIX}_PBMSOLVE_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "arithmetic.h"
#include "{prefix}mtxops.h"

{pbm_sparse}{pbm_cond}
/* #include <{prefix}usefull.h> */  /* TODO implement {prefix}usefull.h */

/* Solve systems of linear equation arising in used primal barrier interior
 * point method.
 */ 

/* External function declarations */

extern void {prefix}_multiply_C_z(
                real_t *product, const real_t *C, const real_t *z,
                const uint32_t n, const uint32_t m, const uint32_t T);

extern void {prefix}_multiply_C_T_v(
                real_t *product, const real_t *C_T, const real_t *v,
                const uint32_t n, const uint32_t m, const uint32_t T);

void solve_sysofleq(real_t delta_z[], real_t delta_v[],
                    const real_t Phi[],
                    const real_t rd[], const real_t rp[],
                    const real_t C[],
                    const real_t *C_T,
                    const real_t A[],
                    const real_t A_T[],
                    const real_t B[],
                    const real_t B_T[],
                    const real_t A_B[],
                    const real_t A_B_T[],
                    const real_t H[],
                    const uint32_t dimA, const uint32_t colsB, const uint32_t horizon,
                    const real_t *eye_nm, const real_t *eye_n,
                    real_t *tmp_1_optvar_seqlen,
                    real_t *tmp_2_optvar_seqlen,
                    real_t *tmp_dual_seqlen,
                    real_t *L_Y, real_t *L_Y_T,
                    real_t *PhiBlock, real_t *PhiBlock_I,
                    real_t *PhiBlock_I_last,
                    real_t *Block_nxn1, real_t *Block_nxn2,
                    real_t *tmp_optvar_veclenxoptvar_veclen,
                    real_t *L_Phi_blocks, real_t *L_Phi_T_blocks);

void form_delta_z(real_t delta_z[],
                  real_t *tmp_optvar_seqlen,
                  const real_t delta_v[],
                  const real_t mtxL_Phi[],
                  const real_t mtxL_Phi_T[],
                  const real_t rd[],
                  const real_t C[],
                  const uint32_t T, const uint32_t n, const uint32_t m);

void form_delta_v(real_t delta_v[],
                  real_t *tmp_dual_seqlen,
                  real_t *tmp_state_veclen,
                  const real_t L_Y[], const real_t L_Y_T[],
                  const uint32_t T, const uint32_t n);

void form_beta(real_t beta[],
               const real_t mtxL_Phi[],
               const real_t mtxL_Phi_T[],
               const real_t rd[], const real_t rp[],
               const uint32_t horizon,
               const real_t mtxA[], const uint32_t dimA,
               /*const real_t mtxB[],*/ const uint32_t colsB,
               real_t *tmp1_optvar_seqlen, real_t *tmp2_optvar_seqlen);

/*returns also mtxL_Phi the cholesky factorization of Phi
 */
void form_Y(real_t mtxL_Y_blocks[], real_t *mtxL_Y_T_blocks,
            real_t mtxL_Phi_blocks[], real_t mtxL_Phi_T_blocks[],
            real_t mtxPhi[],
            const uint32_t horizon,
            const real_t mtxA[], const real_t *A_B, const uint32_t dimA,
            const real_t mtxB[], const real_t *A_T_B_T, const uint32_t colsB,
            const real_t *eye_nm, const real_t *eye_n,
            real_t *tmp_PhiBlock, real_t *PhiBlock_I, real_t *PhiBlock_I_last,
            real_t *bl_nxn1, real_t *bl_nxn2, real_t *tmp);

/* Solves system of linear equation of type mtxA * solution = mtxB
 * mtxA has size dimA x dimA, mtxB and solution have size dimA x colsB
 * additionally returns cholesky factorization of mtxA  */
void solveBlock(real_t *mtxA_Invers, real_t *L_mtxA, real_t *L_mtxA_T,
                real_t *mtxA, const uint32_t dimA, const real_t *mtxB, const uint32_t colsB,
                real_t *tmp_size_dimAxcolsB);

void form_Yii(real_t solution[],
              const real_t A[], const uint32_t rowsA, const uint32_t colsA,
              const real_t B[], const uint32_t rowsB, const uint32_t colsB,
              const real_t C[], const uint32_t rowsC, const uint32_t colsC,
              const real_t Q[],
              real_t *tmp_nxoptvar_veclen);

void form_Y_i_ip1(real_t solution[],
                  const real_t A_T_B_T[],
                  const uint32_t rowsATBT, const uint32_t colsATBT,
                  const real_t Qi_C[]);

void form_Y11(real_t solution[],
              const real_t B[], const real_t *B_T,
              const uint32_t rowsB, const uint32_t colsB,
              const real_t R0_Cholesky[], const real_t R0_Cholesky_T[],
              const real_t Q1[],
              real_t *tmp1_mxn, real_t *tmp2_mxn);

void setBlock(real_t mtx[], const uint32_t dimension, 
              const real_t block[], const uint32_t size_row, const uint32_t size_col,
              const uint32_t first_row, const uint32_t first_col);

void getBlock(real_t block[], 
              const real_t mtx[], const uint32_t dimension,
              const uint32_t first_row, const uint32_t first_col,
              const uint32_t size_row, const uint32_t size_col);

#ifdef __cplusplus
}}
#endif

#endif  /* {PREFIX}_PBMSOLVE_H */