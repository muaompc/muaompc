#ifndef {PREFIX}MTXOPS_H
#define {PREFIX}MTXOPS_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"


/* matrix-matrix multiplication: pout = pmtxA * pmtxB.
 * pmtxA has size (rowsA x colsA),
 * pmtxB has size (colsA x colsB),
 * pout has size (rowsA x colsB) */
extern void {prefix}_mtx_multiply_mtx_mtx(real_t pout[], const real_t pmtxA[],
        const real_t pmtxB[],
        const uint32_t rowsA,
        const uint32_t colsA,
        const uint32_t colsB);

/* matrix-vector multiplication: pout = pmtx * pvec.
 * pmtx has size (rows x columns) */
extern void {prefix}_mtx_multiply_mtx_vec(real_t pout[], const real_t pmtx[],
        const real_t pvec[],
        const uint32_t rows,
        const uint32_t cols);

/* element-wise matrix factorisation: pout = pmtx * factor.
 * pmtx has size (rows x columns) */
extern void {prefix}_mtx_scale(real_t pout[], const real_t pmtx[],
        const real_t factor, const uint32_t rows,
        const uint32_t cols);

/* element-wise matrix factorisation: pout *= factor.
 * pout has size (rows x columns) */
extern void {prefix}_mtx_scale_direct(real_t pout[],
        const real_t factor, const uint32_t rows,
        const uint32_t cols);

/* matrix addition: pmtxc = pmta + pmtxb.
 * pmta, pmtb, pmtc have size (rows x columns) */
extern void {prefix}_mtx_add(real_t pmtxc[], const real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols);

/* matrix addition: pmtxa += pmtxb.
 * pmta, pmtb have size (rows x columns) */
extern void {prefix}_mtx_add_direct(real_t pmtxa[], const real_t pmtxb[],
                    const uint32_t rows, const uint32_t cols);

/* matrix substraction: pmtxc = pmta - pmtxb.
 * pmta, pmtb, pmtc have size (rows x columns) */
extern void {prefix}_mtx_substract(real_t pmtxc[], const real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols);

/* matrix substraction: pmta -= pmtxb.
 * pmta, pmtb have size (rows x columns) */
extern void {prefix}_mtx_substract_direct(real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols);

/* matrix multiplication and addition: pvecc = pvecc + pmtxa*pvecb.
 * pmtxa has shape (rows x columns); pvecc and pvecb are column vectors
 * pvecab = pmtxa * pvecb */
extern void {prefix}_mtx_mul_add(real_t pvecc[], real_t pvecab[], const real_t pmtxa[],
        const real_t pvecb[], const uint32_t rows,
        const uint32_t cols);

/* elementwise 2-sided vector saturation. 
 * lower[i] < vec[i] < upper[i]. pvec has size (rows x 1) */
extern void {prefix}_mtx_saturate_vec(real_t pvec[], const real_t plower[],
        const real_t pupper[],
        const uint32_t rows);

/* matrix transposition. mtxout[j,i] = mtxin[i,j].
 * mtxin has size (rows x columns) */
extern void {prefix}_mtx_transpose(real_t * mtxout, const real_t * mtxin,
        const uint32_t rows, const uint32_t cols);

/* elementwise vector maximum against zero
 * pmax[i] = max(pa[i], 0). pmax has size (rows x 1) */
extern void {prefix}_mtx_max_vec_zero(real_t pmax[],    const uint32_t rows);

/* elementwise vector minimum against zero
 * pmin[i] = min(pa[i], 0). pmin has size (rows x 1) */
extern void {prefix}_mtx_min_vec_zero(real_t pmin[], const uint32_t rows);

/* shift a sequence of seqlen=N*veclen elements, veclen elements to the left
 * outseq[i] = seq[i + veclen], the last veclen elements are set to zero */
extern void {prefix}_mtx_shift_sequence(real_t outseq[],
        const real_t seq[], const uint32_t veclen, const uint32_t seqlen);

/* identity matrix for mtx with size (dimension x dimension) */
extern void {prefix}_mtx_to_eye(real_t mtx[], const uint32_t dimension);

/* zero matrix for mtx with nb_e entries */
extern void {prefix}_mtx_to_zero(real_t mtx[], const uint32_t nb_e);

/* cholesky factorization: chol_fac*chol_fac^T = mtx
 * mtx have size (dimension x dimension),
 * mtx is poitiv definite and symmetric,
 * chol_fac is a lower triangular matrix */
extern void {prefix}_mtx_cholesky(real_t solution[],
              const real_t mtx[], const uint32_t dimension);

/* check if mtx is positiv definite and symmetric.
 * return 1 if condition is fullfilled, 0 otherwise
 * mtx have size (dimension x dimension) */
extern uint32_t {prefix}_mtx_check_chol(real_t chol_fac[],
                    const real_t mtx[], const uint32_t dimension);

/* forward substitution: mtxA * solution = mtxB
 * mtxA has shape (dimA x dimA), lower Dreiecksmatrix
 * solution and mtxB have shape (dimA x columsB)*/
extern void {prefix}_mtx_fwd_subst(real_t solution[],
               const real_t mtxA[], const uint32_t dimensionA,
               const real_t mtxB[], const uint32_t columsB);

/* backward substitution: mtxA * solution = mtxB
 * mtxA has shape (dimA x dimA), upper Dreiecksmatrix
 * solution and mtxB have shape (dimA x columsB)*/
extern void {prefix}_mtx_bwd_subst(real_t solution[],
               const real_t mtxA[], const uint32_t dimensionA,
               const real_t mtxB[], const uint32_t columsB);

/* elementwise comparision between mtxA and mtxB.
 * return -1 if |mtxA[i]-mtxB[i]| <= accuracy,
 * return last i for which the above condition is false otherwise */
extern uint32_t {prefix}_mtx_mtx_cmp(const real_t mtxA[], const real_t mtxB[], real_t dim,
                 real_t accuracy);
                
extern real_t {prefix}_mtx_smpl_exp(real_t exp);

extern real_t {prefix}_mtx_smpl_abs(real_t x);

#ifdef __cplusplus
}}
#endif

#endif
