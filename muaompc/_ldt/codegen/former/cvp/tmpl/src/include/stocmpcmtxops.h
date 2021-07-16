#ifndef STOCMPCMTXOPS_H
#define STOCMPCMTXOPS_H

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"


/* matrix-matrix multiplication: pout = pmtxA * pmtxB.
 * pmtxA has size (rowsA x colsA),
 * pmtxB has size (colsA x colsB),
 * pout has size (rowsA x colsB) */
extern void stocmpc_mtx_multiply_mtx_mtx(real_t pout[], const real_t pmtxA[],
        const real_t pmtxB[],
        const uint32_t rowsA,
        const uint32_t colsA,
        const uint32_t colsB);

/* matrix-vector multiplication: pout = pmtx * pvec.
 * pmtx has size (rows x columns) */
extern void stocmpc_mtx_multiply_mtx_vec(real_t pout[], const real_t pmtx[],
        const real_t pvec[],
        const uint32_t rows,
        const uint32_t cols);

/* element-wise matrix factorisation: pout = pmtx * factor.
 * pmtx has size (rows x columns) */
extern void stocmpc_mtx_scale(real_t pout[], const real_t pmtx[],
        const real_t factor, const uint32_t rows,
        const uint32_t cols);

/* matrix addition: pmtxc = pmta + pmtxb.
 * pmta, pmtb, pmtc have size (rows x columns) */
extern void stocmpc_mtx_add(real_t pmtxc[], const real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols);

/* matrix substraction: pmtxc = pmta - pmtxb.
 * pmta, pmtb, pmtc have size (rows x columns) */
extern void stocmpc_mtx_substract(real_t pmtxc[], const real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols);

/* matrix multiplication and addition: pvecc = pvecc + pmtxa*pvecb.
 * pmtxa has shape (rows x columns); pvecc and pvecb are column vectors
 * pvecab = pmtxa * pvecb */
extern void stocmpc_mtx_mul_add(real_t pvecc[], real_t pvecab[], const real_t pmtxa[],
        const real_t pvecb[], const uint32_t rows,
        const uint32_t cols);

/* elementwise 2-sided vector saturation. 
 * lower[i] < vec[i] < upper[i]. pvec has size (rows x 1) */
extern void stocmpc_mtx_saturate_vec(real_t pvec[], const real_t plower[],
        const real_t pupper[],
        const uint32_t rows);

/* matrix transposition. mtxout[j,i] = mtxin[i,j].
 * mtxin has size (rows x columns) */
extern void stocmpc_mtx_transpose(real_t * mtxout, const real_t * mtxin,
        const uint32_t rows, const uint32_t cols);

/* elementwise vector maximum against zero
 * pmax[i] = max(pa[i], 0). pmax has size (rows x 1) */
extern void stocmpc_mtx_max_vec_zero(real_t pmax[],    const uint32_t rows);

/* elementwise vector minimum against zero
 * pmin[i] = min(pa[i], 0). pmin has size (rows x 1) */
extern void stocmpc_mtx_min_vec_zero(real_t pmin[], const uint32_t rows);

/* shift a sequence of seqlen=N*veclen elements, veclen elements to the left
 * outseq[i] = seq[i + veclen], the last veclen elements are set to zero */
extern void stocmpc_mtx_shift_sequence(real_t outseq[],
        const real_t seq[], const uint32_t veclen, const uint32_t seqlen);
#endif
