#include "include/stocmpcmtxops.h"
#ifdef STOCMPCFIPOPS
#include "stocmpcfipops.h"
#endif


/* matrix-matrix multiplication: pout = pmtxA * pmtxB.
 * pmtxA has size (rowsA x colsA),
 * pmtxB has size (colsA x colsB),
 * pout has size (rowsA x colsB) */
extern void stocmpc_mtx_multiply_mtx_mtx(real_t pout[], const real_t pmtxA[],
        const real_t pmtxB[],
        const uint32_t rowsA,
        const uint32_t colsA,
            const uint32_t colsB)
{
    uint32_t i, j, k; /* loop counters */

    for (i = 0; i < rowsA; i++) {
        for (j = 0; j < colsB; j++) {
                pout[i*colsB+j] = 0.;
            for (k = 0; k < colsA; k++) {
                    pout[i*colsB+j] += pmtxA[i*colsA + k] * pmtxB[k*colsB + j];
            }
        }
    }
    return;
}


/* matrix-vector multiplication: pout = pmtx * pvec.
 * pmtx has size (rows x columns) */
void stocmpc_mtx_multiply_mtx_vec(real_t pout[], const real_t pmtx[],
        const real_t pvec[],
        const uint32_t rows,
        const uint32_t cols)
{
    uint32_t i; /* row number */
    uint32_t j; /* column number */
    uint32_t k = 0; /* matrix index (row * column) */
    for (i = 0; i < rows; i++) {
        pout[i] = 0;
        for (j = 0; j < cols; j++) {
#ifndef STOCMPCFIP_OPS
            pout[i] += pmtx[k] * pvec[j];
#else
            pout[i] += fip_mul(pmtx[k], pvec[j]);
#endif
            k++;
        }
    }
    return;
}

/* element-wise matrix factorisation: pout = pmtx * factor.
 * pmtx has size (rows x columns) */
void stocmpc_mtx_scale(real_t pout[], const real_t pmtx[],
        const real_t factor, const uint32_t rows,
        const uint32_t cols)
{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {
#ifndef STOCMPCFIP_OPS
        pout[k] = pmtx[k] * factor;
#else
        pout[k] = fip_mul(pmtx[k], factor);
#endif
    }
    return;
}

/* matrix addition: pmtxc = pmta + pmtxb.
 * pmta, pmtb, pmtc have size (rows x columns) */
void stocmpc_mtx_add(real_t pmtxc[], const real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols)
{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {
        pmtxc[k] = pmtxa[k] + pmtxb[k];
    }
    return;
}

/* matrix substraction: pmtxc = pmta - pmtxb.
 * pmta, pmtb, pmtc have size (rows x columns) */
void stocmpc_mtx_substract(real_t pmtxc[], const real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols)
{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {
        pmtxc[k] = pmtxa[k] - pmtxb[k];
    }
    return;
}

/* matrix multiplication and addition: pvecc = pvecc + pmtxa*pvecb.
 * pmtxa has shape (rows x columns); pvecc and pvecb are column vectors */
void stocmpc_mtx_mul_add(real_t pvecc[], real_t pvecab[], const real_t pmtxa[],
        const real_t pvecb[], const uint32_t rows,
        const uint32_t cols)
{
    uint32_t k; /* vector index (rows) */
  stocmpc_mtx_multiply_mtx_vec(pvecab, pmtxa, pvecb, rows, cols);
    for (k = 0; k < rows; k++) {
        pvecc[k] += pvecab[k];
    }
    return;
}


/* elementwise 2-sided vector saturation.
 * lower[i] < vec[i] < upper[i]. pvec has size (rows x 1) */
void stocmpc_mtx_saturate_vec(real_t pvec[], const real_t plower[],
        const real_t pupper[],
        const uint32_t rows)
{
    uint32_t i; /* vector index (row number) */

    for (i = 0; i < rows; i++) {
        if (pvec[i] > pupper[i]) {
            pvec[i] = pupper[i];
        } else if (pvec[i] < plower[i]) {
            pvec[i] = plower[i];
        }
    }
}

/* elementwise vector maximum against zero
 * pmax[i] = max(pa[i], 0). pmax has size (rows x 1) */
void stocmpc_mtx_max_vec_zero(real_t pmax[], const uint32_t rows)
{
    uint32_t i; /* vector index (row number) */
    real_t zero = 0.;
    for (i = 0; i < rows; i++) {
        if (pmax[i] < zero) {
            pmax[i] = zero;
        }
    }
}

/* elementwise vector minimum against zero
 *  pmin[i] = min(pa[i], 0). pmin has size (rows x 1) */
void stocmpc_mtx_min_vec_zero(real_t pmin[], const uint32_t rows)
{
    uint32_t i; /* vector index (row number) */
    real_t zero = 0.;
    for (i = 0; i < rows; i++) {
        if (pmin[i] > zero) {
            pmin[i] = zero;
        }
    }
}


/* matrix transposition. mtxout[j,i] = mtxin[i,j].
 * mtxin has size (rows x columns) */
void stocmpc_mtx_transpose(real_t * mtxout, const real_t * mtxin,
        const uint32_t rows, const uint32_t cols)
{
    uint32_t i; /* row number */
    uint32_t j; /* column number */
    uint32_t k = 0; /* matrix index (row * column) */
    uint32_t kT; /* matrix transpose index */

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            kT = j * rows + i;
            mtxout[kT] = mtxin[k];
            k++;
        }
    }
    return;
}

/* shift a sequence of seqlen=N*veclen elements, veclen elements to the left
 * outseq[i] = seq[i + veclen], the last veclen elements are set to zero */
void stocmpc_mtx_shift_sequence(real_t outseq[], const real_t seq[],
        const uint32_t veclen, const uint32_t seqlen)
{
    uint32_t i;
    /* shift vector one horizon step backwards */
    for (i=0; i < (seqlen - veclen); i++) {
        outseq[i] = seq[i + veclen];
    }

    /* set last element to zero */
    for (i=(seqlen - veclen); i < seqlen; i++) {
        outseq[i] = 0.;
    }
    return;
}
