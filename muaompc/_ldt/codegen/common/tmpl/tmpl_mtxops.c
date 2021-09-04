#include "{prefix}mtxops.h"
#ifdef {PREFIX}FIPOPS
#include "{prefix}fipops.h"
#endif


/* matrix-matrix multiplication: pout = pmtxA * pmtxB.
 * pmtxA has size (rowsA x colsA),
 * pmtxB has size (colsA x colsB),
 * pout has size (rowsA x colsB) */
extern void {prefix}_mtx_multiply_mtx_mtx(real_t pout[], const real_t pmtxA[],
        const real_t pmtxB[],
        const uint32_t rowsA,
        const uint32_t colsA,
            const uint32_t colsB)
{{
    uint32_t i, j, k; /* loop counters */

    for (i = 0; i < rowsA; i++) {{
        for (j = 0; j < colsB; j++) {{
                pout[i*colsB+j] = 0.;
            for (k = 0; k < colsA; k++) {{
                    pout[i*colsB+j] += pmtxA[i*colsA + k] * pmtxB[k*colsB + j];
            }}
        }}
    }}
    return;
}}


/* matrix-vector multiplication: pout = pmtx * pvec.
 * pmtx has size (rows x columns) */
void {prefix}_mtx_multiply_mtx_vec(real_t pout[], const real_t pmtx[],
        const real_t pvec[],
        const uint32_t rows,
        const uint32_t cols)
{{
    uint32_t i; /* row number */
    uint32_t j; /* column number */
    uint32_t k = 0; /* matrix index (row * column) */
    for (i = 0; i < rows; i++) {{
        pout[i] = 0;
        for (j = 0; j < cols; j++) {{
#ifndef {PREFIX}FIP_OPS
            pout[i] += pmtx[k] * pvec[j];
#else
            pout[i] += fip_mul(pmtx[k], pvec[j]);
#endif
            k++;
        }}
    }}
    return;
}}

/* element-wise matrix factorisation: pout = pmtx * factor.
 * pmtx has size (rows x columns) */
void {prefix}_mtx_scale(real_t pout[], const real_t pmtx[],
        const real_t factor, const uint32_t rows,
        const uint32_t cols)
{{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {{
#ifndef {PREFIX}FIP_OPS
        pout[k] = pmtx[k] * factor;
#else
        pout[k] = fip_mul(pmtx[k], factor);
#endif
    }}
    return;
}}

/* element-wise matrix factorisation: pout *= factor.
 * pout has size (rows x columns) */
void {prefix}_mtx_scale_direct(real_t pout[],
        const real_t factor, const uint32_t rows,
        const uint32_t cols)
{{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {{
        pout[k] *= factor;
    }}
    return;
}}

/* matrix addition: pmtxc = pmta + pmtxb.
 * pmta, pmtb, pmtc have size (rows x columns) */
void {prefix}_mtx_add(real_t pmtxc[], const real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols)
{{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {{
        pmtxc[k] = pmtxa[k] + pmtxb[k];
    }}
    return;
}}

/* matrix addition: pmtxa += pmtxb.
 * pmta, pmtb have size (rows x columns) */
void {prefix}_mtx_add_direct(real_t pmtxa[], const real_t pmtxb[],
                    const uint32_t rows, const uint32_t cols)
{{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {{
        pmtxa[k] += pmtxb[k];
    }}
    return;
}}

/* matrix substraction: pmtxc = pmta - pmtxb.
 * pmta, pmtb, pmtc have size (rows x columns) */
void {prefix}_mtx_substract(real_t pmtxc[], const real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols)
{{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {{
        pmtxc[k] = pmtxa[k] - pmtxb[k];
    }}
    return;
}}

/* matrix substraction: pmtxa -= pmtxb.
 * pmta, pmtb have size (rows x columns) */
void {prefix}_mtx_substract_direct(real_t pmtxa[],
        const real_t pmtxb[], const uint32_t rows,
        const uint32_t cols)
{{
    uint32_t k; /* matrix index (row * column) */
    for (k = 0; k < rows * cols; k++) {{
        pmtxa[k] -= pmtxb[k];
    }}
    return;
}}

/* matrix multiplication and addition: pvecc = pvecc + pmtxa*pvecb.
 * pmtxa has shape (rows x columns); pvecc and pvecb are column vectors */
void {prefix}_mtx_mul_add(real_t pvecc[], real_t pvecab[], const real_t pmtxa[],
        const real_t pvecb[], const uint32_t rows,
        const uint32_t cols)
{{
    uint32_t k; /* vector index (rows) */
  {prefix}_mtx_multiply_mtx_vec(pvecab, pmtxa, pvecb, rows, cols);
    for (k = 0; k < rows; k++) {{
        pvecc[k] += pvecab[k];
    }}
    return;
}}


/* elementwise 2-sided vector saturation.
 * lower[i] < vec[i] < upper[i]. pvec has size (rows x 1) */
void {prefix}_mtx_saturate_vec(real_t pvec[], const real_t plower[],
        const real_t pupper[],
        const uint32_t rows)
{{
    uint32_t i; /* vector index (row number) */

    for (i = 0; i < rows; i++) {{
        if (pvec[i] > pupper[i]) {{
            pvec[i] = pupper[i];
        }} else if (pvec[i] < plower[i]) {{
            pvec[i] = plower[i];
        }}
    }}
}}

/* elementwise vector maximum against zero
 * pmax[i] = max(pa[i], 0). pmax has size (rows x 1) */
void {prefix}_mtx_max_vec_zero(real_t pmax[], const uint32_t rows)
{{
    uint32_t i; /* vector index (row number) */
    real_t zero = 0.;
    for (i = 0; i < rows; i++) {{
        if (pmax[i] < zero) {{
            pmax[i] = zero;
        }}
    }}
}}

/* elementwise vector minimum against zero
 *  pmin[i] = min(pa[i], 0). pmin has size (rows x 1) */
void {prefix}_mtx_min_vec_zero(real_t pmin[], const uint32_t rows)
{{
    uint32_t i; /* vector index (row number) */
    real_t zero = 0.;
    for (i = 0; i < rows; i++) {{
        if (pmin[i] > zero) {{
            pmin[i] = zero;
        }}
    }}
}}


/* matrix transposition. mtxout[j,i] = mtxin[i,j].
 * mtxin has size (rows x columns) */
void {prefix}_mtx_transpose(real_t * mtxout, const real_t * mtxin,
        const uint32_t rows, const uint32_t cols)
{{
    uint32_t i; /* row number */
    uint32_t j; /* column number */
    uint32_t k = 0; /* matrix index (row * column) */
    uint32_t kT; /* matrix transpose index */

    for (i = 0; i < rows; i++) {{
        for (j = 0; j < cols; j++) {{
            kT = j * rows + i;
            mtxout[kT] = mtxin[k];
            k++;
        }}
    }}
    return;
}}

/* shift a sequence of seqlen=N*veclen elements, veclen elements to the left
 * outseq[i] = seq[i + veclen], the last veclen elements are set to zero */
void {prefix}_mtx_shift_sequence(real_t outseq[], const real_t seq[],
        const uint32_t veclen, const uint32_t seqlen)
{{
    uint32_t i;
    /* shift vector one horizon step backwards */
    for (i=0; i < (seqlen - veclen); i++) {{
        outseq[i] = seq[i + veclen];
    }}

    /* set last element to zero */
    for (i=(seqlen - veclen); i < seqlen; i++) {{
        outseq[i] = 0.;
    }}
    return;
}}

/* identity matrix for mtx with size (dim x dim) */
void {prefix}_mtx_to_eye(real_t mtx[], const uint32_t dim)
{{
    uint32_t i, j;
    for (i = 0; i < dim; i++){{
        for (j = 0; j < dim; j++)
            mtx[i*dim+j] = (i == j) ? 1. : 0.;
    }}
}}

/* zero matrix for mtx with l entries */
void {prefix}_mtx_to_zero(real_t mtx[], const uint32_t l)
{{
    uint32_t i;
    for (i = 0; i < l; i++)
        mtx[i] = 0.;
}}

/* cholesky factorization: fac*fac^T = mtx.
 * mtx have size (dim x dim),
 * mtx is poitiv definite and symmetric,
 * chol_fac is a lower triangular matrix */
void {prefix}_mtx_cholesky(real_t sol[],
              const real_t mtx[], const uint32_t dim)
{{
    uint32_t i, j, k; /* loop counters */
    /* copy mtx to leave it unchanged*/
    for (i = 0; i < dim*dim; i++){{
        sol[i] = mtx[i];
    }}
    /* factorization */
    for (i = 0; i < dim; i++){{
        for (j = 0; j < i; j++){{
            sol[i*dim+i] -= sol[i*dim+j]*sol[i*dim+j];
        }}
#ifdef DEBUG_ALL_MODE
        if (sol[i*dim+i]<0){{
            printf(
                "FAILURE! in cholesky factorization\n(sqrt of a negativ %f)\n",
                sol[i*dim+i]);
        }}
#endif
        /* TODO fix compiler problem with sqrtf */
        sol[i*dim+i] = sqrt(sol[i*dim+i]);
        for (j = i+1; j < dim; j++){{
            for (k = 0; k < i; k++){{
                sol[j*dim+i] -= sol[j*dim+k]*sol[i*dim+k];
            }}
            sol[j*dim+i] /= sol[i*dim+i];             
        }}        
    }}
    /* change upper matrix entries to 0. */
    for (i = 0; i < dim; i++){{
        for (j = i+1; j < dim; j++){{
            sol[i*dim+j] = 0.;
        }}
    }}
}}

/* check if mtx is positiv definite and symmetric.
 * return 1 if condition is fullfilled, 0 otherwise
 * mtx have size (dim x dim) */
uint32_t {prefix}_mtx_check_chol(
                real_t fac[], const real_t mtx[], const uint32_t dim)
{{
    uint32_t ret = 1; /* return value */
    uint32_t i, j, k; /* loop counters */
    /* copy mtx to leave it unchanged*/
    for (i = 0; i < dim*dim; i++){{
        fac[i] = mtx[i];
    }}
    /* check condition "symmetric" */
    for (i = 0; i < dim; i++){{
        for (j = 0; j < dim; j++){{
            ret = (fac[i*dim+j] == fac[j*dim+i]) ? ret : 0;
        }}
    }}
    /* check factorization for condition "positiv definite" */
    for (i = 0; i < dim; i++){{
        for (j = 0; j < i; j++){{
            fac[i*dim+i] -= fac[i*dim+j]*fac[i*dim+j];
        }}
        ret = (fac[i*dim+i]<0) ? 0 : ret; /* cholesky whould fail */
        fac[i*dim+i] = sqrt((fac[i*dim+i]));
        for (j = i+1; j < dim; j++){{
            for (k = 0; k < i; k++){{
                fac[j*dim+i] -= fac[j*dim+k]*fac[i*dim+k];
            }}
            fac[j*dim+i] /= fac[i*dim+i];             
        }}
    }}
    return ret;
}}

/* forward substitution: mtxA * sol = mtxB
 * mtxA has shape (dim x dim), lower triangular matrix
 * sol and mtxB have shape (dim x colsB)*/
void {prefix}_mtx_fwd_subst(real_t sol[],
               const real_t mtxA[], const uint32_t dim,
               const real_t mtxB[], const uint32_t colsB)
{{
    uint32_t i, j, k; /* loop counters */
    for (k = 0; k < colsB; k++){{
        for (i = 0; i < dim; i++){{
            sol[colsB*i+k] = mtxB[colsB*i+k];
            for (j = 0; j < i; j++){{
                sol[colsB*i+k] -= mtxA[i*dim+j]*sol[colsB*j+k];
            }}
            sol[colsB*i+k] /= mtxA[i*dim+j];
        }}
    }}
    
}}

/* backward substitution: mtxA * sol = mtxB
 * mtxA has shape (dim x dim), upper triangular matrix
 * sol and mtxB have shape (dim x colsB)*/
void {prefix}_mtx_bwd_subst(real_t sol[],
               const real_t mtxA[], const uint32_t dim,
               const real_t mtxB[], const uint32_t colsB)
{{
    uint32_t i, j, k; /* loop counters */
    for (k = colsB; k > 0; k--){{
        for (i = 0; i < dim; i++){{
            sol[colsB*(dim-i)-k] = mtxB[colsB*(dim-i)-k];
            for (j = 0; j < i; j++){{
                sol[colsB*(dim-i)-k] -=
                        mtxA[dim*dim-i*dim-j-1]*sol[colsB*(dim-j)-k];
            }}
            sol[colsB*(dim-i)-k] /= mtxA[dim*dim-i*dim-j-1];
        }}
    }}
}}

/* elementwise comparision between mtxA and mtxB.
 * return -1 if |mtxA[i]-mtxB[i]| <= acc,
 * return last i for which the above condition is false otherwise */
uint32_t {prefix}_mtx_mtx_cmp(const real_t mtxA[], const real_t mtxB[], real_t dim,
                 real_t acc)
{{
    uint32_t i, cmp = -1;
    for (i = 0; i<dim; i++){{
        cmp = {prefix}_mtx_smpl_abs((mtxA[i] - mtxB[i])) <= acc ? cmp : i;
    }}
    return cmp;
}}

real_t {prefix}_mtx_smpl_exp(real_t e)
{{
    /*real_t half_pow, e_halbe;
    printf("in sp1\n");
    printf("in sp2 %f\n", e);
    if (e != e){{
         printf("in sp6 %f\n", e);
        return e;}}
    if (e == 0.){{
        return 1.;
        
    }}else if (e < 0){{
        printf("in sp5 %f\n", e_halbe);
        return 1 / smpl_pow(b, -e);
    }}else if (e > 0. && e < 1){{
//         printf("%f\n", (1./((uint32_t)(1/(1./((uint32_t)(1/e)) - e))) - (1./((uint32_t)(1/e)) - e)));
        printf("in sp\n");
        return nth_root(b[0], 1./e) / ( nth_root(b[0], 1./(1./((uint32_t)(1/e)) - e)) / nth_root(b[0], 1./(1./((uint32_t)(1/(1./((uint32_t)(1/e)) - e))) - (1./((uint32_t)(1/e)) - e))) );
    }}else if ((uint32_t)e % 2 == 0){{
        printf("in sp2\n");
        e_halbe = e/2;
        
        half_pow = smpl_pow(b, e_halbe);
        
        return half_pow * half_pow;
    }}else{{
        return b[0] * smpl_pow(b, e - 1);
    }}*/
    return exp(e); /* TODO fix compiler problem with expf */
}}

real_t {prefix}_mtx_smpl_abs(real_t x)
{{
    return (x < 0) ? -x : x;
}}