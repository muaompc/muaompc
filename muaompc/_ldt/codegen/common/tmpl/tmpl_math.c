#define DEBUG_ALL_MODE

#ifdef DEBUG_ALL_MODE
#include "stdio.h"
#endif
#include "{prefix}math.h"


/* cholesky factorization: fac*fac^T = mtx.
 * mtx have size (dim x dim),
 * mtx is poitiv definite and symmetric,
 * chol_fac is a lower triangular matrix */
void cholesky(real_t sol[],
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
uint32_t {prefix}_check_chol(
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
void fwd_subst(real_t sol[],
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
void bwd_subst(real_t sol[],
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
uint32_t mtx_cmp(const real_t mtxA[], const real_t mtxB[], real_t dim,
                 real_t acc)
{{
    uint32_t i, cmp = -1;
    for (i = 0; i<dim; i++){{
        cmp = smpl_abs((mtxA[i] - mtxB[i])) <= acc ? cmp : i;
    }}
    return cmp;
}}

real_t {prefix}_smpl_exp(real_t e)
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

real_t smpl_abs(real_t x)
{{
    return (x < 0) ? -x : x;
}}
