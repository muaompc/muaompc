#ifndef {PREFIX}_MATH_H
#define {PREFIX}_MATH_H

#include <math.h>

#include "arithmetic.h"
#include "{prefix}mtxops.h"

/* External function declarations */

/* cholesky factorization: chol_fac*chol_fac^T = mtx
 * mtx have size (dimension x dimension),
 * mtx is poitiv definite and symmetric,
 * chol_fac is a lower triangular matrix */
void cholesky(real_t solution[],
              const real_t mtx[], const uint32_t dimension);

/* check if mtx is positiv definite and symmetric.
 * return 1 if condition is fullfilled, 0 otherwise
 * mtx have size (dimension x dimension) */
uint32_t check_chol(real_t chol_fac[],
                    const real_t mtx[], const uint32_t dimension);

/* forward substitution: mtxA * solution = mtxB
 * mtxA has shape (dimA x dimA), lower Dreiecksmatrix
 * solution and mtxB have shape (dimA x columsB)*/
void fwd_subst(real_t solution[],
               const real_t mtxA[], const uint32_t dimensionA,
               const real_t mtxB[], const uint32_t columsB);

/* backward substitution: mtxA * solution = mtxB
 * mtxA has shape (dimA x dimA), upper Dreiecksmatrix
 * solution and mtxB have shape (dimA x columsB)*/
void bwd_subst(real_t solution[],
               const real_t mtxA[], const uint32_t dimensionA,
               const real_t mtxB[], const uint32_t columsB);

/* elementwise comparision between mtxA and mtxB.
 * return -1 if |mtxA[i]-mtxB[i]| <= accuracy,
 * return last i for which the above condition is false otherwise */
uint32_t mtx_cmp(const real_t mtxA[], const real_t mtxB[], real_t dim,
                 real_t accuracy);
                
real_t {prefix}_smpl_exp(real_t exp);

real_t smpl_abs(real_t x);

#endif  /* {PREFIX}_MATH_H */