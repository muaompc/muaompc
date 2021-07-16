#include "stocmpcmtxops.h"
#include "stocmpccvp.h"


static void stocmpc_copy_data(struct stocmpc_term *dest, struct stocmpc_term *src);

void stocmpc_cvp_copy_parameters(struct stocmpc_cvp *cvp, 
                struct stocmpc_cvp_parameters *parameters) {
    uint32_t i;
    struct stocmpc_term ctl_par[STOCMPC_PAR_NUM];

    ctl_par[STOCMPC_X_K].data = parameters->x_k;

    for (i=0; i<STOCMPC_PAR_NUM; i++) {
        stocmpc_copy_data(cvp->par[i], &(ctl_par[i]));
    }
    return;
}

void stocmpc_cvp_form_problem(struct stocmpc_cvp *cvp)  {
    uint32_t i, j;
    struct stocmpc_pmetric *pm;

    for (i=0; i<STOCMPC_PMETRIC_NUM; i++) {
        pm = cvp->pmetric[i];
        stocmpc_copy_data(pm->val, pm->fac0);
        for (j=0; j<pm->fac_num[0]; j++) {
            stocmpc_mtx_mul_add(pm->val->data, pm->aux->data,
            pm->fac[j]->data, pm->par[j]->data,
            pm->fac[j]->rows, pm->fac[j]->cols);
        }
    }
    return;
}

void stocmpc_copy_data(struct stocmpc_term *dest, struct stocmpc_term *src)  {
    uint32_t j;
    for (j=0; j<(dest->cols*dest->rows); j++) {
        dest->data[j] = src->data[j];
    }
    return;
}

