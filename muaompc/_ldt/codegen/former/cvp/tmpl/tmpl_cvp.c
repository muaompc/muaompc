#include "{prefix}mtxops.h"
#include "{prefix}cvp.h"


static void {prefix}_cvp_form_pmetric(struct {prefix}_pmetric *pmt);
static void {prefix}_copy_data(struct {prefix}_term *dest, struct {prefix}_term *src);

void {prefix}_cvp_copy_parameters(struct {prefix}_cvp *cvp, 
                struct {prefix}_cvp_parameters *parameters) {{
    uint32_t i;
    struct {prefix}_term ctl_par[{PREFIX}_PAR_NUM];

{copy_parameters}
    for (i=0; i<{PREFIX}_PAR_NUM; i++) {{
        {prefix}_copy_data(cvp->par[i], &(ctl_par[i]));
    }}
    return;
}}

void {prefix}_cvp_form_problem(struct {prefix}_cvp *cvp)  {{
    uint32_t i, j;
    for (i=0; i<{PREFIX}_PMETRIC_NUM; i++) {{
        if (cvp->pmetric[i]->fac_num[0] > 0) {{
            {prefix}_cvp_form_pmetric(cvp->pmetric[i]);
        }}  /* else: pmetric[i] is constant */
    }}
    for (j=0; j<cvp->socc_num[0]; j++) {{
        for (i=0; i<{PREFIX}_SOCC_PMETRIC_NUM; i++) {{
            if (cvp->socc[j]->pmetric[i]->fac_num[0] > 0) {{
                {prefix}_cvp_form_pmetric(cvp->socc[j]->pmetric[i]);
            }}
        }}
    }}
    return;
}}

void {prefix}_cvp_form_pmetric(struct {prefix}_pmetric *pmt)  {{
    uint32_t j;

    {prefix}_copy_data(pmt->val, pmt->fac0);
    for (j=0; j<pmt->fac_num[0]; j++) {{
        {prefix}_mtx_mul_add(pmt->val->data, pmt->aux->data,
        pmt->fac[j]->data, pmt->par[j]->data,
        pmt->fac[j]->rows, pmt->fac[j]->cols);
    }}
    return;
}}

void {prefix}_copy_data(struct {prefix}_term *dest, struct {prefix}_term *src)  {{
    uint32_t j;
    for (j=0; j<(dest->cols*dest->rows); j++) {{
        dest->data[j] = src->data[j];
    }}
    return;
}}

