#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */

#include "mc04types.h"
#include "cjson.h"
#include "{prefix}pbmdynmem.h"
#include "{prefix}mtxops.h"

/* Static functions declarations */

static {prefix}_dynmem_error_t {prefix}_pbm_parse_elements(
                struct {prefix}_pbm *pbm, cJSON *data);

static {prefix}_dynmem_error_t {prefix}_pbm_init_zeros(real_t **mtx, uint32_t rows, uint32_t cols);
static {prefix}_dynmem_error_t {prefix}_pbm_malloc(real_t **mtx, uint32_t rows, uint32_t cols);

/* External function definitions */

struct {prefix}_pbm *{prefix}_pbm_allocate_solver(void)
{{
    /* TODO compare to structure in tmpl_fgmdynmem.c.
     * TODO test changes of 22.02.2016. */
    struct {prefix}_pbm *pbm;
    struct {prefix}_pbm_conf *conf;
    struct {prefix}_pbm_P_hat *P;
    
    pbm = (struct {prefix}_pbm*)malloc(sizeof(struct {prefix}_pbm));
    if (NULL == pbm) {{return NULL;}}
    conf = (struct {prefix}_pbm_conf*)malloc(sizeof(struct {prefix}_pbm_conf));
    if (NULL == conf) {{return NULL;}}
    pbm->conf = conf;
    P = (struct {prefix}_pbm_P_hat*)malloc(sizeof(struct {prefix}_pbm_P_hat));
    if (NULL == P) {{return NULL;}}
    pbm->P_of_z = P;
    pbm->kappa = (real_t *)malloc(sizeof(real_t));
    if (NULL == pbm->kappa) {{return NULL;}}
    pbm->roh = (real_t *)malloc(sizeof(real_t));
    if (NULL == pbm->roh) {{return NULL;}}
    
    return pbm;
}}

/*
{prefix}_dynmem_error_t {prefix}_pbm_setup_solver(struct {prefix}_pbm *pbm,
                                            struct {prefix}_socp_prb *prb,
                                            char *fname)
*/
{prefix}_dynmem_error_t {prefix}_pbm_setup_solver(struct {prefix}_pbm *pbm,
                                            struct {prefix}_cvp_prb *prb,
                                            char *fname)
{{
    {prefix}_dynmem_error_t ret;
    cJSON *data;
    uint32_t i;
    real_t *tmp;
    
    struct {prefix}_pbm_block *bl;
    struct {prefix}_pbm_block **bls;
    
    data = {prefix}_dynmem_get_data(fname);
    if (NULL == data) {{return {PREFIX}_DYNMEM_FAIL;}}
    ret = {prefix}_pbm_parse_elements(pbm, data);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    
#if STOCMPC_CVP_PRB_SOCC
    pbm->b = prb->b->data;
    /* pbm->h = prb->h->data; */
    pbm->hsoft = prb->hsoft->data;
    pbm->Fusoft = prb->Fusoft->data;
    pbm->Fxsoft = prb->Fxsoft->data;
    pbm->Ffsoft = prb->Ffsoft->data;
    pbm->g = prb->g->data;
    pbm->A = prb->A->data;
    pbm->B = prb->B->data;
    pbm->C = prb->C->data;
    pbm->H = prb->H->data;
    pbm->z_ini = prb->z_ini->data;
    pbm->v_ini = prb->v_ini->data;
    pbm->zref = prb->zref->data;

    pbm->optvar_seqlen = pbm->optvar_veclen * pbm->horizon;
    pbm->dual_seqlen = pbm->state_veclen * pbm->horizon;
    pbm->control_veclen = pbm->optvar_veclen - pbm->state_veclen;
    pbm->sizeof_optvar_seqlen = sizeof(real_t) * pbm->optvar_seqlen;
    pbm->sizeof_dual_seqlen = sizeof(real_t) * pbm->dual_seqlen;
/*    
    pbm->z_ini = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->z_ini) {{return {PREFIX}_DYNMEM_FAIL;}}*/
    pbm->z_opt = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->z_opt) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->delta_z = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->delta_z) {{return {PREFIX}_DYNMEM_FAIL;}}
/*    
    pbm->v_ini = (real_t *)malloc(pbm->sizeof_dual_seqlen);
    if (NULL == pbm->v_ini) {{return {PREFIX}_DYNMEM_FAIL;}}*/
    pbm->v_opt = (real_t *)malloc(pbm->sizeof_dual_seqlen);
    if (NULL == pbm->v_opt) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->delta_v = (real_t *)malloc(pbm->sizeof_dual_seqlen);
    if (NULL == pbm->delta_v) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    /* u_k points on first entries in z_opt*/
    prb->u_k->data = pbm->z_opt;
    prb->u_k->rows = pbm->control_veclen;
    prb->u_k->cols = 1;
    
    pbm->P_of_z->h = prb->h->data;
    pbm->P_of_z->P = prb->P->data;
    pbm->P_of_z->nb_lin_constr = prb->P->rows;
    pbm->P_of_z->nb_socc = prb->nb_socc;
    pbm->P_of_z->nb_qc = prb->nb_qc;
    pbm->nb_of_ueq_constr = prb->P->rows+prb->nb_socc+prb->nb_qc;
    pbm->nb_of_soft_constr =
            (prb->Fusoft->rows)*pbm->horizon + prb->Ffsoft->rows;
           /* printf("nb of s%d\n", pbm->nb_of_soft_constr); */
    pbm->rowsFusoft = prb->Fusoft->rows;
    pbm->rowsFfsoft = prb->Ffsoft->rows;
    pbm->Psoft = prb->Psoft->data;
    pbm->Psoft_T =
            (real_t *)malloc(sizeof(real_t) * pbm->optvar_seqlen*pbm->nb_of_soft_constr);
    if (NULL == pbm->Psoft_T) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_transpose(pbm->Psoft_T, pbm->Psoft, prb->Psoft->rows, prb->Psoft->cols);
    
    pbm->P_of_z->h_hat =
            (real_t*)malloc(sizeof(real_t) * pbm->nb_of_ueq_constr);
    if (NULL == pbm->P_of_z->h_hat) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->P_of_z->socc = (struct {prefix}_pbm_socc**)calloc(prb->nb_socc, sizeof(struct {prefix}_pbm_socc*));
    if (NULL == pbm->P_of_z->socc) {{return {PREFIX}_DYNMEM_FAIL;}}
    for (i = 0; i < prb->nb_socc; i++){{
        pbm->P_of_z->socc[i] = (struct {prefix}_pbm_socc*)malloc(sizeof(struct {prefix}_pbm_socc));
        if (NULL == pbm->P_of_z->socc[i]) {{return {PREFIX}_DYNMEM_FAIL;}}
        pbm->P_of_z->socc[i]->rowsA = prb->socc[i]->A->rows;
        pbm->P_of_z->socc[i]->colsA = prb->socc[i]->A->cols;
        pbm->P_of_z->socc[i]->A = prb->socc[i]->A->data;
        pbm->P_of_z->socc[i]->A_T =
                (real_t *)malloc(sizeof(real_t) * prb->socc[i]->A->cols*prb->socc[i]->A->rows);
        if (NULL == pbm->P_of_z->socc[i]->A_T) {{return {PREFIX}_DYNMEM_FAIL;}}
        {prefix}_mtx_transpose(pbm->P_of_z->socc[i]->A_T, pbm->P_of_z->socc[i]->A,
                             prb->socc[i]->A->rows, prb->socc[i]->A->cols);
        pbm->P_of_z->socc[i]->b = prb->socc[i]->b->data;
        pbm->P_of_z->socc[i]->c = prb->socc[i]->c->data;
        pbm->P_of_z->socc[i]->d = prb->socc[i]->d->data;
        pbm->P_of_z->socc[i]->AAmcc=
                (real_t *)malloc(sizeof(real_t) * prb->socc[i]->A->cols*prb->socc[i]->A->cols);
        if (NULL == pbm->P_of_z->socc[i]->AAmcc) {{return {PREFIX}_DYNMEM_FAIL;}}
        {prefix}_mtx_multiply_mtx_mtx(pbm->P_of_z->socc[i]->AAmcc,
                                    pbm->P_of_z->socc[i]->A_T,
                                    pbm->P_of_z->socc[i]->A,
                                    pbm->P_of_z->socc[i]->colsA,
                                    pbm->P_of_z->socc[i]->rowsA,
                                    pbm->P_of_z->socc[i]->colsA);
        real_t cc[pbm->P_of_z->socc[i]->colsA*pbm->P_of_z->socc[i]->colsA];
        {prefix}_mtx_multiply_mtx_mtx(cc,
                                    pbm->P_of_z->socc[i]->c,
                                    pbm->P_of_z->socc[i]->c,
                                    pbm->P_of_z->socc[i]->colsA, 1,
                                    pbm->P_of_z->socc[i]->colsA);
        {prefix}_mtx_substract_direct(pbm->P_of_z->socc[i]->AAmcc, cc,
                                    pbm->P_of_z->socc[i]->colsA,
                                    pbm->P_of_z->socc[i]->colsA);
        
        pbm->P_of_z->socc[i]->par = pbm->z_opt+prb->socc[i]->par_0;
        pbm->P_of_z->socc[i]->par_0 = prb->socc[i]->par_0;
        pbm->P_of_z->socc[i]->par_l = prb->socc[i]->par_l;
        tmp = pbm->P_of_z->h_hat+prb->P->rows+prb->nb_qc+i;
        {prefix}_mtx_multiply_mtx_vec(tmp, pbm->P_of_z->socc[i]->b, pbm->P_of_z->socc[i]->b,
                                    1, pbm->P_of_z->socc[i]->rowsA);
        tmp[0] = pbm->P_of_z->socc[i]->d[0]*pbm->P_of_z->socc[i]->d[0] - tmp[0];
    }}
    pbm->P_of_z->qc = (struct {prefix}_pbm_qc**)calloc(prb->nb_qc, sizeof(struct {prefix}_pbm_qc*));
    if (NULL == pbm->P_of_z->qc) {{return {PREFIX}_DYNMEM_FAIL;}}
    for (i = 0; i < prb->nb_qc; i++){{
        pbm->P_of_z->qc[i] = (struct {prefix}_pbm_qc*)malloc(sizeof(struct {prefix}_pbm_qc));
        if (NULL == pbm->P_of_z->qc[i]) {{return {PREFIX}_DYNMEM_FAIL;}}
        pbm->P_of_z->qc[i]->dimGamma = prb->qc[i]->Gamma->rows;
        pbm->P_of_z->qc[i]->Gamma = prb->qc[i]->Gamma->data;
        pbm->P_of_z->qc[i]->beta = prb->qc[i]->beta->data;
        pbm->P_of_z->qc[i]->alpha = prb->qc[i]->alpha->data;
        pbm->P_of_z->qc[i]->par = pbm->z_opt+prb->qc[i]->par_0;
        pbm->P_of_z->qc[i]->par_0 = prb->qc[i]->par_0;
        pbm->P_of_z->qc[i]->par_l = prb->qc[i]->par_l;
        tmp = pbm->P_of_z->h_hat+prb->P->rows+i;
        tmp[0] = pbm->P_of_z->qc[i]->alpha[0];
    }}
    pbm->P_of_z->P_hat =
            (real_t *)malloc(sizeof(real_t) * pbm->nb_of_ueq_constr*pbm->optvar_seqlen);
    if (NULL == pbm->P_of_z->P_hat) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_to_zero(pbm->P_of_z->P_hat,
                         pbm->optvar_seqlen*pbm->nb_of_ueq_constr);
    pbm->P_of_z->P_hat_T =
            (real_t *)malloc(sizeof(real_t) * pbm->optvar_seqlen*pbm->nb_of_ueq_constr);
    if (NULL == pbm->P_of_z->P_hat_T) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->P_of_z->P2_hat =
            (real_t *)malloc(sizeof(real_t) * pbm->nb_of_ueq_constr*pbm->optvar_seqlen);
    if (NULL == pbm->P_of_z->P2_hat) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_to_zero(pbm->P_of_z->P2_hat,
                         pbm->optvar_seqlen*pbm->nb_of_ueq_constr);
    pbm->P_of_z->P2_hat_T =
            (real_t *)malloc(sizeof(real_t) * pbm->optvar_seqlen*pbm->nb_of_ueq_constr);
    if (NULL == pbm->P_of_z->P2_hat_T) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->h = pbm->P_of_z->h_hat;
    pbm->P = pbm->P_of_z->P_hat;
    pbm->P_T = pbm->P_of_z->P_hat_T;
    pbm->P2 = pbm->P_of_z->P2_hat;
    pbm->P2_T = pbm->P_of_z->P2_hat_T;
    
    pbm->j_in = &(pbm->conf->in_iter);
    pbm->reg = &(pbm->conf->reg);
    
    pbm->d = (real_t *)malloc(sizeof(real_t) * pbm->nb_of_ueq_constr);
    if (NULL == pbm->d) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->dsoft = (real_t *)malloc(sizeof(real_t) * pbm->nb_of_soft_constr);
    if (NULL == pbm->dsoft) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->diag_d_sq =
            (real_t *)malloc(sizeof(real_t) * pbm->nb_of_ueq_constr*pbm->nb_of_ueq_constr);
    if (NULL == pbm->diag_d_sq) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->diag_d_soft =
            (real_t *)malloc(sizeof(real_t) * pbm->nb_of_soft_constr*pbm->nb_of_soft_constr);
    if (NULL == pbm->diag_d_soft) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->Phi = (real_t *)malloc(sizeof(real_t) * pbm->optvar_seqlen*pbm->optvar_seqlen);
    if (NULL == pbm->Phi) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_to_zero(pbm->Phi, pbm->optvar_seqlen*pbm->optvar_seqlen);
    pbm->Phi_soft = (real_t *)malloc(sizeof(real_t) * pbm->optvar_seqlen*pbm->optvar_seqlen);
    if (NULL == pbm->Phi_soft) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_to_zero(pbm->Phi_soft,
                         pbm->optvar_seqlen*pbm->optvar_seqlen);
/*    
    pbm->P_T = (real_t *)malloc(sizeof(real_t) * prb->P->rows*prb->P->cols);
    if (NULL == pbm->P_T) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_transpose(pbm->P_T, pbm->P, prb->P->rows, prb->P->cols);
*/    
    pbm->A_T = (real_t *)malloc(sizeof(real_t) * prb->A->rows*prb->A->cols);
    if (NULL == pbm->A_T) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_transpose(pbm->A_T, pbm->A, prb->A->rows, prb->A->cols);
    
    pbm->B_T = (real_t *)malloc(sizeof(real_t) * prb->B->rows*prb->B->cols);
    if (NULL == pbm->B_T) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_transpose(pbm->B_T, pbm->B, prb->B->rows, prb->B->cols);
    
    pbm->A_B_T =
            (real_t *)malloc(sizeof(real_t) * (prb->A->cols+prb->B->cols)*prb->A->rows);
    if (NULL == pbm->A_B_T) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    for (i = 0; i < pbm->state_veclen*pbm->state_veclen; i++)
        pbm->A_B_T[i] = pbm->A_T[i];
    for (i = 0; i < pbm->state_veclen*pbm->control_veclen; i++)
        pbm->A_B_T[pbm->state_veclen*pbm->state_veclen+i] = pbm->B_T[i];
    pbm->A_B =
            (real_t *)malloc(sizeof(real_t) * prb->A->rows*(prb->A->cols+prb->B->cols));
    if (NULL == pbm->A_B) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_transpose(pbm->A_B, pbm->A_B_T, pbm->optvar_veclen, pbm->state_veclen);
    
    pbm->C_T = (real_t *)malloc(sizeof(real_t) * prb->C->rows*prb->C->cols);
    if (NULL == pbm->C_T) {{return {PREFIX}_DYNMEM_FAIL;}}
    {prefix}_mtx_transpose(pbm->C_T, pbm->C, prb->C->rows, prb->C->cols);

    pbm->r_d = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->r_d) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->r_d_soft = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->r_d_soft) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->r_p = (real_t *)malloc(pbm->sizeof_dual_seqlen);
    if (NULL == pbm->r_p) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->st_size = (real_t *)malloc(sizeof(real_t));
    if (NULL == pbm->st_size) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp1_optvar_seqlen = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->tmp1_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp2_optvar_seqlen = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->tmp2_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp2_dual_seqlen = (real_t *)malloc(pbm->sizeof_dual_seqlen);
    if (NULL == pbm->tmp2_dual_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp3_mtx_optvar_nb_of_ueq = (real_t *)malloc(sizeof(real_t) * pbm->optvar_seqlen*pbm->nb_of_ueq_constr);
    if (NULL == pbm->tmp3_mtx_optvar_nb_of_ueq) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp3_mtx_optvar_nb_of_soft = (real_t *)malloc(sizeof(real_t) * pbm->optvar_seqlen*pbm->nb_of_soft_constr);
    if (NULL == pbm->tmp3_mtx_optvar_nb_of_ueq) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp3_state_veclen = (real_t *)malloc(sizeof(real_t) * pbm->state_veclen);
    if (NULL == pbm->tmp3_state_veclen) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp4_mtx_optvar_optvar = (real_t *)malloc(sizeof(real_t) * pbm->optvar_seqlen*pbm->optvar_seqlen);
    if (NULL == pbm->tmp4_mtx_optvar_optvar) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp5_nb_of_constr = (real_t *)malloc(sizeof(real_t) * pbm->nb_of_ueq_constr);
    if (NULL == pbm->tmp5_nb_of_constr) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp6_optvar_seqlen = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->tmp6_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp7_dual_seqlen = (real_t *)malloc(pbm->sizeof_dual_seqlen);
    if (NULL == pbm->tmp7_dual_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp8_L_Y =
            (real_t *)malloc(sizeof(real_t) * (2*pbm->horizon-1)*pbm->state_veclen*pbm->state_veclen);
    if (NULL == pbm->tmp8_L_Y) {{return {PREFIX}_DYNMEM_FAIL;}}  /*T*n*n* big enough*/
    
    pbm->tmp9_L_Y_T =
            (real_t *)malloc(sizeof(real_t) * (2*pbm->horizon-1)*pbm->state_veclen*pbm->state_veclen);
    if (NULL == pbm->tmp9_L_Y_T) {{return {PREFIX}_DYNMEM_FAIL;}}  /*T*n*n* big enough*/
        
    pbm->tmp8_L_Phi =
            (real_t *)malloc(sizeof(real_t) * pbm->horizon*pbm->optvar_veclen*pbm->optvar_veclen);
    if (NULL == pbm->tmp8_L_Phi) {{return {PREFIX}_DYNMEM_FAIL;}}  /*m*m+(T-1)*(n+m)*(n+m)+n*n big enough*/
    
    pbm->tmp9_L_Phi_T =
            (real_t *)malloc(sizeof(real_t) * pbm->horizon*pbm->optvar_veclen*pbm->optvar_veclen);
    if (NULL == pbm->tmp9_L_Phi_T) {{return {PREFIX}_DYNMEM_FAIL;}}  /*m*m+(T-1)*(n+m)*(n+m)+n*n big enough*/
    
    pbm->tmp_phibl1 = (real_t *)malloc(sizeof(real_t) * pbm->optvar_veclen*pbm->optvar_veclen);
    if (NULL == pbm->tmp_phibl1) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp_phibl2 = (real_t *)malloc(sizeof(real_t) * pbm->optvar_veclen*pbm->optvar_veclen);
    if (NULL == pbm->tmp_phibl2) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp_phibl3 = (real_t *)malloc(sizeof(real_t) * pbm->optvar_veclen*pbm->optvar_veclen);
    if (NULL == pbm->tmp_phibl3) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp10 = (real_t *)malloc(sizeof(real_t) * pbm->optvar_veclen*pbm->optvar_veclen);
    if (NULL == pbm->tmp10) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmpYbl = (real_t *)malloc(sizeof(real_t) * pbm->state_veclen*pbm->state_veclen);
    if (NULL == pbm->tmpYbl) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmpQbl = (real_t *)malloc(sizeof(real_t) * pbm->state_veclen*pbm->state_veclen);
    if (NULL == pbm->tmpQbl) {{return {PREFIX}_DYNMEM_FAIL;}}
    
    pbm->tmp1_res_os = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->tmp1_res_os) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->tmp2_res_os = (real_t *)malloc(pbm->sizeof_optvar_seqlen);
    if (NULL == pbm->tmp2_res_os) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->tmp3_res_ds = (real_t *)malloc(pbm->sizeof_dual_seqlen);
    if (NULL == pbm->tmp3_res_ds) {{return {PREFIX}_DYNMEM_FAIL;}}
        
    pbm->eye_optvar_veclen = (real_t *)malloc(sizeof(real_t) * pbm->optvar_veclen*pbm->optvar_veclen);
    if (NULL == pbm->eye_optvar_veclen) {{return {PREFIX}_DYNMEM_FAIL;}}
     {prefix}_mtx_to_eye(pbm->eye_optvar_veclen, pbm->optvar_veclen); 
    
    pbm->eye_state_veclen = (real_t *)malloc(sizeof(real_t) * pbm->state_veclen*pbm->state_veclen);
    if (NULL == pbm->eye_state_veclen) {{return {PREFIX}_DYNMEM_FAIL;}}
     {prefix}_mtx_to_eye(pbm->eye_state_veclen, pbm->state_veclen); 
    
    pbm->sizeof_optvar_veclen = sizeof(real_t) * pbm->optvar_veclen;
    
    bl = (struct {prefix}_pbm_block*)malloc(sizeof(struct {prefix}_pbm_block));
    if (NULL == bl) {{return {PREFIX}_DYNMEM_FAIL;}}
    bl->rows = pbm->optvar_veclen;
    bl->cols = pbm->nb_of_soft_constr;
    bl->data = (real_t *)malloc(pbm->sizeof_optvar_veclen * bl->cols);
    if (NULL == bl->data) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->tmp_Phi_sft_blk = bl;
    bls = (struct {prefix}_pbm_block**)calloc(pbm->horizon + 1, sizeof(struct {prefix}_pbm_block*));
    if (NULL == bls) {{return {PREFIX}_DYNMEM_FAIL;}}
    bls[0] = (struct {prefix}_pbm_block*)malloc(sizeof(struct {prefix}_pbm_block));
    if (NULL == bls[0]) {{return {PREFIX}_DYNMEM_FAIL;}}
    bls[0]->rows = pbm->control_veclen;
    bls[0]->cols = pbm->control_veclen;
    bls[0]->data = (real_t *)malloc(sizeof(real_t) * bls[0]->rows * bls[0]->cols);
    if (NULL == bls[0]->data) {{return {PREFIX}_DYNMEM_FAIL;}}
    for (i = 1; i < pbm->horizon; i++){{
        bls[i] = (struct {prefix}_pbm_block*)malloc(sizeof(struct {prefix}_pbm_block));
        if (NULL == bls[i]) {{return {PREFIX}_DYNMEM_FAIL;}}
        bls[i]->rows = pbm->optvar_veclen;
        bls[i]->cols = pbm->optvar_veclen;
        bls[i]->data = (real_t *)malloc(sizeof(real_t) * bls[i]->rows * bls[i]->cols);
        if (NULL == bls[i]->data) {{return {PREFIX}_DYNMEM_FAIL;}}
    }}
    bls[i] = (struct {prefix}_pbm_block*)malloc(sizeof(struct {prefix}_pbm_block));
    if (NULL == bls[i]) {{return {PREFIX}_DYNMEM_FAIL;}}
    bls[i]->rows = pbm->state_veclen;
    bls[i]->cols = pbm->state_veclen;
    bls[i]->data = (real_t *)malloc(sizeof(real_t) * bls[i]->rows * bls[i]->cols);
    if (NULL == bls[i]->data) {{return {PREFIX}_DYNMEM_FAIL;}}
    pbm->Phi_sft_blks = bls;

    #endif
    
    return {PREFIX}_DYNMEM_OK;
}}

/* Static function definitions */

{prefix}_dynmem_error_t {prefix}_pbm_parse_elements(struct {prefix}_pbm *pbm, cJSON *data)
{{
    cJSON *kappa, *roh, *optvar, *veclen, *horizon, *state_veclen;
    
    kappa = cJSON_GetObjectItem(data, "kappa");
    if (NULL == kappa) {{
        printf("ERROR: could not parse item %s \n", "kappa");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    *(pbm->kappa) = (real_t)kappa->valuedouble;
    
    roh = cJSON_GetObjectItem(data, "roh");
    if (NULL == roh) {{
        printf("ERROR: could not parse item %s \n", "roh");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    *(pbm->roh) = (real_t)roh->valuedouble;
    
    optvar = cJSON_GetObjectItem(data, "optvar");
    if (NULL == optvar) {{
        printf("ERROR: could not parse item %s \n", "optvar");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    veclen = cJSON_GetObjectItem(optvar, "veclen");
    if (NULL == veclen) {{
        printf("ERROR: could not parse item %s \n", "veclen");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    pbm->optvar_veclen = (uint32_t)veclen->valueint;
    
    horizon = cJSON_GetObjectItem(optvar, "horizon");
    if (NULL == horizon) {{
        printf("ERROR: could not parse item %s \n", "horizon");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    pbm->horizon = (uint32_t)horizon->valueint;
    
    state_veclen = cJSON_GetObjectItem(optvar, "state_veclen");
    if (NULL == state_veclen) {{
        printf("ERROR: could not parse item %s \n", "state_veclen");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    pbm->state_veclen = (uint32_t)state_veclen->valueint;
    
    return {PREFIX}_DYNMEM_OK;
}}

static {prefix}_dynmem_error_t {prefix}_pbm_init_zeros(real_t **mtx, uint32_t rows, uint32_t cols)
{{
    if ({prefix}_pbm_malloc(mtx, rows, cols)) {{
        return {PREFIX}_DYNMEM_FAIL;
    }}
    mpc_mtx_to_zero(*mtx, rows*cols);
    return {PREFIX}_DYNMEM_OK;
}}

static {prefix}_dynmem_error_t {prefix}_pbm_malloc(real_t **mtx, uint32_t rows, uint32_t cols)
{{
    *mtx = (real_t *)malloc(sizeof(real_t) * rows*cols);
    if (*mtx == NULL) {{
        return {PREFIX}_DYNMEM_FAIL;
    }}
    return {PREFIX}_DYNMEM_OK;
}}