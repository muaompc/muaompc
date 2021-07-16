#include "{prefix}{dataname}pbmdata.h"

/* #include "../../src/static_data.h" */
#if {PREFIX}_CVP_PRB_SPARSE
uint32_t {prefix}{dataname}_horizon = {pbm_horizon};
uint32_t {prefix}{dataname}_state_veclen = {state_veclen};
uint32_t {prefix}{dataname}_optvar_veclen = HHMPC_OV;
uint32_t {prefix}{dataname}_optvar_seqlen = HHMPC_OS;
uint32_t {prefix}{dataname}_nb_qc =  HHMPC_NB_QC;
uint32_t {prefix}{dataname}_nb_socc = HHMPC_NB_SOCC;
#else
/* TODO Generate Values. */
uint32_t {prefix}{dataname}_horizon = {pbm_horizon};
uint32_t {prefix}{dataname}_state_veclen = {state_veclen};
uint32_t {prefix}{dataname}_optvar_veclen = {pbm_optvar_veclen}; 
uint32_t {prefix}{dataname}_optvar_seqlen = {optvar_seqlen};
uint32_t {prefix}{dataname}_nb_qc =  0;  /* TODO Set nb_qc. */
uint32_t {prefix}{dataname}_nb_socc = {number_socc};
#endif
struct {prefix}_pbm_conf {prefix}{dataname}_pbm_conf = {{1, 1, 0.01}};
struct {prefix}_pbm_P_hat {prefix}{dataname}_pbm_P;
real_t {prefix}{dataname}_kappa = 90.;
real_t {prefix}{dataname}_roh = .1;

real_t {prefix}{dataname}_zr_data[] = {zeros_optvar_seqlen};
struct {prefix}_term {prefix}{dataname}_zr_term =
                {{{optvar_seqlen}, 1, {prefix}{dataname}_zr_data}};

real_t {prefix}{dataname}_zini_data[] = {zeros_optvar_seqlen};
struct {prefix}_term {prefix}{dataname}_zini_term =
                {{{optvar_seqlen}, 1, {prefix}{dataname}_zini_data}};
                
real_t {prefix}{dataname}_vini_data[] = {zeros_dual_seqlen};
struct {prefix}_term {prefix}{dataname}_vini_term =
                {{{dual_seqlen}, 1, {prefix}{dataname}_vini_data}};

real_t {prefix}{dataname}_z_opt[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_delta_z[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_v_opt[] = {zeros_dual_seqlen};
real_t {prefix}{dataname}_delta_v[] = {zeros_dual_seqlen};

/* TODO Generate pmetric b in cvp. */
real_t {prefix}{dataname}_pmetric_b_data[] = {zeros_dual_seqlen};
struct {prefix}_term {prefix}{dataname}_pmetric_b_term =
                {{{dual_seqlen}, 1, {prefix}{dataname}_pmetric_b_data}};

struct {prefix}_term {prefix}{dataname}_gcond_term;
struct {prefix}_term {prefix}{dataname}_Hcond_term;
real_t {prefix}{dataname}_A_data[] = {pbm_a_data};
struct {prefix}_term {prefix}{dataname}_A_term =
                {{{state_veclen}, {state_veclen}, {prefix}{dataname}_A_data}};
real_t {prefix}{dataname}_B_data[] = {pbm_b_data};
struct {prefix}_term {prefix}{dataname}_B_term =
                {{{state_veclen}, {control_veclen}, {prefix}{dataname}_B_data}};
real_t {prefix}{dataname}_C_data[] = {pbm_c_data};
struct {prefix}_term {prefix}{dataname}_C_term =
                {{{dual_seqlen}, {optvar_seqlen}, {prefix}{dataname}_C_data}};

/* TODO Generate b correctly in cvp. */
real_t {prefix}{dataname}_h_val_data[] = {zeros_affine};
struct {prefix}_term {prefix}{dataname}_h_val_term =
                {{{number_affine}, 1, {prefix}{dataname}_h_val_data}};
#if {PREFIX}_CVP_PRB_AFFINE
real_t {prefix}{dataname}_P_data[] = {pbm_p_data};
#else
real_t {prefix}{dataname}_P_data[] ={{0}};  /* TODO HHMPC_NB_LCONSTR*HHMPC_OS*/
#endif
struct {prefix}_term {prefix}{dataname}_P_term =
                {{{number_affine}, {optvar_seqlen}, {prefix}{dataname}_P_data}};
real_t {prefix}{dataname}_h_hat[] = {zeros_ineq};

#if {PREFIX}_CVP_PRB_SOCC
struct {prefix}_pbm_socc *{prefix}{dataname}_psocc[{number_socc}];
struct {prefix}_pbm_socc {prefix}{dataname}_pbm_socc[{number_socc}]; 
real_t soccA_T1[2*2];  /* TODO Autogenerate size number of socc variables. */
real_t soccA_T2[2*2];
real_t *soccA_T[{number_socc}];
real_t soccAAmcc1[2*2];
real_t soccAAmcc2[2*2];
real_t *soccAAmcc[{number_socc}];
#endif

real_t {prefix}{dataname}_P_hat[] = {zeros_optvar_seqlen_x_number_ineq};
real_t {prefix}{dataname}_P_hat_T[] = {zeros_optvar_seqlen_x_number_ineq};
real_t {prefix}{dataname}_P2_hat[] = {zeros_optvar_seqlen_x_number_ineq};
real_t {prefix}{dataname}_P2_hat_T[] = {zeros_optvar_seqlen_x_number_ineq};

real_t {prefix}{dataname}_d[] = {zeros_ineq};
real_t {prefix}{dataname}_dsoft[] = {zeros_soft};
real_t {prefix}{dataname}_diag_d_sq[] = {zeros_ineq_x_ineq};
real_t {prefix}{dataname}_diag_d_soft[] = {zeros_soft_x_soft};
real_t {prefix}{dataname}_Phi[] = {zeros_optvar_seqlen_x_optvar_seqlen};
real_t {prefix}{dataname}_Phi_soft[] = {zeros_optvar_seqlen_x_optvar_seqlen};

real_t {prefix}{dataname}_A_T[] = {zeros_state_veclen_x_state_veclen};
real_t {prefix}{dataname}_B_T[] = {zeros_state_veclen_x_control_veclen};
real_t {prefix}{dataname}_A_B_T[] = {zeros_state_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_A_B[] = {zeros_state_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_C_T[] = {zeros_optvar_seqlen_x_dual_seqlen};

real_t {prefix}{dataname}_r_d[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_r_d_soft[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_r_p[] = {zeros_dual_seqlen};

real_t {prefix}{dataname}_st_size[] = {{0.}};

real_t {prefix}{dataname}_tmp1_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_tmp2_optvar_seqlen[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_tmp2_dual_seqlen[] = {zeros_dual_seqlen};
real_t {prefix}{dataname}_tmp3_mtx_optvar_nb_of_ueq[] = {zeros_optvar_seqlen_x_number_ineq};
real_t {prefix}{dataname}_tmp3_mtx_optvar_nb_of_soft[] = {zeros_optvar_seqlen_x_soft};
real_t {prefix}{dataname}_tmp3_state_veclen[] = {zeros_state_veclen};
real_t {prefix}{dataname}_tmp4_mtx_optvar_optvar[] = {zeros_optvar_seqlen_x_optvar_seqlen};
real_t {prefix}{dataname}_tmp5_nb_of_constr[] = {zeros_ineq};
real_t {prefix}{dataname}_tmp6_optvar_seqlen[]  = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_tmp7_dual_seqlen[] = {zeros_dual_seqlen};
real_t {prefix}{dataname}_tmp8_L_Y[] = {zeros_2_x_horizon_x_state_veclen_x_state_veclen};
real_t {prefix}{dataname}_tmp9_L_Y_T[] = {zeros_2_x_horizon_x_state_veclen_x_state_veclen};
real_t {prefix}{dataname}_tmp8_L_Phi[] = {pbm_zeros_horizon_x_optvar_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_tmp9_L_Phi_T[] = {pbm_zeros_horizon_x_optvar_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_tmp_phibl1[] = {pbm_zeros_optvar_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_tmp_phibl2[] = {pbm_zeros_optvar_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_tmp_phibl3[] = {pbm_zeros_optvar_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_tmp10[] = {pbm_zeros_optvar_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_tmpYbl[] = {zeros_state_veclen_x_state_veclen};
real_t {prefix}{dataname}_tmpQbl[] = {zeros_state_veclen_x_state_veclen};
real_t {prefix}{dataname}_tmp1_res_os[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_tmp2_res_os[] = {zeros_optvar_seqlen};
real_t {prefix}{dataname}_tmp3_res_ds[] = {zeros_dual_seqlen};

real_t {prefix}{dataname}_eye_optvar_veclen[] = {pbm_zeros_optvar_veclen_x_optvar_veclen};
real_t {prefix}{dataname}_eye_state_veclen[] = {zeros_state_veclen_x_state_veclen};
                
              
void {prefix}{dataname}_zeroes(real_t mtx[], const uint32_t l)
{{
    uint32_t i;
    for (i = 0; i < l; i++)
        mtx[i] = 0.;
}}

/* Setup pbm. */
void {prefix}{dataname}_pbm_setup_solver(
                struct {prefix}_pbm *pbm,
                struct {prefix}_cvp_prb *prb)
{{
    uint32_t i, j;
    real_t *tmp, *cc;
    pbm->rowsFusoft = 0;  /* TODO Right number. */
    pbm->rowsFfsoft = 0;
    
    pbm->horizon = {prefix}{dataname}_horizon;
    pbm->state_veclen = {prefix}{dataname}_state_veclen;
    pbm->optvar_veclen = {prefix}{dataname}_optvar_veclen;
    pbm->optvar_seqlen = {prefix}{dataname}_optvar_seqlen;
    
    pbm->dual_seqlen = pbm->state_veclen * pbm->horizon;
    pbm->control_veclen = pbm->optvar_veclen - pbm->state_veclen;
    pbm->sizeof_optvar_seqlen = sizeof(real_t) * pbm->optvar_seqlen;
    pbm->sizeof_dual_seqlen = sizeof(real_t) * pbm->dual_seqlen;
    
    pbm->conf = &{prefix}{dataname}_pbm_conf;
    pbm->kappa = &{prefix}{dataname}_kappa;
    pbm->roh = &{prefix}{dataname}_roh;
    pbm->P_of_z = &{prefix}{dataname}_pbm_P;
    
    pbm->P_of_z->h = {prefix}{dataname}_h_val_term.data;
    pbm->P_of_z->P = {prefix}{dataname}_P_term.data;
    /* Next line is now calulated in pbm/codegen.py
    {prefix}_mtx_scale_direct(pbm->P_of_z->P, -1, HHMPC_NB_LCONSTR/2, HHMPC_OS);*/
    pbm->P_of_z->nb_lin_constr = {prefix}{dataname}_P_term.rows;
    pbm->P_of_z->nb_socc = {prefix}{dataname}_nb_socc;
    pbm->P_of_z->nb_qc = {prefix}{dataname}_nb_qc;
    
    pbm->nb_of_ueq_constr = {prefix}{dataname}_P_term.rows +
                    {prefix}{dataname}_nb_socc + {prefix}{dataname}_nb_qc;
    pbm->nb_of_soft_constr = 0;/*
                (prb->Fusoft->rows)*ipm->horizon + prb->Ffsoft->rows;*/
                
    pbm->P_of_z->h_hat = {prefix}{dataname}_h_hat;    
    
    /*TODO zeroes */
    pbm->P_of_z->P_hat = {prefix}{dataname}_P_hat;
    {prefix}{dataname}_zeroes(pbm->P_of_z->P_hat, pbm->optvar_seqlen*pbm->nb_of_ueq_constr);
    pbm->P_of_z->P_hat_T = {prefix}{dataname}_P_hat_T;
    {prefix}{dataname}_zeroes(pbm->P_of_z->P_hat_T, pbm->nb_of_ueq_constr*pbm->optvar_seqlen);
    pbm->P_of_z->P2_hat = {prefix}{dataname}_P2_hat;
    {prefix}{dataname}_zeroes(pbm->P_of_z->P2_hat, pbm->optvar_seqlen*pbm->nb_of_ueq_constr);
    pbm->P_of_z->P2_hat_T = {prefix}{dataname}_P2_hat_T;
    {prefix}{dataname}_zeroes(pbm->P_of_z->P2_hat_T, pbm->nb_of_ueq_constr*pbm->optvar_seqlen);
    
    pbm->h = pbm->P_of_z->h_hat;
    pbm->P = pbm->P_of_z->P_hat;
    pbm->P_T = pbm->P_of_z->P_hat_T;
    pbm->P2 = pbm->P_of_z->P2_hat;
    pbm->P2_T = pbm->P_of_z->P2_hat_T;
    
    pbm->j_in = &(pbm->conf->in_iter);
    pbm->reg = &(pbm->conf->reg);
    
    pbm->z_ini = {prefix}{dataname}_zini_term.data;
    pbm->v_ini = {prefix}{dataname}_vini_term.data;
    pbm->zref = {prefix}{dataname}_zr_term.data;
    
    pbm->z_opt = {prefix}{dataname}_z_opt;
    pbm->delta_z = {prefix}{dataname}_delta_z;
    pbm->v_opt = {prefix}{dataname}_v_opt;
    pbm->delta_v = {prefix}{dataname}_delta_v;
    
    pbm->b = {prefix}{dataname}_pmetric_b_term.data;
    /* TODO */
    /* ipm->hsoft = prb->hsoft->data;
    ipm->Fusoft = prb->Fusoft->data;
    ipm->Fxsoft = prb->Fxsoft->data;
    ipm->Ffsoft = prb->Ffsoft->data; */
    
    pbm->g = {prefix}{dataname}_gcond_term.data;
    pbm->g = prb->g->data;
    pbm->H = {prefix}{dataname}_Hcond_term.data;
    pbm->H = prb->H->data;
    
    pbm->A = {prefix}{dataname}_A_term.data;
    pbm->B = {prefix}{dataname}_B_term.data;
    pbm->C = {prefix}{dataname}_C_term.data;
    
    /*TODO zeroes */
    pbm->d = {prefix}{dataname}_d;
    pbm->dsoft = {prefix}{dataname}_dsoft;
    pbm->diag_d_sq = {prefix}{dataname}_diag_d_sq;
    {prefix}{dataname}_zeroes(pbm->diag_d_sq, pbm->nb_of_ueq_constr*pbm->nb_of_ueq_constr);
    pbm->diag_d_soft = {prefix}{dataname}_diag_d_soft;
    pbm->Phi = {prefix}{dataname}_Phi;
    {prefix}{dataname}_zeroes(pbm->Phi, pbm->optvar_seqlen*pbm->optvar_seqlen);
    pbm->Phi_soft = {prefix}{dataname}_Phi_soft;
    {prefix}{dataname}_zeroes(pbm->Phi_soft, pbm->optvar_seqlen*pbm->optvar_seqlen);
    pbm->A_T = {prefix}{dataname}_A_T;
    {prefix}_mtx_transpose(pbm->A_T, pbm->A, {prefix}{dataname}_A_term.rows,
                           {prefix}{dataname}_A_term.cols);
    pbm->B_T = {prefix}{dataname}_B_T;
    {prefix}_mtx_transpose(pbm->B_T, pbm->B, {prefix}{dataname}_B_term.rows,
                           {prefix}{dataname}_B_term.cols);
    pbm->A_B_T = {prefix}{dataname}_A_B_T;
    for (i = 0; i < pbm->state_veclen*pbm->state_veclen; i++)
        pbm->A_B_T[i] = pbm->A_T[i];
    for (i = 0; i < pbm->state_veclen*pbm->control_veclen; i++)
        pbm->A_B_T[pbm->state_veclen*pbm->state_veclen+i] = pbm->B_T[i];
    pbm->A_B = {prefix}{dataname}_A_B;
    {prefix}_mtx_transpose(pbm->A_B, pbm->A_B_T, pbm->optvar_veclen, pbm->state_veclen);
    pbm->C_T = {prefix}{dataname}_C_T;
    {prefix}_mtx_transpose(pbm->C_T, pbm->C, {prefix}{dataname}_C_term.rows,
                           {prefix}{dataname}_C_term.cols);
    
#if STOCMPC_CVP_PRB_SOCC
    pbm->P_of_z->socc = {prefix}{dataname}_psocc;
    
    soccA_T[0] = soccA_T1;
    soccA_T[1] = soccA_T2;
    soccAAmcc[0] = soccAAmcc1;
    soccAAmcc[1] = soccAAmcc2;
    
    for (i = 0; i < pbm->P_of_z->nb_socc; i++){{
        pbm->P_of_z->socc[i] = &{prefix}{dataname}_pbm_socc[i];
        pbm->P_of_z->socc[i]->rowsA = prb->socc[i]->Wm->rows;
        pbm->P_of_z->socc[i]->colsA = prb->socc[i]->Wm->cols;
        pbm->P_of_z->socc[i]->b = prb->socc[i]->wn->data;
        pbm->P_of_z->socc[i]->c = prb->socc[i]->wvT->data;
        pbm->P_of_z->socc[i]->d = prb->socc[i]->ws->data;
        
        /* TODO Define correct parameter. */
        pbm->P_of_z->socc[i]->A = prb->socc[i]->Wm->data;
        pbm->P_of_z->socc[i]->par = pbm->z_opt+0;
        pbm->P_of_z->socc[i]->par_0 = 0;
        pbm->P_of_z->socc[i]->par_l = {prefix}{dataname}_optvar_seqlen;
        pbm->P_of_z->socc[i]->A_T = soccA_T[i];
        {prefix}_mtx_transpose(pbm->P_of_z->socc[i]->A_T, pbm->P_of_z->socc[i]->A,
                             prb->socc[i]->Wm->rows, prb->socc[i]->Wm->cols);
        pbm->P_of_z->socc[i]->AAmcc = soccAAmcc[i];
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
        tmp = pbm->P_of_z->h_hat+pbm->P_of_z->nb_lin_constr+{prefix}{dataname}_nb_qc+i;
        {prefix}_mtx_multiply_mtx_vec(tmp, pbm->P_of_z->socc[i]->b, pbm->P_of_z->socc[i]->b,
                                    1, pbm->P_of_z->socc[i]->rowsA);
        tmp[0] = pbm->P_of_z->socc[i]->d[0]*pbm->P_of_z->socc[i]->d[0] - tmp[0];
    }}
    /*TODO TODOTODO*/
#endif    
    
    pbm->r_d = {prefix}{dataname}_r_d;
    pbm->r_d_soft = {prefix}{dataname}_r_d_soft;
    pbm->r_p = {prefix}{dataname}_r_p;
    pbm->st_size = {prefix}{dataname}_st_size;
    
    pbm->tmp1_optvar_seqlen = {prefix}{dataname}_tmp1_optvar_seqlen;
    pbm->tmp2_optvar_seqlen = {prefix}{dataname}_tmp2_optvar_seqlen;
    pbm->tmp2_dual_seqlen = {prefix}{dataname}_tmp2_dual_seqlen;
    pbm->tmp3_mtx_optvar_nb_of_ueq = {prefix}{dataname}_tmp3_mtx_optvar_nb_of_ueq;
    pbm->tmp3_mtx_optvar_nb_of_soft = {prefix}{dataname}_tmp3_mtx_optvar_nb_of_soft;
    pbm->tmp3_state_veclen = {prefix}{dataname}_tmp3_state_veclen;
    pbm->tmp4_mtx_optvar_optvar = {prefix}{dataname}_tmp4_mtx_optvar_optvar;
    pbm->tmp5_nb_of_constr = {prefix}{dataname}_tmp5_nb_of_constr;
    pbm->tmp6_optvar_seqlen = {prefix}{dataname}_tmp6_optvar_seqlen;
    pbm->tmp7_dual_seqlen = {prefix}{dataname}_tmp7_dual_seqlen;
    pbm->tmp8_L_Y = {prefix}{dataname}_tmp8_L_Y;
    pbm->tmp9_L_Y_T = {prefix}{dataname}_tmp9_L_Y_T;
    pbm->tmp8_L_Phi = {prefix}{dataname}_tmp8_L_Phi;
    pbm->tmp9_L_Phi_T = {prefix}{dataname}_tmp9_L_Phi_T;
    pbm->tmp_phibl1 = {prefix}{dataname}_tmp_phibl1;
    pbm->tmp_phibl2 = {prefix}{dataname}_tmp_phibl2;
    pbm->tmp_phibl3 = {prefix}{dataname}_tmp_phibl3;
    pbm->tmp10 = {prefix}{dataname}_tmp10;
    pbm->tmpYbl = {prefix}{dataname}_tmpYbl;
    pbm->tmpQbl = {prefix}{dataname}_tmpQbl;
    pbm->tmp1_res_os = {prefix}{dataname}_tmp1_res_os;
    pbm->tmp2_res_os = {prefix}{dataname}_tmp2_res_os;
    pbm->tmp3_res_ds = {prefix}{dataname}_tmp3_res_ds;
    
    pbm->eye_optvar_veclen = {prefix}{dataname}_eye_optvar_veclen;
    for (i = 0; i < pbm->optvar_veclen; i++){{
        for (j = 0; j < pbm->optvar_veclen; j++)
            pbm->eye_optvar_veclen[i*pbm->optvar_veclen+j] = (i == j) ? 1. : 0.;
    }}
    /*eye(pbm->eye_optvar_veclen, pbm->optvar_veclen);*/
    pbm->eye_state_veclen = {prefix}{dataname}_eye_state_veclen;
    for (i = 0; i < pbm->state_veclen; i++){{
        for (j = 0; j < pbm->state_veclen; j++)
            pbm->eye_state_veclen[i*pbm->state_veclen+j] = (i == j) ? 1. : 0.;
    }}
    return;
}}