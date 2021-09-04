#define {PREFIX}_DEBUG_ALL_MODE_NO
/** Second order cone program solver based on a primal barrier interior point
  * method.
  */ 

#ifdef {PREFIX}_DEBUG_ALL_MODE
#include "stdio.h"
#endif

#include <string.h>  /* TODO For? */

#include "{prefix}mtxops.h"
#include "{prefix}pbm.h"

/* #include <{prefix}usefull.h> */  /* TODO implement {prefix}usefull.h */
#include "{prefix}pbmsolve.h"
#if {PREFIX}_PBM_PRB_COND
#define HHMPC_OS 10 /* TODO */
#else
#include "static_data.h" /* TODO implement static_data.h */
#endif
/* static functions declaration */

static void {prefix}_pbm_get_valid(real_t *z_valid, real_t *tmp_otpseq_optseq,
                                   real_t *tmp_optseq,
                                   const real_t *P, const real_t *P_T, const real_t *h,
                                   const uint32_t optseq, const uint32_t nb_ineq,
                                   real_t *tmp3_optseq, real_t *tmp4_optseq,
                                   real_t *tmp1_optvec_optvec, real_t *tmp2_optvec_optvec,
                                   real_t *tmp5_nb_of_constr);

static void form_d(
                real_t *d, const real_t *P, const real_t *h, const real_t *z,
                const uint32_t rowsP, const uint32_t colsP);

static void form_dsoft(
                real_t *dsoft, real_t *diag_d_soft, real_t *rd_soft,
                real_t *Phi_soft,
                real_t *tmp,
                const real_t *roh, const real_t *z,
                const real_t *Psoft, const real_t *Psoft_T, const real_t *hsoft,
                const real_t *Fusoft, const real_t *Fxsoft, const real_t *Ffsoft,
                const uint32_t rowsFusoft, const uint32_t control_veclen,
                const uint32_t rowsFfsoft, const uint32_t state_veclen,
                const uint32_t horizon, const struct {prefix}_pbm *pbm);

static void form_diag_d_sq(
                real_t *diag_d_sq, const real_t *d, const uint32_t dim);

static void form_Phi(
                real_t *Phi, real_t *help, real_t *tmp_Phi,
                const struct {prefix}_pbm *pbm,
                const real_t *H, const real_t *P_T, const real_t *P,
                const struct {prefix}_pbm_P_hat *P_hat,
                const real_t *d, const real_t *diag_d_sq,
                const real_t kappa,
                const uint32_t optvar_seqlen, const uint32_t nb_of_ueq_constr);

static void {prefix}_residual_norm(
                real_t *f,
                const real_t *r_d, const real_t* r_p,
                const uint32_t optvar_seqlen, const uint32_t dual_seqlen);

static void {prefix}_bt_line_search_new(
                real_t *good_step, const struct {prefix}_pbm *pbm);

static void {prefix}_bt_line_search(
                real_t *good_step, const struct {prefix}_pbm *pbm);

static void {prefix}_pbm_warm_start(
                const struct {prefix}_pbm *pbm);

static void {prefix}_multiply_H_z(
                real_t *product, const real_t *H, const real_t *z,
                const struct {prefix}_pbm *pbm);

static void {prefix}_multiply_P_z(
                real_t *product, const real_t *P, const real_t *z,
                const struct {prefix}_pbm *pbm);

static void {prefix}_multiply_P_T_d(
                real_t *product, const real_t *P_T, const real_t *z,
                const struct {prefix}_pbm *pbm);

static uint32_t {prefix}_better_step_size(
                const struct {prefix}_pbm *pbm, real_t f1, real_t f2);

/* external functions definition */

void {prefix}_pbm_test_get_valid(const struct {prefix}_pbm *pbm)
{{
    {prefix}_pbm_get_valid(pbm->z_opt, pbm->tmp4_mtx_optvar_optvar,
                           pbm->tmp1_optvar_seqlen, pbm->P, pbm->P_T,
                           pbm->h, pbm->optvar_seqlen,
                           pbm->nb_of_ueq_constr,
                           pbm->tmp2_optvar_seqlen, pbm->tmp6_optvar_seqlen,
                           pbm->tmp8_L_Phi, pbm->tmp9_L_Phi_T,
                           pbm->tmp5_nb_of_constr);
}}

void {prefix}_pbm_get_valid(real_t *z_valid, real_t *P_T_P, real_t *P_T_h,
                            const real_t *P, const real_t *P_T, const real_t *h,
                            const uint32_t optseq, const uint32_t nb_ineq,
                            real_t *tmp3_os, real_t *tmp4_os,
                            real_t *tmp1_ov_ov, real_t *tmp2_ov_ov,
                            real_t *tmp5_nb_c)
{{
    uint32_t i;
    /* TODO tmp_os_x_os for sparse date or use block structure. */
    real_t tmp1[HHMPC_OS*HHMPC_OS]; /* = tmp1_ov_ov; TODO */
    real_t tmp2[HHMPC_OS*HHMPC_OS]; /* = tmp2_ov_ov; TODO */
    real_t *tmp3 = tmp3_os;
    real_t *z_delta = tmp4_os;
    real_t *tmp4 = tmp5_nb_c;
    real_t *v = tmp1;  /* Pointer to first entry of tmp1. */
    real_t *tmp5 = tmp1+1;  /* Pointer to second entry of tmp1,
                            * used as array of lenght HHMPC_OS. */
    
    /* Get violation lines of the inequality constraints. */
    {prefix}_mtx_multiply_mtx_vec(tmp4, P, z_valid, nb_ineq, optseq);
    {prefix}_mtx_substract_direct(tmp4, h, nb_ineq, 1);
    for (i = 0; i < nb_ineq; i++){{/*
        printf("%f\n", tmp4[i]);*/
        tmp4[i] = (tmp4[i] < 0.) ? 0. : tmp4[i] + 0.0001;/*
        printf("%f\n", tmp4[i]);*/
    }}
    
    /* */
    {prefix}_mtx_multiply_mtx_vec(P_T_h, P_T, tmp4, optseq, nb_ineq);
    {prefix}_mtx_multiply_mtx_mtx(P_T_P, P_T, P, optseq, nb_ineq, optseq);
    for (i = 0; i < optseq; i++){{
            P_T_P[i*optseq+i] += 0.0001;
    }}
    solveBlock(z_delta, tmp1, tmp2, P_T_P,
                optseq, P_T_h, 1, tmp3);
    {prefix}_mtx_multiply_mtx_vec(tmp5, P, z_delta, nb_ineq, optseq);
    
    v[0] = 0.;
    for (i = 0; i < nb_ineq; i++){{
        v[0] = (tmp4[i]/tmp5[i] > v[0]) ? tmp4[i]/tmp5[i] : v[0];
    }}/*
    printf("v = %f\n", v[0]);*/
    {prefix}_mtx_scale_direct(z_delta, v[0], optseq, 1);
    
    {prefix}_mtx_substract_direct(z_valid, z_delta, optseq, 1);
}}

/* Minimize MPC qp/qcqp/socp using primal barrier interior method */
void {prefix}_pbm_solve_problem(const struct {prefix}_pbm *pbm)
{{
    uint32_t i, j;
    real_t *t_solve_optvar_seqlen = pbm->tmp1_optvar_seqlen;
    real_t *t_optvar_seqlen = pbm->tmp2_optvar_seqlen;
    real_t *t_solve_dual_seqlen = pbm->tmp2_dual_seqlen;
    real_t *t_L_Y = pbm->tmp8_L_Y;
    real_t *t_L_Y_T = pbm->tmp9_L_Y_T;
    real_t *eye_nm = pbm->eye_optvar_veclen;
    real_t *eye_n = pbm->eye_state_veclen;
    real_t f;
#if {PREFIX}_PBM_PRB_COND
    real_t delta_to_zero = 0.01 /* pos */, delta_to_zero_lin = 0.05;
#endif
/*//     printf("g = ");
//     print_mtx(pbm->g, 5, 1);
//     printf("H = ");
//     print_mtx(pbm->H, 5, 5);*/
#if {PREFIX}_PBM_PRB_COND/*
//     for (i = 0; i < 5; i++){{
//         for (j = 0; j < 5; j++){{
//             pbm->H[i*5+j] *= 500.;
//         }}
//     }}*//* TODO Daten korrigieren? */
    for (i = 0; i < HHMPC_OS; i++){{
            pbm->g[i] *= 2.;
    }}
#endif
/*     printf(" inner step %d ", inner); */
    /*Check if initial value is valid*/
/*     printf("%d\n", {prefix}_pbm_check_valid(pbm, pbm->z_ini)); */
    
    /*Take initial value*/
    memcpy(pbm->z_opt, pbm->z_ini, pbm->sizeof_optvar_seqlen);
    memcpy(pbm->v_opt, pbm->v_ini, pbm->sizeof_dual_seqlen);
    
    /* Calculate Kappa in every time_step */
    {prefix}_pbm_calc_kappa(pbm->kappa, pbm, pbm->z_opt);
#ifdef {PREFIX}_DEBUG_ALL_MODE
    printf("calculated kappa = %.20f\n", pbm->kappa[0]);
#endif    
      
    /* Update h for new xk */
    memcpy(pbm->P_of_z->h_hat, pbm->P_of_z->h, sizeof(real_t) * pbm->P_of_z->nb_lin_constr);
#if 0
#ifdef HHMPC_SOCPCONDTEST    
/*     {prefix}_pbm_get_valid_trick(pbm);*/  /* Trick we used before */
#ifdef HHMPC_SOCPCONDTEST5
    {prefix}_pbm_get_positiv(pbm, &delta_to_zero);  /* Only ensure, that c*z +d > 0 */
#endif
#ifndef HHMPC_SOCPCONDTEST5
    {prefix}_pbm_get_valid_lin_constr(pbm, &delta_to_zero_lin);
#endif
#endif
#endif/*

#ifdef HHMPC_SOCPCONDTEST    */
    /*
    {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
    {prefix}_pbm_get_positiv(pbm, &delta_to_zero);
    */
    {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
    {prefix}_pbm_get_valid(pbm->z_opt, pbm->tmp4_mtx_optvar_optvar,
                           pbm->tmp1_optvar_seqlen, pbm->P, pbm->P_T,
                           pbm->h, pbm->optvar_seqlen,
                           pbm->nb_of_ueq_constr,
                           pbm->tmp2_optvar_seqlen, pbm->tmp6_optvar_seqlen,
                           pbm->tmp8_L_Phi, pbm->tmp9_L_Phi_T,
                           pbm->tmp5_nb_of_constr);/*
#endif*/
/*     print_mtx(pbm->h, 1, pbm->nb_of_ueq_constr); */
    /*Improve z for a fixed number of steps j_in*/
    for (i = 0; i < *(pbm->j_in); i++) {{  /*  TODO here error, compiling for MC */
        {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
        /*Check if initial value is valid*/
        if ({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1){{
            printf("%d return, j = %d", {prefix}_pbm_check_valid(pbm, pbm->z_opt), i);
/*//             while ({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1){{
//                 pbm->z_opt[{prefix}_pbm_check_valid(pbm, pbm->z_opt)] *= 1.1;
//                 {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
//                        t_solve_optvar_seqlen, t_optvar_seqlen);
//             }}
//             uint32_t z = {prefix}_pbm_check_valid(pbm, pbm->z_opt);
//             if (pbm->z_opt[z]>0){{
//                 pbm->z_opt[z] = 0.04;
//                 pbm->z_ini[z] = 0.04; 
//             }}
//             if (pbm->z_opt[z]<0){{
//                 pbm->z_opt[z] = 0.04;
//                 pbm->z_ini[z] = 0.04; 
//             }}
//             real_t h_low[5];
//             real_t z_low[5];
//             for (j = 0; j < 5; j++){{
//                 h_low[j] = pbm->h[j+pbm->P_of_z->nb_lin_constr] - 0.001;
//             }}
//             fwd_subst(z_low, pbm->P+5*pbm->P_of_z->nb_lin_constr, 5, h_low, 1);
// //             print_mtx(z_low, 5, 1);
//             printf("%d return", {prefix}_pbm_check_valid(pbm, z_low));
//             HIER
//             return;*/
        }}
/*         printf("%d\n", {prefix}_pbm_check_valid(pbm, pbm->z_opt)); */
        {prefix}_form_pbm_matrices(pbm->d, pbm->diag_d_sq, pbm->dsoft,
                                   pbm->diag_d_soft, pbm->Phi, pbm);
        
        /* Calculate the residual */
        residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
        {prefix}_residual_norm(&f, pbm->r_d, pbm->r_p, pbm->optvar_seqlen, pbm->dual_seqlen);
/*//         print_mtx(pbm->r_p, pbm->dual_seqlen, 1);
//         printf("res_norm = %f\n", f);
//         print_mtx(pbm->Phi, pbm->optvar_seqlen, pbm->optvar_seqlen); */
        /* Solve system of linear equations to obtain the step direction */
        
        /* regularization (+epsilon*I) */
        for (j = 0; j < pbm->optvar_seqlen; j++){{
            pbm->Phi[j*pbm->optvar_seqlen+j] += pbm->reg[0];
        }}
        solve_sysofleq(pbm->delta_z, pbm->delta_v, pbm->Phi, pbm->r_d, pbm->r_p,
                       pbm->C, pbm->C_T, pbm->A, pbm->A_T, pbm->B, pbm->B_T,
                       pbm->A_B, pbm->A_B_T, pbm->H,
                       pbm->state_veclen, pbm->optvar_veclen-pbm->state_veclen, /* TODO m einführen */
                       pbm->horizon,
                       eye_nm, eye_n,
                       t_solve_optvar_seqlen,
                       pbm->tmp2_optvar_seqlen,
                       t_solve_dual_seqlen,
                       t_L_Y, t_L_Y_T,
                       pbm->tmp_phibl1, pbm->tmp_phibl2,
                       pbm->tmp_phibl3,
                       pbm->tmpQbl, pbm->tmpYbl,
                       pbm->tmp10,
                       pbm->tmp8_L_Phi, pbm->tmp9_L_Phi_T);
        /* undo regularization (-epsilon*I) */
        for (j = 0; j < pbm->optvar_seqlen; j++){{
            pbm->Phi[j*pbm->optvar_seqlen+j] -= pbm->reg[0];
        }}

/*         print_mtx(pbm->delta_z, pbm->optvar_seqlen, 1); */
        {prefix}_pbm_iterative_refinement(pbm);

        /* Find best step size (0...1] */
/*         real_t test = 0.; */
         {prefix}_bt_line_search_new(pbm->st_size, pbm);
/*//         test = pbm->st_size[0];
//         printf("new st_size = %.8f\n", pbm->st_size[0]);
//         bt_line_search(pbm->st_size, pbm);
//         if(smpl_abs(test-pbm->st_size[0]) > 0.00){{
//             HIER
//             printf("new st_size = %.8f\n", test);
//             printf("old st_size = %.8f\n", pbm->st_size[0]);
//         }}
//         pbm->st_size[0] = test;
//         printf("old st_size = %.8f\n", pbm->st_size[0]); */
        
        /* Update z */
        {prefix}_mtx_scale_direct(pbm->delta_z, pbm->st_size[0],
                                pbm->optvar_seqlen, 1);
        {prefix}_mtx_add_direct(pbm->z_opt, pbm->delta_z,
                              pbm->optvar_seqlen, 1);
        {prefix}_mtx_scale_direct(pbm->delta_v, pbm->st_size[0],
                                pbm->dual_seqlen, 1);
        {prefix}_mtx_add_direct(pbm->v_opt, pbm->delta_v,
                              pbm->dual_seqlen, 1);
        {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
           t_solve_optvar_seqlen, t_optvar_seqlen);
        if ({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1){{
            printf("No valid step possible: return");
            return;
        }}
/*//         print_mtx(pbm->delta_z, pbm->optvar_seqlen, 1);
//         print_mtx(pbm->z_opt, pbm->optvar_veclen, 1);
//         print_mtx(pbm->v_opt, pbm->dual_seqlen, 1); */
        if (f <= 1e-12){{
/*//             printf("break, res_norm = %f\n", f);
//             break; */
        }}
/*         pbm->kappa[0] *= 0.8; */
    }}
/*     pbm->kappa[0] = 95.; */
    {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
           t_solve_optvar_seqlen, t_optvar_seqlen);
    form_d(pbm->d, pbm->P, pbm->h, pbm->z_opt,
           pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
    form_dsoft(pbm->dsoft, pbm->diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
                   pbm->tmp3_mtx_optvar_nb_of_soft,
                   pbm->roh, pbm->z_opt, 
                   pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
                   pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
                   pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
                   pbm->state_veclen, pbm->horizon, pbm);
    residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
    {prefix}_residual_norm(&f, pbm->r_d, pbm->r_p, pbm->optvar_seqlen, pbm->dual_seqlen);
   /* printf("res_norm = %.11f\n", f);*/
    /* Update x_k (und andere Parameter) */
    if (pbm->conf->warm_start) {{
        {prefix}_pbm_warm_start(pbm);
    }}
/*     print_mtx(pbm->z_opt, pbm->optvar_veclen, 1);
//     print_mtx(pbm->r_p, pbm->state_veclen, 1); */
    /*
    memcpy(pbm->z_ini, pbm->z_opt, pbm->sizeof_optvar_seqlen);
    memcpy(pbm->v_ini, pbm->v_opt, pbm->sizeof_dual_seqlen);*/
}}

/* TODO Describe external function definitions */
uint32_t {prefix}_pbm_check_positiv(
                const struct {prefix}_pbm *pbm, const real_t *z_check)
{{
    uint32_t i;
    real_t *tmp_os = pbm->tmp1_optvar_seqlen;
    
    for (i = 0; i < pbm->P_of_z->nb_socc; i++){{
        {prefix}_mtx_multiply_mtx_vec(tmp_os, pbm->P_of_z->socc[i]->c, z_check,
                                    1, pbm->optvar_seqlen);
        if (tmp_os[0] + pbm->P_of_z->socc[i]->d[0] < 0){{
            /*
            printf("c = %f %f, d = %f \n", pbm->P_of_z->socc[i]->c[0], pbm->P_of_z->socc[i]->c[1], pbm->P_of_z->socc[i]->d[0]);
            printf("A = %f %f %f %f, b = %f %f\n", pbm->P_of_z->socc[i]->A[0], pbm->P_of_z->socc[i]->A[1],pbm->P_of_z->socc[i]->A[2], pbm->P_of_z->socc[i]->A[3], pbm->P_of_z->socc[i]->b[0], pbm->P_of_z->socc[i]->b[1]);
            */
            return i;    
        }}
    }}
    
    return -1;
}}

uint32_t {prefix}_pbm_check_valid(
                const struct {prefix}_pbm *pbm, const real_t *z_check)
{{
    uint32_t i;
    real_t *help1 = pbm->tmp5_nb_of_constr;
    
    {prefix}_multiply_P_z(help1, pbm->P, z_check, pbm);

    {prefix}_mtx_substract_direct(help1, pbm->h, pbm->nb_of_ueq_constr, 1);
/*//     print_mtx(help1, pbm->nb_of_ueq_constr, 1);
//     printf("check:\n");*/
#ifdef HHMPC_QPTEST
    for (i = 0; i < 2; i++){{
/*         printf("%f\n", help1[i]); */
        if (help1[i] > 0) {{return i;}}
    }}
    for (i = 3; i < pbm->nb_of_ueq_constr; i++){{
/*         printf("%f\n", help1[i]); */
        if (help1[i] > 0) {{return i;}}
    }}
#endif
#ifdef HHMPC_QPSMALLTEST
    for (i = 0; i < 2; i++){{
/*         printf("%f\n", help1[i]); */
        if (help1[i] >= 0.) {{return i;}}
    }}
/* //         printf("%f\n", help1[i]);
//     if (help1[i] > 0.) {{
//         print_mtx(help1, pbm->nb_of_ueq_constr, 1);
//         return i;
//     }} */
    for (i = 3; i < pbm->nb_of_ueq_constr; i++){{
/*         printf("%f\n", help1[i]); */
        if (help1[i] >= 0.) {{return i;}}
    }}
#endif
#if {PREFIX}_PBM_PRB_COND

    /*
//     
//     print_mtx(pbm->P, pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
//     print_mtx(pbm->h, 1, pbm->nb_of_ueq_constr);
//     
//     print_mtx(pbm->z_opt, 1, pbm->optvar_seqlen);
// 
//     real_t help01[pbm->P_of_z->nb_socc];
//     real_t help02[pbm->P_of_z->nb_socc];
//     for (i = 0; i < pbm->P_of_z->nb_socc; i++){{
// //         printf("%d\n", i); 
//         {prefix}_mtx_multiply_mtx_vec(help01, pbm->P_of_z->socc[i]->A, pbm->z_opt, pbm->P_of_z->nb_socc, 5);
//         {prefix}_mtx_add_direct(help01, pbm->P_of_z->socc[i]->b, pbm->P_of_z->nb_socc, 1);
//         {prefix}_mtx_multiply_mtx_vec(help02, help01, help01, 1, pbm->P_of_z->nb_socc);
// //         printf("Au+b = %f", help02[0]);
// //         print_mtx(help01, 5, 1);
//         {prefix}_mtx_multiply_mtx_vec(help01, pbm->P_of_z->socc[i]->c, pbm->z_opt, 1, pbm->P_of_z->nb_socc);
//         {prefix}_mtx_add_direct(help01, pbm->P_of_z->socc[i]->d, pbm->P_of_z->nb_socc, 1);
// //         printf("cu+d = %f \n", help01[0]);*/
/*         if (help01[0] < help02[0]) {{*//*print_mtx(help1, pbm->nb_of_ueq_constr, 1); printf("%d\n", i);*/ /*return i;}}*/
/*// //         print_mtx(help01, 5, 1);
//     }} */
    if ({prefix}_pbm_check_positiv(pbm, z_check) + 1){{
        return pbm->nb_of_ueq_constr;
    }}else{{
        for (i = 0; i < pbm->nb_of_ueq_constr; i++){{/*
            printf("%f\n", help1[i]);*/
            if (help1[i] >= 0.0) {{return i;}}
        }}
    }}
#endif
#ifdef HHMPC_SOCPTEST
    for (i = 0; i < pbm->nb_of_ueq_constr; i++){{
//         printf("%f\n", help1[i]);
        if (help1[i] >= 0.) {{return i;}}
    }}
#endif
#ifdef HHMPC_SOCPSOFTTEST
    for (i = 0; i < pbm->nb_of_ueq_constr; i++){{
//         printf("%f\n", help1[i]);
        if (help1[i] >= 0.) {{return i;}}
    }}
#endif
#ifdef HHMPC_SOCPHARDTEST
    for (i = 0; i < 2; i++){{
//         printf("%f\n", help1[i]);
        if (help1[i] > 0.) {{return i;}}
    }}
    for (i = 3; i < pbm->nb_of_ueq_constr; i++){{
//         printf("%f\n", help1[i]);
        if (help1[i] > 0.) {{return i;}}
    }}
#endif
    return -1;
}}

/* Ensure, that c*z +d > 0. So we can use generalized inequalities*/
void {prefix}_pbm_get_positiv(const struct {prefix}_pbm *pbm, real_t *delta)
{{
    real_t *t_solve_optvar_seqlen = pbm->tmp1_optvar_seqlen;
    real_t *t_optvar_seqlen = pbm->tmp2_optvar_seqlen;
    uint32_t i, k;
    {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
        t_solve_optvar_seqlen, t_optvar_seqlen);
    
    if ({prefix}_pbm_check_positiv(pbm, pbm->z_opt)+1){{
/*//         print_mtx(pbm->z_opt, 5, 1);
//         printf("corrected pos z_opt[0] %d\n", {prefix}_pbm_check_positiv(pbm, pbm->z_opt));*/
        for (k = 0; k < pbm->P_of_z->nb_socc; k++){{
            if (k == {prefix}_pbm_check_positiv(pbm, pbm->z_opt)){{
                pbm->z_opt[k] = delta[0] - pbm->P_of_z->socc[k]->d[0];
                for (i = 0; i < k; i++){{
                    pbm->z_opt[k] -= pbm->P_of_z->socc[k]->c[i]*pbm->z_opt[i];
                }}
                pbm->z_opt[k] /= pbm->P_of_z->socc[k]->c[k];
            }}
        }}
/*//         print_mtx(pbm->z_opt, 5, 1);
// //         if ({prefix}_pbm_check_positiv(pbm, pbm->z_opt)+1){{
// //             printf("not corrected pos z_opt[0] enough %d\n", {prefix}_pbm_check_positiv(pbm, pbm->z_opt));
// //             return;
// //         }}*/
    }}
}}

/* Ensure, that P*z > h */
void {prefix}_pbm_get_valid_lin_constr(
                const struct {prefix}_pbm *pbm, real_t *delta)
{{
    /* P für condensiertes Problem hat hinsichtlich der Zeilen, die für die
     * Einhaltung der UNB für die Zustände verantwortlich sind Dreiecksgestalt 
     * Deswegen hier jede zweite Zeile der zweiten Hälfte von P (upper bounds)
     * verwenden, um sicherzustellen, dass z valide ist */
    
    real_t *t_solve_optvar_seqlen = pbm->tmp1_optvar_seqlen;
    real_t *t_optvar_seqlen = pbm->tmp2_optvar_seqlen;
    uint32_t i, k;
    {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
        t_solve_optvar_seqlen, t_optvar_seqlen);
/*//     printf("1corrected pos z_opt[0] %d\n", {prefix}_pbm_check_valid(pbm, pbm->z_opt));
//     pbm->z_opt[0] = pbm->P_of_z->h[60] - delta[0];
//     pbm->z_opt[0] /= pbm->P_of_z->P[60*HHMPC_OS];
//     printf("2corrected pos z_opt[0] %d\n", {prefix}_pbm_check_valid(pbm, pbm->z_opt));
//     delta[0] = (pbm->P_of_z->h[pbm->P_of_z->nb_lin_constr/2] + pbm->P_of_z->h[0])/2.;*/

    if ({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1){{
/*         print_mtx(pbm->z_opt, 30, 1); */
        printf("corrected valid z_opt[0] %d\n", {prefix}_pbm_check_valid(pbm, pbm->z_opt));
        for (k = 0; k < 1; k++){{
/*             if (pbm->P_of_z->nb_lin_constr/2 + 2*k == {prefix}_pbm_check_valid(pbm, pbm->z_opt)){{ */
                pbm->z_opt[k] = pbm->P_of_z->h[pbm->P_of_z->nb_lin_constr/2 + 2*k] - delta[0];
                for (i = 0; i < k; i = i+2){{
                    pbm->z_opt[k] -= pbm->P_of_z->P[(pbm->P_of_z->nb_lin_constr/2 + 2*k)*HHMPC_OS + i]*pbm->z_opt[i];
                }}
                pbm->z_opt[k] /= pbm->P_of_z->P[(pbm->P_of_z->nb_lin_constr/2 + 2*k)*HHMPC_OS + k];
/*             }} */
        }}
    }}
    if ({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1){{
/*         print_mtx(pbm->z_opt, 30, 1); */
        printf("corrected valid z_opt[0] %d\n", {prefix}_pbm_check_valid(pbm, pbm->z_opt));
        for (k = 0; k < pbm->optvar_seqlen; k++){{
/*             if ((pbm->P_of_z->nb_lin_constr/2 + 2*k) == {prefix}_pbm_check_valid(pbm, pbm->z_opt)){{ */
                pbm->z_opt[k] = pbm->P_of_z->h[pbm->P_of_z->nb_lin_constr/2 + 2*k] - delta[0];
                for (i = 0; i < k; i = i+2){{
                    pbm->z_opt[k] -= pbm->P_of_z->P[(pbm->P_of_z->nb_lin_constr/2 + 2*k)*HHMPC_OS + i]*pbm->z_opt[i];
                }}
                pbm->z_opt[k] /= pbm->P_of_z->P[(pbm->P_of_z->nb_lin_constr/2 + 2*k)*HHMPC_OS + k];
/*             }} */
        }}
    
/*//         print_mtx(pbm->z_opt, 30, 1);
// //         if ({prefix}_pbm_check_positiv(pbm, pbm->z_opt)+1){{
// //             printf("not corrected pos z_opt[0] enough %d\n", {prefix}_pbm_check_positiv(pbm, pbm->z_opt));
// //             return;
// //         }} */
    
    }}
}}

void {prefix}_pbm_get_valid_trick(const struct {prefix}_pbm* pbm)
{{
    real_t *t_solve_optvar_seqlen = pbm->tmp1_optvar_seqlen;
    real_t *t_optvar_seqlen = pbm->tmp2_optvar_seqlen;
    
    
/*// //     i = 1;
// //     {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
// //             t_solve_optvar_seqlen, t_optvar_seqlen);
// //     printf("%d\n", {prefix}_pbm_check_valid(pbm, pbm->z_opt)+1);
// //     while ({prefix}_pbm_check_valid(pbm, pbm->z_opt) >= pbm->P_of_z->nb_lin_constr){{
// // //         if (pbm->P_of_z->socc[0]->c[0] * pbm->z_opt[0] <
// // //                 (0.002*i*i - pbm->P_of_z->socc[0]->d[0])){{
// //             printf("corrected %d z_opt[0]\n", i);
// //             pbm->z_opt[0] = 
// //                 (0.001*i*i - pbm->P_of_z->socc[0]->d[0]) / pbm->P_of_z->socc[0]->c[0];
// //             print_mtx(pbm->z_opt, 5, 1);
// // //         }}
// //         {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
// //             t_solve_optvar_seqlen, t_optvar_seqlen);
// //         i++;
// //     }} */
    
        {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
            t_solve_optvar_seqlen, t_optvar_seqlen);
        if ({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1){{
        
/*//             if ({prefix}_pbm_check_positiv(pbm, pbm->z_opt)+1){{
//                 print_mtx(pbm->z_opt, 5, 1);
//                 printf("corrected pos z_opt[0]\n");
//                 pbm->z_opt[0] = (-1 * pbm->P_of_z->socc[0]->d[0]) / pbm->P_of_z->socc[0]->c[0];
//                 print_mtx(pbm->z_opt, 5, 1);
//             }}
//             real_t tmp[5];
//             print_mtx(pbm->P_of_z->socc[4]->A, 5, 5);
//             print_mtx(pbm->P_of_z->socc[4]->b, 1, 5);
//             {prefix}_mtx_multiply_mtx_vec(tmp, pbm->P_of_z->socc[4]->A, pbm->z_opt, 5, 5);
//             print_mtx(tmp, 5, 1);
//             {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
//                t_solve_optvar_seqlen, t_optvar_seqlen);
// //             print_mtx(pbm->P, pbm->nb_of_ueq_constr, 5);
//             printf("%f\n", pbm->P[pbm->P_of_z->nb_lin_constr*pbm->optvar_seqlen]);
//             pbm->z_opt[0] = (pbm->h[pbm->P_of_z->nb_lin_constr] - 0.0001) /
//                             pbm->P[pbm->P_of_z->nb_lin_constr*pbm->optvar_seqlen];
//             print_mtx(pbm->z_opt, 5, 1); */
                            
            if (pbm->P_of_z->socc[0]->c[0] * pbm->z_opt[0] <
                    (0.01 - pbm->P_of_z->socc[0]->d[0])){{
                /* printf("corrected 1 z_opt[0]\n"); */
                pbm->z_opt[0] = 
                    (0.01 - pbm->P_of_z->socc[0]->d[0]) / pbm->P_of_z->socc[0]->c[0];
                /* print_mtx(pbm->z_opt, 5, 1); */           
            }}
        }}
        {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
            t_solve_optvar_seqlen, t_optvar_seqlen);
        if ({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1){{
            if (pbm->P_of_z->socc[0]->c[0] * pbm->z_opt[0] <
                    (0.05 - pbm->P_of_z->socc[0]->d[0])){{
                /* printf("corrected 2 z_opt[0]\n"); */
                pbm->z_opt[0] = 
                    (0.05 - pbm->P_of_z->socc[0]->d[0]) / pbm->P_of_z->socc[0]->c[0];
                /* print_mtx(pbm->z_opt, 5, 1); */            
            }}
        }}
        {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
            t_solve_optvar_seqlen, t_optvar_seqlen);
        if ({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1){{
            if (pbm->P_of_z->socc[0]->c[0] * pbm->z_opt[0] <
                    (0.2 - pbm->P_of_z->socc[0]->d[0])){{
                printf("corrected 2 z_opt[0]\n");
                pbm->z_opt[0] = 
                    (0.2 - pbm->P_of_z->socc[0]->d[0]) / pbm->P_of_z->socc[0]->c[0];
                /* print_mtx(pbm->z_opt, 5, 1); */           
            }}
        }}
}}

void {prefix}_pbm_iterative_refinement(const struct {prefix}_pbm *pbm)
{{
#if {PREFIX}_PBM_PRB_SPARSE
    real_t *delta_rd = pbm->tmp1_optvar_seqlen;
    real_t *delta_rp = pbm->tmp2_dual_seqlen;
    real_t *delta_delta_v = pbm->tmp7_dual_seqlen;
    real_t *delta_delta_z = pbm->tmp4_mtx_optvar_optvar+pbm->dual_seqlen;
    real_t *tmp1 = pbm->tmp2_optvar_seqlen;
    real_t *tmp2 = pbm->tmp6_optvar_seqlen;
    real_t *tmp3 = pbm->tmp3_state_veclen;
    real_t *tmp4 = pbm->tmp4_mtx_optvar_optvar;
    
    real_t *L_Phi_blocks = pbm->tmp8_L_Phi;
    real_t *L_Phi_T_blocks = pbm->tmp9_L_Phi_T;
    real_t *L_Y = pbm->tmp8_L_Y;
    real_t *L_Y_T = pbm->tmp9_L_Y_T;

    {prefix}_multiply_C_T_v(tmp1, pbm->C_T, pbm->delta_v,
                            pbm->state_veclen, pbm->control_veclen, pbm->horizon);

    {prefix}_multiply_H_z(delta_rd, pbm->Phi, pbm->delta_z, pbm);

    {prefix}_mtx_add_direct(delta_rd, tmp1, pbm->optvar_seqlen, 1);
    {prefix}_mtx_add_direct(delta_rd, pbm->r_d, pbm->optvar_seqlen, 1);
    {prefix}_mtx_scale_direct(delta_rd, -1., pbm->optvar_seqlen, 1);

    {prefix}_multiply_C_z(delta_rp, pbm->C, pbm->delta_z,
                          pbm->state_veclen, pbm->control_veclen, pbm->horizon);

    {prefix}_mtx_add_direct(delta_rp, pbm->r_p, pbm->dual_seqlen, 1);
    {prefix}_mtx_scale_direct(delta_rp, -1., pbm->dual_seqlen, 1);
   
/*     print_mtx(delta_rp, pbm->dual_seqlen, 1); */
        
    form_beta(delta_delta_v, L_Phi_blocks, L_Phi_T_blocks, delta_rd, delta_rp,
              pbm->horizon, pbm->C,
              pbm->state_veclen, pbm->optvar_veclen-pbm->state_veclen,
              tmp1, tmp2);
    form_delta_v(delta_delta_v, tmp4, tmp3, 
                 L_Y, L_Y_T, pbm->horizon, pbm->state_veclen);
    form_delta_z(delta_delta_z, tmp1, delta_delta_v,
                 L_Phi_blocks, L_Phi_T_blocks, delta_rd, pbm->C_T, pbm->horizon,
                 pbm->state_veclen, pbm->optvar_veclen-pbm->state_veclen);
    
    
/*     print_mtx(delta_delta_z, pbm->optvar_seqlen, 1); */
    {prefix}_mtx_substract_direct(pbm->delta_z, delta_delta_z, pbm->optvar_seqlen, 1);
    {prefix}_mtx_substract_direct(pbm->delta_v, delta_delta_v, pbm->dual_seqlen, 1);
/*     print_mtx(pbm->delta_z, pbm->optvar_seqlen, 1); */
#endif
}}

void {prefix}_pbm_update_P(struct {prefix}_pbm_P_hat *P,
                           const uint32_t optvar_seqlen,
                           real_t *tmp1, real_t *tmp2)
{{
    uint32_t i;
    /* Update matrix P(z) */
    struct {prefix}_pbm_qc *qc_i;
    struct {prefix}_pbm_socc *socc_i;
    /*  TODO P vorher 0 setzen */
    memcpy(P->P_hat, P->P, sizeof(real_t) * P->nb_lin_constr*optvar_seqlen);
    P->P_hat += P->nb_lin_constr*optvar_seqlen;  /* Pointer wird nicht nur lokal verändern */
    
    /* Determine rows for qc */
    for (i = 0; i < P->nb_qc; i++){{
        qc_i = P->qc[i];
        {prefix}_mtx_multiply_mtx_mtx(P->P_hat+qc_i->par_0,
                                    qc_i->par, qc_i->Gamma,
                                    1, qc_i->dimGamma, qc_i->dimGamma);
        {prefix}_mtx_add_direct(P->P_hat+qc_i->par_0,
                              qc_i->beta, 1, qc_i->dimGamma);
        P->P_hat += optvar_seqlen;
    }}
    /* Determine rows for socc */
    for (i = 0; i < P->nb_socc; i++){{
        socc_i = P->socc[i];
        {prefix}_mtx_multiply_mtx_vec(tmp1, socc_i->A, socc_i->par,
                                    socc_i->rowsA, socc_i->colsA);
        {prefix}_mtx_add_direct(tmp1, socc_i->b, 1, socc_i->rowsA);
        {prefix}_mtx_add_direct(tmp1, socc_i->b, 1, socc_i->rowsA);

        {prefix}_mtx_multiply_mtx_mtx(P->P_hat+socc_i->par_0,
                                    tmp1, socc_i->A,
                                    1, socc_i->rowsA, socc_i->colsA);
/*         print_mtx(P->P_hat+socc_i->par_0, 1, socc_i->colsA); */
        {prefix}_mtx_multiply_mtx_vec(tmp1, socc_i->c, socc_i->par,
                                    1, socc_i->colsA);
        tmp1[0] += 2*socc_i->d[0];
        {prefix}_mtx_scale(tmp2, socc_i->c, tmp1[0], socc_i->colsA, 1);
        {prefix}_mtx_substract_direct(P->P_hat+socc_i->par_0,
                                    tmp2, 1, socc_i->colsA);
        P->P_hat += optvar_seqlen;
    }}

    P->P_hat -= (P->nb_lin_constr+P->nb_qc + P->nb_socc)*optvar_seqlen;
    {prefix}_mtx_transpose(P->P_hat_T, P->P_hat,
                         P->nb_lin_constr + P->nb_qc + P->nb_socc,
                         optvar_seqlen);
    
/**************************/
    /* Update matrix P(2*z) */
    memcpy(P->P2_hat, P->P, sizeof(real_t) * P->nb_lin_constr*optvar_seqlen);
    P->P2_hat += P->nb_lin_constr*optvar_seqlen;  /* Pointer wird nicht nur lokal verändern */
    
    /* Determine rows for qc */
    for (i = 0; i < P->nb_qc; i++){{
        qc_i = P->qc[i];
        {prefix}_mtx_multiply_mtx_mtx(P->P2_hat+qc_i->par_0,
                                    qc_i->par, qc_i->Gamma,
                                    1, qc_i->dimGamma, qc_i->dimGamma);
        /* 2*z */
        {prefix}_mtx_scale_direct(P->P2_hat+qc_i->par_0, 2., 1, qc_i->dimGamma);
        {prefix}_mtx_add_direct(P->P2_hat+qc_i->par_0,
                              qc_i->beta, 1, qc_i->dimGamma);
        P->P2_hat += optvar_seqlen;
    }}
    /* Determine rows for socc */
    for (i = 0; i < P->nb_socc; i++){{
        socc_i = P->socc[i];
        {prefix}_mtx_multiply_mtx_vec(tmp1, socc_i->A, socc_i->par,
                                    socc_i->rowsA, socc_i->colsA);
        /* 2*z */
        {prefix}_mtx_scale_direct(tmp1, 2., 1, socc_i->rowsA);
        {prefix}_mtx_add_direct(tmp1, socc_i->b, 1, socc_i->rowsA);
        {prefix}_mtx_add_direct(tmp1, socc_i->b, 1, socc_i->rowsA);
        {prefix}_mtx_multiply_mtx_mtx(P->P2_hat+socc_i->par_0,
                                    tmp1, socc_i->A,
                                    1, socc_i->rowsA, socc_i->colsA);
        {prefix}_mtx_multiply_mtx_vec(tmp1, socc_i->c, socc_i->par,
                                    1, socc_i->colsA);
        /* 2*z */
        tmp1[0] *= 2.;
        tmp1[0] += 2*socc_i->d[0];
        {prefix}_mtx_scale(tmp2, socc_i->c, tmp1[0], socc_i->colsA, 1);
        {prefix}_mtx_substract_direct(P->P2_hat+socc_i->par_0,
                                    tmp2, 1, socc_i->colsA);
        P->P2_hat += optvar_seqlen;
    }}
    
    P->P2_hat -= (P->nb_lin_constr+P->nb_qc + P->nb_socc)*optvar_seqlen;
    {prefix}_mtx_transpose(P->P2_hat_T, P->P2_hat,
                         P->nb_lin_constr + P->nb_qc + P->nb_socc,
                         optvar_seqlen);
}}

/* Calculate the residum r = (r_d, r_p) */
void residual(const struct {prefix}_pbm *pbm,
              const real_t *z, const real_t *v, const real_t *d,
              const real_t kappa)
{{
    real_t *tmp1_os = pbm->tmp1_res_os;
    real_t *tmp2_os = pbm->tmp2_res_os;
    real_t *tmp3_ds = pbm->tmp3_res_ds;
    /* Term kappa*P^T*d */  
    {prefix}_multiply_P_T_d(tmp1_os, pbm->P2_T, d, pbm);

    {prefix}_mtx_scale(pbm->r_d, tmp1_os, kappa, pbm->optvar_seqlen, 1);
    /* Add term C^T*v */
    {prefix}_multiply_C_T_v(tmp1_os, pbm->C_T, v,
                            pbm->state_veclen, pbm->control_veclen, pbm->horizon);

    {prefix}_mtx_add_direct(pbm->r_d, tmp1_os, pbm->optvar_seqlen, 1);
    /* Add term g */
    {prefix}_mtx_add_direct(pbm->r_d, pbm->g,
                          pbm->optvar_seqlen, 1);
    /* Add term 2*H*(z-z_ref) */
    /* Calculate distance to the reference (z-z_ref) */
    {prefix}_mtx_substract(tmp2_os, z, pbm->zref, pbm->optvar_seqlen, 1);

    {prefix}_multiply_H_z(tmp1_os, pbm->H, tmp2_os, pbm);
    
    {prefix}_mtx_scale_direct(tmp1_os, 2., pbm->optvar_seqlen, 1);
    {prefix}_mtx_add_direct(pbm->r_d, tmp1_os, pbm->optvar_seqlen, 1);
    
    {prefix}_mtx_add_direct(pbm->r_d, pbm->r_d_soft, pbm->optvar_seqlen, 1);

    {prefix}_multiply_C_z(pbm->r_p, pbm->C, z,
                          pbm->state_veclen, pbm->control_veclen, pbm->horizon);

    {prefix}_mtx_substract_direct(pbm->r_p, pbm->b, pbm->dual_seqlen, 1);
}}

void {prefix}_form_pbm_matrices(real_t *d, real_t *diag_d_sq,
                                real_t *dsoft, real_t *diag_d_soft,
                                real_t *Phi,
                                const struct {prefix}_pbm *pbm)
{{
    form_d(d, pbm->P, pbm->h, pbm->z_opt,
           pbm->nb_of_ueq_constr, pbm->optvar_seqlen);

    form_dsoft(dsoft, diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
               pbm->tmp3_mtx_optvar_nb_of_soft,
               pbm->roh, pbm->z_opt, 
               pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
               pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
               pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
               pbm->state_veclen, pbm->horizon, pbm);
    
    form_diag_d_sq(diag_d_sq, d, pbm->nb_of_ueq_constr);
    
    form_Phi(Phi, pbm->tmp3_mtx_optvar_nb_of_ueq, pbm->tmp4_mtx_optvar_optvar,
             pbm, pbm->H,
             pbm->P2_T, pbm->P2, pbm->P_of_z , pbm->d, pbm->diag_d_sq,
             pbm->kappa[0], pbm->optvar_seqlen, pbm->nb_of_ueq_constr);
}}

void form_Phi(real_t *Phi, real_t *help, real_t *t_Phi,
              const struct {prefix}_pbm *pbm,
              const real_t *H, const real_t *P_T, const real_t *P,
              const struct {prefix}_pbm_P_hat *P_hat,
              const real_t *d, const real_t *diag_d_sq,
              const real_t kappa,
              const uint32_t optvar, const uint32_t nb_of_ueq)
{{    
    uint32_t i, j, k, pos_d, pos_Phi;
    struct {prefix}_pbm_qc *qc_i;
    struct {prefix}_pbm_socc *socc_i;

#if {PREFIX}_PBM_PRB_SPARSE
    /* Der erste Teilschritt sollte auch bei condensed Problemen funktionieren */
    /* Multiply P_T mit diag_d_sq */
    for (j = 0; j < optvar; j++){{
        for (i = 0; i < nb_of_ueq; i++){{
            help[j*nb_of_ueq+i] = P_T[j*nb_of_ueq+i]*diag_d_sq[i*nb_of_ueq+i];
        }}
    }}
    /* Fast multication of above product with P by use of its structure */
/*//     for (i = 0; i < pbm->control_veclen; i++){{  
//         for (k = 0; k < pbm->control_veclen; k++){{ 
//             {prefix}_mtx_multiply_mtx_vec(Phi+i*optvar+k, help+i*nb_of_ueq, P_T+k*nb_of_ueq, 1, nb_of_ueq);
//         }}
//     }}
//     for (j = 0; j < pbm->horizon-1; j++){{
//         for (i = pbm->control_veclen + j*pbm->optvar_veclen; i < pbm->control_veclen + (j+1)*pbm->optvar_veclen; i++){{ 
//             for (k = pbm->control_veclen + j*pbm->optvar_veclen; k < pbm->control_veclen + (j+1)*pbm->optvar_veclen; k++){{ /
//                 {prefix}_mtx_multiply_mtx_vec(Phi+i*optvar+k, help+i*nb_of_ueq, P_T+k*nb_of_ueq, 1, nb_of_ueq);
//             }}
//         }}
//     }}
//     for (i = pbm->control_veclen + j*pbm->optvar_veclen; i < (j+1)*pbm->optvar_veclen; i++){{  
//         for (k = pbm->control_veclen + j*pbm->optvar_veclen; k < (j+1)*pbm->optvar_veclen; k++){{ 
//             {prefix}_mtx_multiply_mtx_vec(Phi+i*optvar+k, help+i*nb_of_ueq, P_T+k*nb_of_ueq, 1, nb_of_ueq);
//         }}
//     }} */
    
    
    for (i = 0; i < pbm->optvar_seqlen; i++){{
        {prefix}_multiply_P_T_d(Phi+i*optvar, P_T, help+i*nb_of_ueq, pbm);
    }}
    
    
    
    
#endif
#if {PREFIX}_PBM_PRB_COND
    {prefix}_mtx_multiply_mtx_mtx(help, P_T, diag_d_sq, optvar, nb_of_ueq, nb_of_ueq);
    {prefix}_mtx_multiply_mtx_mtx(Phi, help, P, optvar, nb_of_ueq, optvar);
#endif
    /* Additional terms for Phi resulting of the second derivative of qc */
    for (i = 0; i < P_hat->nb_qc; i++){{
        qc_i = P_hat->qc[i];
        pos_d = P_hat->nb_lin_constr+i;
        {prefix}_mtx_scale(t_Phi, qc_i->Gamma, 2*d[pos_d],
                         qc_i->dimGamma, qc_i->dimGamma);
        
        for (j = 0; j < qc_i->dimGamma; j++){{
            pos_Phi = (qc_i->par_0 + j)*optvar + qc_i->par_0;
            {prefix}_mtx_add_direct(Phi+pos_Phi, t_Phi+j*qc_i->dimGamma, 1, qc_i->dimGamma);
        }}
    }}
    /* Additional terms for Phi resulting of the second derivative of socc */
    for (i = 0; i < P_hat->nb_socc; i++){{
        socc_i = P_hat->socc[i];
        pos_d = P_hat->nb_lin_constr + P_hat->nb_qc + i;
        {prefix}_mtx_scale(t_Phi, socc_i->AAmcc, 2*d[pos_d],
                         socc_i->colsA, socc_i->colsA);
        
        for (j = 0; j < socc_i->colsA; j++){{
            pos_Phi = (socc_i->par_0 + j)*optvar + socc_i->par_0;
            {prefix}_mtx_add_direct(Phi+pos_Phi, t_Phi+j*socc_i->colsA, 1, socc_i->colsA);
        }}
    }}
    
    {prefix}_mtx_scale_direct(Phi, kappa, optvar, optvar);
    {prefix}_mtx_add_direct(Phi, H, optvar, optvar);  /* statt 2*H */
    {prefix}_mtx_add_direct(Phi, H, optvar, optvar);
    
    {prefix}_mtx_add_direct(Phi, pbm->Phi_soft, optvar, optvar);
}}

void form_d(real_t *d, const real_t *P, const real_t *h, const real_t *z,
            const uint32_t rowsP, const uint32_t colsP)
{{
    uint32_t i, j;
    for (i = 0; i < rowsP; i++){{
        d[i] = h[i];
        for (j = 0; j < colsP; j++){{
            d[i] -= P[i*colsP + j]*z[j];
        }}
        d[i] = 1/d[i];
    }}
}}

void form_dsoft(real_t *ds, real_t *diags, real_t *rd_soft, real_t *Phi_soft,
                real_t *tmp,
                const real_t *roh, const real_t *z,
                const real_t *Psoft, const real_t *Psoft_T, const real_t *hs,
                const real_t *Fus, const real_t *Fxs, const real_t *Ffs,
                const uint32_t rowsFus, const uint32_t c_veclen,
                const uint32_t rowsFfs, const uint32_t s_veclen,
                const uint32_t T, const struct {prefix}_pbm *pbm)
{{
/*//     real_t rd_soft[T*(c_veclen+s_veclen)];
//     real_t tmp[T*(c_veclen+s_veclen)*(rowsFus*T+rowsFfs)];
//     real_t Phi_soft[T*(c_veclen+s_veclen)*T*(c_veclen+s_veclen)]; */

    
    uint32_t k, i, j;
    for (i = 0; i < rowsFus; i++){{
        ds[i] = hs[i];
        for (j = 0; j < c_veclen; j++){{
            ds[i] -= (Fus+i*c_veclen)[j]*z[j];
        }}
        ds[i] = {prefix}_mtx_smpl_exp(roh[0]*ds[i]);
        diags[i*(rowsFus*T+rowsFfs)+i] = ds[i];
        ds[i] = 1 / (1 + ds[i]);
        diags[i*(rowsFus*T+rowsFfs)+i] = diags[i*(rowsFus*T+rowsFfs)+i] * ds[i] * ds[i];
    }}
    for (k = 1; k < T; k++){{
        for (i = k*rowsFus; i < (k+1)*rowsFus; i++){{
            ds[i] = hs[i];
            for (j = 0; j < s_veclen; j++){{
                ds[i] -= (Fxs+(i-k*rowsFus)*s_veclen)[j]*(z+k*c_veclen+(k-1)*s_veclen)[j];
            }}
            for (j = 0; j < c_veclen; j++){{
                ds[i] -= (Fus+(i-k*rowsFus)*c_veclen)[j]*(z+k*c_veclen+(k)*s_veclen)[j];
            }}
            ds[i] = {prefix}_mtx_smpl_exp(roh[0]*ds[i]);
            diags[i*(rowsFus*T+rowsFfs)+i] = ds[i];
            ds[i] = 1 / (1 + ds[i]);
            diags[i*(rowsFus*T+rowsFfs)+i] = diags[i*(rowsFus*T+rowsFfs)+i] * ds[i] * ds[i];
        }}
    }}
    for (i = T*rowsFus; i < T*rowsFus + rowsFfs; i++){{
        ds[i] = hs[i];
        for (j = 0; j < s_veclen; j++){{
            ds[i] -= (Ffs+(i-T*rowsFus)*s_veclen)[j]*(z+T*c_veclen+(T-1)*s_veclen)[j];
        }}
        ds[i] = {prefix}_mtx_smpl_exp(roh[0]*ds[i]);
        diags[i*(rowsFus*T+rowsFfs)+i] = ds[i];
        ds[i] = 1 / (1 + ds[i]);
        diags[i*(rowsFus*T+rowsFfs)+i] = diags[i*(rowsFus*T+rowsFfs)+i] * ds[i] * ds[i];
    }}
/*//     print_mtx(ds, T*rowsFus+rowsFfs, 1);
//     print_mtx(diags, T*rowsFus+rowsFfs, T*rowsFus+rowsFfs); */
    
    for (i = 0; i < c_veclen; i++){{
        rd_soft[i] = 0;
        for (j = 0; j < rowsFus; j++){{
            rd_soft[i] += Fus[i+j*c_veclen]*ds[j];
        }}
    }}
    for (k = 1; k < T; k++){{
        for (i = k*c_veclen+(k-1)*s_veclen; i < k*c_veclen+k*s_veclen; i++){{
            rd_soft[i] = 0;
            for (j = 0; j < rowsFus; j++){{
                rd_soft[i] += Fxs[i-(k*c_veclen+(k-1)*s_veclen)+j*s_veclen]*(ds+k*rowsFus)[j];
            }}
        }}
        for (i = k*c_veclen+(k)*s_veclen; i < (k+1)*c_veclen+k*s_veclen; i++){{
            rd_soft[i] = 0;
            for (j = 0; j < rowsFus; j++){{
                rd_soft[i] += Fus[i-(k*c_veclen+(k)*s_veclen)+j*c_veclen]*(ds+k*rowsFus)[j];
            }}
        }}
    }}
    for (i = T*c_veclen+(T-1)*s_veclen; i < T*(c_veclen+s_veclen); i++){{
        rd_soft[i] = 0.;
        for (j = 0; j < rowsFfs; j++){{
            rd_soft[i] += Ffs[i-(T*c_veclen+(T-1)*s_veclen)+j*s_veclen] * (ds+T*rowsFus)[j];
        }}
    }}
/*//     print_mtx(rd_soft, T*(c_veclen + s_veclen), 1);
    
//     printf("rd_soft stimmt wenn -1 = %d\n", mtx_cmp(rd_soft, rd_soft_ref, 155, 1e-7));
//     printf("%.17f\t%.17f\n", rd_soft[149], rd_soft_ref[149]); */
    
    
    {prefix}_mtx_multiply_mtx_mtx(tmp, Psoft_T, diags, T*(c_veclen+s_veclen),
                                rowsFus*T+rowsFfs, rowsFus*T+rowsFfs);
    {prefix}_mtx_multiply_mtx_mtx(Phi_soft, tmp, Psoft, T*(c_veclen+s_veclen),
                                rowsFus*T+rowsFfs, T*(c_veclen+s_veclen));
    {prefix}_mtx_scale_direct(Phi_soft, roh[0], T*(c_veclen+s_veclen),
                            T*(c_veclen+s_veclen));
    
/*//     print_mtx(Phi_soft, T*(c_veclen+s_veclen),
//                             T*(c_veclen+s_veclen));
  
    
//     printf("Phi_soft stimmt wenn -1 = %d\n", mtx_cmp(Phi_soft, Phi_soft_ref, 155*155, 1e-6));
//     printf("%.17f\t%.17f\n", Phi_soft[23244], Phi_soft_ref[23244]);
    
    
    
//     {prefix}_mtx_multiply_mtx_mtx(pbm->tmp_Phi_sft_blk->data, diags, Fus,
//                                 rowsFus, rowsFus, c_veclen);
    
//     {prefix}_mtx_multiply_mtx_mtx();
//     print_mtx(pbm->tmp_Phi_sft_blk->data, rowsFus, c_veclen); */
    
}}

void form_diag_d_sq(real_t *diag_d_sq, const real_t *d, const uint32_t dim)
{{
    uint32_t i;
    for (i = 0; i < dim; i++){{
            diag_d_sq[i*dim+i] = d[i]*d[i];
    }}
}}

void {prefix}_pbm_calc_kappa(real_t *kappa, const struct {prefix}_pbm *pbm, const real_t *z)
{{
    real_t *tmp1 = pbm->tmp1_optvar_seqlen;
    real_t *tmp3 = pbm->tmp2_optvar_seqlen;
    real_t *tmp2 = pbm->tmp3_state_veclen; 
    
    {prefix}_mtx_substract(tmp3, pbm->z_opt, pbm->zref, pbm->optvar_seqlen, 1);

    {prefix}_multiply_H_z(tmp1, pbm->H, pbm->z_opt, pbm);

    {prefix}_mtx_multiply_mtx_vec(tmp2, tmp1, pbm->z_opt, 1, pbm->optvar_seqlen);
    {prefix}_mtx_multiply_mtx_vec(kappa, pbm->g, z, 1, pbm->optvar_seqlen);
    kappa[0] = tmp2[0];
    kappa[0] *= 0.01*30/pbm->optvar_seqlen;  /* TODO auf optvar_seqlen umstellen*/
#ifdef HHMPC_QPSMALLTEST
    kappa[0] /= 24;
#endif
#if {PREFIX}_PBM_PRB_COND
    kappa[0] /= 24*3;
/*     kappa[0] = 0.00008;*/
#endif
    kappa[0] += 0;    
/*    printf("calculated kappa = %.20f\n", pbm->kappa[0]);*/

    /* -5 N=5, höchstens -6 N=20, -8 N= 30 */
    kappa[0] = (kappa[0] >5e-8)? kappa[0] : 5e-8;  /* -3 statt -5 für QP cond */
 /* N=30: cond 5*1e-9, uncond 5*1e-6*/
#if {PREFIX}_PBM_PRB_COND
 kappa[0] *=225.15;
#endif
}}

/* static functions definition */

void {prefix}_pbm_warm_start(const struct {prefix}_pbm *pbm)
{{
    real_t *tmp = pbm->tmp3_state_veclen;
    {prefix}_mtx_shift_sequence(pbm->z_ini, pbm->z_opt, pbm->optvar_veclen,
            pbm->optvar_seqlen);
    {prefix}_mtx_multiply_mtx_vec((pbm->z_ini)+pbm->optvar_seqlen-pbm->state_veclen,
                                pbm->A, (pbm->z_opt)+pbm->optvar_seqlen-pbm->state_veclen,
                                pbm->state_veclen, pbm->state_veclen);
    {prefix}_mtx_mul_add((pbm->z_ini)+pbm->optvar_seqlen-pbm->state_veclen, tmp,
                                pbm->B, (pbm->z_opt)+pbm->optvar_seqlen-pbm->optvar_veclen,
                                pbm->state_veclen, pbm->control_veclen);    
    {prefix}_mtx_scale_direct(pbm->z_ini, 1., pbm->optvar_seqlen, 1);
    {prefix}_mtx_shift_sequence(pbm->v_ini, pbm->v_opt, pbm->state_veclen,
            pbm->dual_seqlen);
    
    return;
}}

uint32_t {prefix}_better_step_size(
                const struct {prefix}_pbm *pbm, real_t f_p_g, real_t f_p_alpha)
{{
    if (({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1) ||
        (f_p_g > f_p_alpha) ||
        ({prefix}_pbm_check_positiv(pbm, pbm->z_opt)+1) ){{
            return 0;
        }}
        return -1;
}}

void {prefix}_bt_line_search(real_t *st_size, const struct {prefix}_pbm *pbm)
{{
    const real_t g_step = 1e-6;
    const real_t alpha = 0.25; /* 0.15; */  /* [0.4] Measure for reduction of function value  */
    const real_t beta = 0.6; /*0.5; */  /* [0.6] Factor to decrease step ervery iteration */
    real_t *help_z = pbm->tmp6_optvar_seqlen;
    real_t *help_v = pbm->tmp7_dual_seqlen;
    real_t *t_solve_optvar_seqlen = pbm->tmp1_optvar_seqlen;
    real_t *t_optvar_seqlen = pbm->tmp2_optvar_seqlen;
    real_t st = 1.;
    
    real_t f_p;
    real_t f_p_g;
    real_t g_in_dir;
    
    st_size[0] = st;
    
    /*Save maintain value*/
    memcpy(help_z, pbm->z_opt, pbm->sizeof_optvar_seqlen);
    memcpy(help_v, pbm->v_opt, pbm->sizeof_dual_seqlen);
    
    /* update P matrices */
        {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
        form_d(pbm->d, pbm->P, pbm->h, pbm->z_opt,
               pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
        form_dsoft(pbm->dsoft, pbm->diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
                   pbm->tmp3_mtx_optvar_nb_of_soft,
                   pbm->roh, pbm->z_opt, 
                   pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
                   pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
                   pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
                   pbm->state_veclen, pbm->horizon, pbm);
        residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
    {prefix}_residual_norm(&f_p, pbm->r_d, pbm->r_p,
                  pbm->optvar_seqlen, pbm->dual_seqlen);
    
    {prefix}_mtx_scale(pbm->z_opt, pbm->delta_z, g_step, pbm->optvar_seqlen, 1);
    {prefix}_mtx_add_direct(pbm->z_opt, help_z, pbm->optvar_seqlen, 1);
    {prefix}_mtx_scale(pbm->v_opt, pbm->delta_v, g_step, pbm->dual_seqlen, 1);
    {prefix}_mtx_add_direct(pbm->v_opt, help_v, pbm->dual_seqlen, 1);
    /* update P matrices */
    {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
    form_d(pbm->d, pbm->P, pbm->h, pbm->z_opt, pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
    form_dsoft(pbm->dsoft, pbm->diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
                   pbm->tmp3_mtx_optvar_nb_of_soft,
                   pbm->roh, pbm->z_opt, 
                   pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
                   pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
                   pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
                   pbm->state_veclen, pbm->horizon, pbm);
    residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
    {prefix}_residual_norm(&f_p_g, pbm->r_d, pbm->r_p,
                  pbm->optvar_seqlen, pbm->dual_seqlen);
    g_in_dir = (f_p_g - f_p)/g_step;
/*//     printf("Grad in dir = %.8f\n", g_in_dir);
//     g_in_dir = g_in_dir <= 0 ? g_in_dir : 0.;
//     printf("Grad in dir = %.8f\n", g_in_dir);
//     printf("st size inner = %f\n", st_size[0]); */
    {prefix}_mtx_scale(pbm->z_opt, pbm->delta_z, st_size[0],
                     pbm->optvar_seqlen, 1);
    {prefix}_mtx_add_direct(pbm->z_opt, help_z,
                          pbm->optvar_seqlen, 1);
    {prefix}_mtx_scale(pbm->v_opt, pbm->delta_v, st_size[0],
                     pbm->dual_seqlen, 1);
    {prefix}_mtx_add_direct(pbm->v_opt, help_v,
                          pbm->dual_seqlen, 1);
    /* update P matrices */
    {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
    {prefix}_pbm_check_valid(pbm, pbm->z_opt);
    form_d(pbm->d, pbm->P, pbm->h, pbm->z_opt, pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
    form_dsoft(pbm->dsoft, pbm->diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
                   pbm->tmp3_mtx_optvar_nb_of_soft,
                   pbm->roh, pbm->z_opt, 
                   pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
                   pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
                   pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
                   pbm->state_veclen, pbm->horizon, pbm);
    residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
    {prefix}_residual_norm(&f_p_g, pbm->r_d, pbm->r_p,
                  pbm->optvar_seqlen, pbm->dual_seqlen);
    while (({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1) || (f_p_g > (f_p + alpha*st*g_in_dir)) || ({prefix}_pbm_check_positiv(pbm, pbm->z_opt)+1) )
    {{
        st *= beta;
        if (st < g_step){{st = 0.; break;}}
        
        {prefix}_mtx_scale(pbm->z_opt, pbm->delta_z, st,
                         pbm->optvar_seqlen, 1);
        {prefix}_mtx_add_direct(pbm->z_opt, help_z,
                              pbm->optvar_seqlen, 1);
        {prefix}_mtx_scale(pbm->v_opt, pbm->delta_v, st,
                         pbm->dual_seqlen, 1);
        {prefix}_mtx_add_direct(pbm->v_opt, help_v,
                              pbm->dual_seqlen, 1);
        /* update P matrices */
        {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
        form_d(pbm->d, pbm->P, pbm->h, pbm->z_opt,
               pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
        form_dsoft(pbm->dsoft, pbm->diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
                   pbm->tmp3_mtx_optvar_nb_of_soft,
                   pbm->roh, pbm->z_opt, 
                   pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
                   pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
                   pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
                   pbm->state_veclen, pbm->horizon, pbm);
        residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
        {prefix}_residual_norm(&f_p_g, pbm->r_d, pbm->r_p,
                      pbm->optvar_seqlen, pbm->dual_seqlen);
    }}
    
    /*Load back maintain value*/
    memcpy(pbm->z_opt, help_z, pbm->sizeof_optvar_seqlen);
    memcpy(pbm->v_opt, help_v, pbm->sizeof_dual_seqlen);
    st_size[0] = st;
}}

void {prefix}_bt_line_search_new(
                real_t *st_size, const struct {prefix}_pbm *pbm)
{{
    const real_t g_step = 1e-6;
    const real_t alpha = 0.25; /* 0.15; */  /* [0.4] Measure for reduction of function value  */
    const real_t beta = 0.6; /* 0.5; */  /* [0.6] Factor to decrease step ervery iteration */
    real_t *help_z = pbm->tmp6_optvar_seqlen;
    real_t *help_v = pbm->tmp7_dual_seqlen;
    real_t *t_solve_optvar_seqlen = pbm->tmp1_optvar_seqlen;
    real_t *t_optvar_seqlen = pbm->tmp2_optvar_seqlen;
    real_t st = 1.;
    
    real_t f_p;
    real_t f_p_g;
    real_t g_in_dir;
    
    uint32_t i;
    
    st_size[0] = st;
    
    /*Save maintain value*/
    memcpy(help_z, pbm->z_opt, pbm->sizeof_optvar_seqlen);
    memcpy(help_v, pbm->v_opt, pbm->sizeof_dual_seqlen);
    
    
    {prefix}_residual_norm(&f_p, pbm->r_d, pbm->r_p,
                  pbm->optvar_seqlen, pbm->dual_seqlen);
    
    {prefix}_mtx_scale(pbm->z_opt, pbm->delta_z, g_step, pbm->optvar_seqlen, 1);
    {prefix}_mtx_add_direct(pbm->z_opt, help_z, pbm->optvar_seqlen, 1);
    {prefix}_mtx_scale(pbm->v_opt, pbm->delta_v, g_step, pbm->dual_seqlen, 1);
    {prefix}_mtx_add_direct(pbm->v_opt, help_v, pbm->dual_seqlen, 1);
    /* update P matrices */
    {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
    form_d(pbm->d, pbm->P, pbm->h, pbm->z_opt, pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
    form_dsoft(pbm->dsoft, pbm->diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
                   pbm->tmp3_mtx_optvar_nb_of_soft,
                   pbm->roh, pbm->z_opt, 
                   pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
                   pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
                   pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
                   pbm->state_veclen, pbm->horizon, pbm);
    residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
    {prefix}_residual_norm(&f_p_g, pbm->r_d, pbm->r_p,
                  pbm->optvar_seqlen, pbm->dual_seqlen);
    g_in_dir = (f_p_g - f_p)/g_step;
/* //     printf("Grad in dir = %.8f\n", g_in_dir);
//     g_in_dir = g_in_dir <= 0 ? g_in_dir : 0.;
//     printf("Grad in dir = %.8f\n", g_in_dir);
//     printf("st size inner = %f\n", st_size[0]);
//     {prefix}_mtx_scale(pbm->z_opt, pbm->delta_z, st_size[0],
//                      pbm->optvar_seqlen, 1);
//     {prefix}_mtx_add_direct(pbm->z_opt, help_z,
//                           pbm->optvar_seqlen, 1);
//     {prefix}_mtx_scale(pbm->v_opt, pbm->delta_v, st_size[0],
//                      pbm->dual_seqlen, 1);
//     {prefix}_mtx_add_direct(pbm->v_opt, help_v,
//                           pbm->dual_seqlen, 1);
//      *//* update P matrices *//*
//     {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
//                t_solve_optvar_seqlen, t_optvar_seqlen);
//     {prefix}_pbm_check_valid(pbm, pbm->z_opt);
//     form_d(pbm->d, pbm->P, pbm->h, pbm->z_opt, pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
//     form_dsoft(pbm->dsoft, pbm->diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
//                    pbm->tmp3_mtx_optvar_nb_of_soft,
//                    pbm->roh, pbm->z_opt, 
//                    pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
//                    pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
//                    pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
//                    pbm->state_veclen, pbm->horizon, pbm);
//     residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
//     {prefix}_residual_norm(&f_p_g, pbm->r_d, pbm->r_p,
//                   pbm->optvar_seqlen, pbm->dual_seqlen);
//     while (({prefix}_pbm_check_valid(pbm, pbm->z_opt)+1) || (f_p_g > (f_p + alpha*st*g_in_dir)) || ({prefix}_pbm_check_positiv(pbm, pbm->z_opt)+1) )
//     while ({prefix}_better_step_size(pbm, f_p_g, (f_p + alpha*st*g_in_dir))) */
    st = 1.0234904e-6;
    st = 0.0060466176;
    st = 0.0279936;
    st_size[0] = g_step;
    for (i = 0; i < 8; i++)
    {{
        {prefix}_mtx_scale(pbm->z_opt, pbm->delta_z, st,
                         pbm->optvar_seqlen, 1);
        {prefix}_mtx_add_direct(pbm->z_opt, help_z,
                              pbm->optvar_seqlen, 1);
        {prefix}_mtx_scale(pbm->v_opt, pbm->delta_v, st,
                         pbm->dual_seqlen, 1);
        {prefix}_mtx_add_direct(pbm->v_opt, help_v,
                              pbm->dual_seqlen, 1);
        /* update P matrices */
        {prefix}_pbm_update_P(pbm->P_of_z, pbm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
        form_d(pbm->d, pbm->P, pbm->h, pbm->z_opt,
               pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
        form_dsoft(pbm->dsoft, pbm->diag_d_soft, pbm->r_d_soft, pbm->Phi_soft,
                   pbm->tmp3_mtx_optvar_nb_of_soft,
                   pbm->roh, pbm->z_opt, 
                   pbm->Psoft, pbm->Psoft_T, pbm->hsoft,
                   pbm->Fusoft, pbm->Fxsoft, pbm->Ffsoft,
                   pbm->rowsFusoft, pbm->control_veclen, pbm->rowsFfsoft,
                   pbm->state_veclen, pbm->horizon, pbm);
        residual(pbm, pbm->z_opt, pbm->v_opt, pbm->d, pbm->kappa[0]);
        {prefix}_residual_norm(&f_p_g, pbm->r_d, pbm->r_p,
                      pbm->optvar_seqlen, pbm->dual_seqlen);
/*         printf("step %f\n", st); */
        st_size[0] = {prefix}_better_step_size(pbm, f_p_g, (f_p + alpha*st*g_in_dir)) ? st : st_size[0];
        st /= beta;
    }}
    
    /*Load back maintain value*/
    memcpy(pbm->z_opt, help_z, pbm->sizeof_optvar_seqlen);
    memcpy(pbm->v_opt, help_v, pbm->sizeof_dual_seqlen);
/*     st_size[0] = st; */
}}

void {prefix}_residual_norm(
                real_t *f, const real_t *r_d, const real_t *r_p,
                const uint32_t ld, const uint32_t lp)
{{
    uint32_t i;
    f[0] = 0.;
    /* Better solution if only norm(r_p) is calculated */
/*//     for (i = 0; i < ld; i++){{
//         f[0] += r_d[i]*r_d[i];
//     }} */
    for (i = 0; i < lp; i++){{
        f[0] += r_p[i]*r_p[i];
    }}
}}

void {prefix}_multiply_P_T_d(real_t *product, const real_t *P_T, const real_t *d,
                        const struct {prefix}_pbm* pbm)
{{
#if {PREFIX}_PBM_PRB_SPARSE
    uint32_t i, j, k;
    for (i = 0; i < pbm->control_veclen; i++){{
        product[i] = 0.;
        for (j = 0; j < HHMPC_NB_ROWSFU; j++){{
            product[i] += P_T[i*HHMPC_NB_LCONSTR+j]*d[j];
        }}
    }}
    for (k = 0; k < pbm->horizon-1; k++){{
        for (i = pbm->control_veclen + k*pbm->optvar_veclen; i < pbm->control_veclen + (k+1)*pbm->optvar_veclen; i++){{
            product[i] = 0.;
            for (j = HHMPC_NB_ROWSFU + k*HHMPC_NB_ROWSFU; j < HHMPC_NB_ROWSFU + (k+1)*HHMPC_NB_ROWSFU; j++){{
                product[i] += P_T[i*HHMPC_NB_LCONSTR+j]*d[j];
            }}
        }}
    }}
    for (i = pbm->control_veclen + k*pbm->optvar_veclen; i < (k+1)*pbm->optvar_veclen; i++){{
        product[i] = 0.;
        for (j = HHMPC_NB_ROWSFU + k*HHMPC_NB_ROWSFU; j < HHMPC_NB_ROWSFU + k*HHMPC_NB_ROWSFU + HHMPC_NB_ROWSFF; j++){{
            product[i] += P_T[i*HHMPC_NB_LCONSTR+j]*d[j];
        }}
    }}
#endif
#if {PREFIX}_PBM_PRB_COND
    {prefix}_mtx_multiply_mtx_vec(product, P_T, d,
                                  pbm->optvar_seqlen, pbm->nb_of_ueq_constr);
#endif
}}

void {prefix}_multiply_P_z(real_t *product, const real_t *P, const real_t *z,
                        const struct {prefix}_pbm* pbm)
{{
#if {PREFIX}_PBM_PRB_SPARSE
    uint32_t i, j, k;
    for (i = 0; i < HHMPC_NB_ROWSFU; i++){{
        product[i] = 0.;
        for (j = 0; j < pbm->control_veclen; j++){{
            product[i] += P[i*pbm->optvar_seqlen+j]*z[j];
        }}
    }}
    for (k = 0; k < pbm->horizon-1; k++){{
        for (i = HHMPC_NB_ROWSFU + k*HHMPC_NB_ROWSFU; i < HHMPC_NB_ROWSFU + (k+1)*HHMPC_NB_ROWSFU; i++){{
            product[i] = 0.;
            for (j = pbm->control_veclen + k*pbm->optvar_veclen; j < pbm->control_veclen + (k+1)*pbm->optvar_veclen; j++){{
                product[i] += P[i*pbm->optvar_seqlen+j]*z[j];
            }}
        }}
    }}
    for (i = HHMPC_NB_ROWSFU + k*HHMPC_NB_ROWSFU; i < HHMPC_NB_ROWSFU + k*HHMPC_NB_ROWSFU + HHMPC_NB_ROWSFF; i++){{
        product[i] = 0.;
        for (j = pbm->control_veclen + k*pbm->optvar_veclen; j < (k+1)*pbm->optvar_veclen; j++){{
            product[i] += P[i*pbm->optvar_seqlen+j]*z[j];
        }}
    }}
#endif
#if {PREFIX}_PBM_PRB_COND
  {prefix}_mtx_multiply_mtx_vec(product, P, z,
                                pbm->nb_of_ueq_constr, pbm->optvar_seqlen);
#endif
}}

void {prefix}_multiply_H_z(real_t *product, const real_t *H, const real_t *z,
                        const struct {prefix}_pbm* pbm)
{{
#if {PREFIX}_PBM_PRB_SPARSE
    uint32_t i, j, k;
    for (i = 0; i < pbm->control_veclen; i++){{
        product[i] = 0.;
        for (j = 0; j < pbm->control_veclen; j++){{
            product[i] += H[i*pbm->optvar_seqlen+j]*z[j];
        }}
    }}
    for (k = 0; k < pbm->horizon-1; k++){{
        for (i = pbm->control_veclen + k*pbm->optvar_veclen; i < pbm->control_veclen + (k+1)*pbm->optvar_veclen; i++){{
            product[i] = 0.;
            for (j = pbm->control_veclen + k*pbm->optvar_veclen; j < pbm->control_veclen + (k+1)*pbm->optvar_veclen; j++){{
                product[i] += H[i*pbm->optvar_seqlen+j]*z[j];
            }}
        }}
    }}
    for (i = pbm->control_veclen + k*pbm->optvar_veclen; i < (k+1)*pbm->optvar_veclen; i++){{
        product[i] = 0.;
        for (j = pbm->control_veclen + k*pbm->optvar_veclen; j < (k+1)*pbm->optvar_veclen; j++){{
            product[i] += H[i*pbm->optvar_seqlen+j]*z[j];
        }}
    }}
#endif    
#if {PREFIX}_PBM_PRB_COND
    {prefix}_mtx_multiply_mtx_vec(product, H, z,
                                  pbm->optvar_seqlen, pbm->optvar_seqlen);
#endif
}}
