/* Solve systems of linear equation arising in used primal barrier interior
 * point method.
 */ 

#include "{prefix}pbmsolve.h"/*
#include "static_data.h"   *//* TODO implement static_data.h */

/* external functions definition */

void {prefix}_multiply_C_z(
                real_t *product, const real_t *C, const real_t *z,
                const uint32_t n, const uint32_t m, const uint32_t T)
{{
    uint32_t optvar_seqlen = T*(n + m);
#if {PREFIX}_PBM_PRB_SPARSE
    uint32_t i, j, k;
    uint32_t optvar_veclen = (n + m);
    for (i = 0; i < n; i++){{
        product[i] = 0.;
        for (j = 0; j < optvar_seqlen; j++){{
            product[i] += C[i*optvar_seqlen+j]*z[j];
        }}
    }}
    for (k = 1; k < T; k++){{
        for (i = k*n; i < (k+1)*n; i++){{
            product[i] = 0.;
            for (j = m + (k-1)*optvar_veclen; j < (k+1)*optvar_veclen; j++){{
                product[i] += C[i*optvar_seqlen+j]*z[j];
            }}
        }}
    }}
#endif
#if {PREFIX}_PBM_PRB_COND
    uint32_t dual_seqlen = T*n;
    {prefix}_mtx_multiply_mtx_vec(product, C, z, dual_seqlen, optvar_seqlen);
#endif
}}

void {prefix}_multiply_C_T_v(
                real_t *product, const real_t *C_T, const real_t *v,
                const uint32_t n, const uint32_t m, const uint32_t T)
{{
    
    uint32_t dual_seqlen = T*n;
#if {PREFIX}_PBM_PRB_SPARSE
    uint32_t i, j, k;
    uint32_t optvar_veclen = (n + m);
    for (k = 0; k < T-1; k++){{
        for (i = k*optvar_veclen; i < (k+1)*optvar_veclen; i++){{
            product[i] = 0.;
            for (j = k*n; j < (k+2)*n; j++){{
                product[i] += C_T[i*dual_seqlen+j]*v[j];
            }}
        }}
    }}
    for (i = k*optvar_veclen; i < (k+1)*optvar_veclen; i++){{
        product[i] = 0.;
        for (j = k*n; j < (k+1)*n; j++){{
            product[i] += C_T[i*dual_seqlen+j]*v[j];
        }}
    }}
#endif
#if {PREFIX}_PBM_PRB_COND
    uint32_t optvar_seqlen = T*(n + m);
    {prefix}_mtx_multiply_mtx_vec(product, C_T, v, optvar_seqlen, dual_seqlen);
#endif
}}

void solve_sysofleq(real_t delta_z[], real_t delta_v[],
                    const real_t Phi[],
                    const real_t rd[], const real_t rp[],
                    const real_t C[], const real_t *C_T,
                    const real_t A[], const real_t A_T[],
                    const real_t B[], const real_t B_T[],
                    const real_t A_B[], const real_t A_B_T[],
                    const real_t H[],
                    const uint32_t n, const uint32_t m, const uint32_t T,
                    const real_t *eye_nm, const real_t *eye_n,
                    real_t *tmp1_optvar_seqlen,
                    real_t *tmp2_optvar_seqlen,
                    real_t *tmp_dual_seqlen,
                    real_t *L_Y, real_t *L_Y_T,
                    real_t *PhiBlock, real_t *PhiBlock_I,
                    real_t *PhiBlock_I_last,
                    real_t *Block_nxn1, real_t *Block_nxn2,
                    real_t *tmp_optvar_veclenxoptvar_veclen,
                    real_t *L_Phi_blocks, real_t *L_Phi_T_blocks)
{{
    real_t *tmp3_state_veclen = tmp2_optvar_seqlen; /* In a other function */
    uint32_t optvar_seqlen;
    optvar_seqlen = (n+m)*T;
    /* TODO Change algorithm, that it works without extra part for horizon=1 */
    /* If only one block for control variables without states. */
    if (optvar_seqlen == m){{
        {prefix}_mtx_scale(tmp2_optvar_seqlen, rd, -1, optvar_seqlen, 1);
        /*
        cholesky(L_Phi_blocks, Phi, HHMPC_OS);
        */
        /* Check if cholesky factorization failed. */
        /*
        if (L_Phi_blocks[optvar_seqlen*optvar_seqlen-1] !=
                L_Phi_blocks[optvar_seqlen*optvar_seqlen-1]){{
            cholesky(L_Phi_blocks, H, HHMPC_OS);*/  /* TODO Find a better way to ensure taht the algorithm does not fail. */
            /* Idea : Solve System of linear eqations only for the unconstrained
            * problem. */
            /*
            {prefix}_mtx_transpose(L_Phi_T_blocks, L_Phi_blocks, HHMPC_OS, HHMPC_OS);
            fwd_subst(tmp1_optvar_seqlen, L_Phi_blocks, HHMPC_OS, tmp2_optvar_seqlen, 1);
        */
            /* Return delta_z = 0. This should be replaced by a better solution,
            * as decribed above. */
            /*
            {prefix}_mtx_to_zero(delta_z, optvar_seqlen);
            return;
        }}
        {prefix}_mtx_transpose(L_Phi_T_blocks, L_Phi_blocks, HHMPC_OS, HHMPC_OS);
        fwd_subst(tmp1_optvar_seqlen, L_Phi_blocks, HHMPC_OS, tmp2_optvar_seqlen, 1);
        */
    /*//     print_mtx(rd, 5, 1);
    //     print_mtx(delta_z, 5, 1); */
    /*
        bwd_subst(delta_z, L_Phi_T_blocks, HHMPC_OS, tmp1_optvar_seqlen, 1);
        */
    /*     print_mtx(delta_z, 5 ,1); */
        solveBlock(delta_z, L_Phi_blocks, L_Phi_T_blocks, Phi,
                optvar_seqlen, tmp2_optvar_seqlen, 1, tmp1_optvar_seqlen);

    }}else{{
        form_Y(L_Y, L_Y_T, L_Phi_blocks, L_Phi_T_blocks,
            Phi, T, A_B, A_B_T, n, B, B_T, m, eye_nm, eye_n,
            PhiBlock, PhiBlock_I, PhiBlock_I_last,
            Block_nxn1, Block_nxn2, tmp_optvar_veclenxoptvar_veclen);

        /* temporären Zeiger für beta sparen */
        form_beta(delta_v, L_Phi_blocks, L_Phi_T_blocks, rd, rp, T, C, n, m,
                tmp1_optvar_seqlen, tmp2_optvar_seqlen);
        form_delta_v(delta_v, tmp_dual_seqlen, tmp3_state_veclen, 
                    L_Y, L_Y_T, T, n);
        form_delta_z(delta_z, tmp1_optvar_seqlen, delta_v,
                    L_Phi_blocks, L_Phi_T_blocks, rd, C_T, T, n, m);
    }}
}}

void form_delta_z(real_t delta_z[],
                  real_t *tmp_optvar_seqlen,
                  const real_t delta_v[],
                  const real_t L_Phi[],
                  const real_t L_Phi_T[],
                  const real_t rd[],
                  const real_t C_T[],
                  const uint32_t T, const uint32_t n, const uint32_t m)
{{
    uint32_t i;
    
    
    {prefix}_multiply_C_T_v(delta_z, C_T, delta_v, n, m, T);
/*     {prefix}_mtx_multiply_mtx_vec(delta_z, C_T, delta_v, T*(n+m), T*n); */
    {prefix}_mtx_add_direct(delta_z, rd, T*(n+m), 1);
    {prefix}_mtx_scale_direct(delta_z, -1, T*(n+m), 1);
/*     print_mtx(delta_z, T*(n+m), 1); */
    
/*//     fwd_subst(tmp_optvar_seqlen, L_Phi, T*(n+m), delta_z, 1);
//     print_mtx(tmp_optvar_seqlen, T*(n+m), 1);
//     bwd_subst(delta_z, L_Phi_T, T*(n+m), tmp_optvar_seqlen, 1); */
    
    {prefix}_mtx_fwd_subst(tmp_optvar_seqlen, L_Phi, m, delta_z, 1);
    for (i = 0; i < T-1; i++){{
        {prefix}_mtx_fwd_subst(tmp_optvar_seqlen+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n+m, delta_z+m+i*(n+m), 1);
    }}
    {prefix}_mtx_fwd_subst(tmp_optvar_seqlen+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n, delta_z+m+i*(n+m), 1);
/*     print_mtx(tmp_optvar_seqlen, T*(n+m), 1); */
    
    {prefix}_mtx_bwd_subst(delta_z+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n, tmp_optvar_seqlen+m+i*(n+m), 1);
    for (i = T-1; i > 0; i--){{
        {prefix}_mtx_bwd_subst(delta_z+m+(i-1)*(n+m), L_Phi_T+m*m+(i-1)*(n+m)*(n+m), n+m, tmp_optvar_seqlen+m+(i-1)*(n+m), 1);
    }}
        /*bwd_subst(delta_z+m+(i-1)*(n+m), L_Phi_T+m*m+(i-1)*(n+m)*(n+m), n+m, tmp_optvar_seqlen+m+(i-1)*(n+m), 1);*/
    {prefix}_mtx_bwd_subst(delta_z, L_Phi_T, m, tmp_optvar_seqlen, 1);
/*     print_mtx(delta_z, T*(n+m), 1); */
}}

void form_beta(real_t beta[],
               const real_t L_Phi[],
               const real_t L_Phi_T[],
               const real_t rd[], const real_t rp[],
               const uint32_t T,
               const real_t C[], const uint32_t n, const uint32_t m,
               real_t *tmp1, real_t *tmp2
               /*const real_t A[], const uint32_t n,
               const real_t B[], const uint32_t m*/)
{{
    uint32_t i;
    /* TODO beta lässt sich sicher auch parallel zu Y formen */
    
/*//     fwd_subst(tmp1, L_Phi, T*(n+m), rd, 1);
//     bwd_subst(tmp2, L_Phi_T, T*(n+m), tmp1, 1); */
    
    {prefix}_mtx_fwd_subst(tmp1, L_Phi, m, rd, 1);
    for (i = 0; i < T-1; i++){{
        {prefix}_mtx_fwd_subst(tmp1+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n+m, rd+m+i*(n+m), 1);
    }}
    {prefix}_mtx_fwd_subst(tmp1+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n, rd+m+i*(n+m), 1);
    
    {prefix}_mtx_bwd_subst(tmp2+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n, tmp1+m+i*(n+m), 1);
    for (i = T-1; i > 0; i--){{
        {prefix}_mtx_bwd_subst(tmp2+m+(i-1)*(n+m), L_Phi_T+m*m+(i-1)*(n+m)*(n+m), n+m, tmp1+m+(i-1)*(n+m), 1);
    }}
    /*bwd_subst(tmp2+m+(i-1)*(n+m), L_Phi_T+m*m+(i-1)*(n+m)*(n+m), n+m, tmp1+m+(i-1)*(n+m), 1);*/
    {prefix}_mtx_bwd_subst(tmp2, L_Phi_T, m, tmp1, 1);
    
    {prefix}_multiply_C_z(tmp1, C, tmp2, n, m, T);
/*     {prefix}_mtx_multiply_mtx_vec(tmp1, C, tmp2, T*n, T*(n+m)); */
    {prefix}_mtx_substract(beta, tmp1, rp, T*n, 1);
/*     print_mtx(beta, T*n, 1); */
    {prefix}_mtx_scale_direct(beta, -1., T*n, 1);  /* return -beta */
}}

void form_delta_v(real_t delta_v[],
                  real_t *tmp_dual_seqlen, real_t *tmp_n,
                  const real_t L_Y[], const real_t L_Y_T[],
                  const uint32_t T, const uint32_t n)
{{
    uint32_t i;
    /* {prefix}_mtx_scale(delta_v, beta, -1., T*n, 1); */
    
    for (i = 0; i < T-1; i++) {{
        {prefix}_mtx_fwd_subst(tmp_dual_seqlen+i*n, L_Y+2*i*n*n, n, delta_v+i*n, 1);
        {prefix}_mtx_multiply_mtx_vec(tmp_n, L_Y+2*i*n*n+n*n, tmp_dual_seqlen+i*n, n, n);
        {prefix}_mtx_substract_direct(delta_v+i*n+n, tmp_n, n, 1);
    }}
    {prefix}_mtx_fwd_subst(tmp_dual_seqlen+i*n, L_Y+2*i*n*n, n, delta_v+i*n, 1); /*i=T-1*/
    
    for (i = T-1; i > 0; i--) {{
        {prefix}_mtx_bwd_subst(delta_v+i*n, L_Y_T+2*i*n*n, n, tmp_dual_seqlen+i*n, 1);
        {prefix}_mtx_multiply_mtx_vec(tmp_n, L_Y_T+2*i*n*n-n*n, delta_v+i*n, n, n);
        {prefix}_mtx_substract_direct(tmp_dual_seqlen+i*n-n, tmp_n, n, 1);
    }}
    {prefix}_mtx_bwd_subst(delta_v+i*n, L_Y_T+2*i*n*n, n, tmp_dual_seqlen+i*n, 1); /*i=0*/
}}

void form_Y(real_t L_Y[], real_t *L_Y_T, real_t L_Phi[], real_t *L_Phi_T,
            real_t Phi[],
            const uint32_t T,
            const real_t *A_B, const real_t *A_T_B_T, const uint32_t n,
            const real_t B[], const real_t *B_T, const uint32_t m,
            const real_t *eye_nm, const real_t *eye_n,
            real_t *PhiBlock, real_t *PhiBlock_I, real_t *last_PhiBlock_I,
            real_t *Qi_tilde, real_t *Y_bl, real_t *hilf1)
{{
    uint32_t i, j, ri, bl;
    /* hilf1 auch mehrmal als temporäre Variable für [n*n] und andere Größen verwendet */
    
    for (i = 0; i < T; i++){{
        for (j = 0; j < (n+m)*(n+m); j++){{  /* all blocks */
            last_PhiBlock_I[j] = PhiBlock_I[j];
        }}
        if ( i == T-1){{  /* last block i = T-1 */
            getBlock(PhiBlock, Phi, T*(n+m), m+i*(n+m), m+i*(n+m), n, n);
            bl = m*m + i*(n+m)*(n+m);
/*//            cholesky(L_Phi+bl, PhiBlock, n);
//            {prefix}_mtx_transpose(L_Phi_T+bl, L_Phi+bl, n, n);
//            fwd_subst(hilf1, L_Phi+bl, n, eye_n, n);
//            bwd_subst(Qi_tilde, L_Phi_T+bl, n, hilf1, n); */
            solveBlock(Qi_tilde, L_Phi+bl, L_Phi_T+bl,
                       PhiBlock, n, eye_n, n, hilf1);
            form_Yii(Y_bl, A_B, n, n+m, last_PhiBlock_I, n+m, n+m, A_T_B_T,
                     n+m, n, Qi_tilde, hilf1);
             /* regularization (-epsilon*I) */
/*//             for (ri = 0; ri < n; ri++){{
//                 (Y_bl+ri*n+ri)[0] += reg[0];
//             }}
//             setBlock(Y, T*n, Y_bl, n, n, i*n, i*n); */
            
            {prefix}_mtx_multiply_mtx_mtx(hilf1, L_Y+2*i*(n*n)-(n*n), L_Y_T+2*i*(n*n)-(n*n), n, n, n);
            {prefix}_mtx_scale_direct(hilf1, -1, n, n);
            {prefix}_mtx_add_direct(hilf1, Y_bl, n, n);
            {prefix}_mtx_cholesky(L_Y+2*i*(n*n), hilf1, n);
            {prefix}_mtx_transpose(L_Y_T+2*i*(n*n), L_Y+2*i*(n*n), n, n);
        }}
        
        if (i < T-1){{
            getBlock(PhiBlock, Phi, T*(n+m), m+i*(n+m), m+i*(n+m), n+m, n+m);
            bl = m*m + i*(n+m)*(n+m);
/*//            cholesky(L_Phi+bl, PhiBlock, n+m);
//            {prefix}_mtx_transpose(L_Phi_T+bl, L_Phi+bl, n+m, n+m);
//            fwd_subst(hilf1, L_Phi+bl, n+m, eye_nm, n+m);
//            bwd_subst(PhiBlock_I, L_Phi_T+bl, n+m, hilf1, n+m); */
            solveBlock(PhiBlock_I, L_Phi+bl, L_Phi_T+bl,
                       PhiBlock, n+m, eye_nm, n+m, hilf1);
            getBlock(Qi_tilde, PhiBlock_I, n+m, 0, 0, n, n);
            
            if (i == 0){{  /* first Block: i = 0 */
                getBlock(PhiBlock, Phi, T*(n+m), 0, 0, m, m);
/*//                cholesky(L_Phi, PhiBlock, m);
//                {prefix}_mtx_transpose(L_Phi_T, L_Phi, m, m);
//                fwd_subst(hilf1, L_Phi, m, B_T, n); 
//                bwd_subst(hilf1+(m*n), L_Phi_T, m, hilf1, n); */
                solveBlock(hilf1+(m*n), L_Phi, L_Phi_T,
                       PhiBlock, m, B_T, n, hilf1);
                form_Y11(Y_bl, B, B_T, n, m, L_Phi, L_Phi_T, Qi_tilde,
                    hilf1, hilf1+(m*n));
                 /* regularization (-epsilon*I) */
/*//                 for (ri = 0; ri < n; ri++){{
//                     (Y_bl+ri*n+ri)[0] += reg[0];
//                 }}
//                 setBlock(Y, T*n, Y_bl, n, n, i*n, i*n); */
                /* hilf1 has size (n+m)*(n+m), so there is enough space for all */
                
                {prefix}_mtx_cholesky(L_Y+2*i*(n*n), Y_bl, n);
                {prefix}_mtx_transpose(L_Y_T+2*i*(n*n), L_Y+2*i*(n*n), n, n);
            }}
            
            if (i > 0) {{  /* not first block i != 0 */
                form_Yii(Y_bl, A_B, n, n+m, last_PhiBlock_I, n+m, n+m, A_T_B_T,
                         n+m, n, Qi_tilde, hilf1);
                 /* regularization (-epsilon*I) */
/*//                 for (ri = 0; ri < n; ri++){{
//                     (Y_bl+ri*n+ri)[0] += reg[0];
//                 }}
//                 setBlock(Y, T*n, Y_bl, n, n, i*n, i*n); */
                
                {prefix}_mtx_multiply_mtx_mtx(hilf1, L_Y+2*i*(n*n)-(n*n), L_Y_T+2*i*(n*n)-(n*n), n, n, n);
                {prefix}_mtx_scale_direct(hilf1, -1, n, n);
                {prefix}_mtx_add_direct(hilf1, Y_bl, n, n);
                {prefix}_mtx_cholesky(L_Y+2*i*(n*n), hilf1, n);
                {prefix}_mtx_transpose(L_Y_T+2*i*(n*n), L_Y+2*i*(n*n), n, n);
            }}
            
            form_Y_i_ip1(Y_bl, A_T_B_T, n+m, n, PhiBlock_I);
/*//             setBlock(Y, T*n, Y_bl, n, n, i*n, (i+1)*n);
//             {prefix}_mtx_transpose(hilf1, Y_bl, n, n);
//             setBlock(Y, T*n, hilf1, n, n, (i+1)*n, i*n); */
            
            {prefix}_mtx_fwd_subst(L_Y_T+2*i*(n*n)+(n*n), L_Y+2*i*(n*n), n, Y_bl, n);
            {prefix}_mtx_transpose(L_Y+2*i*(n*n)+(n*n), L_Y_T+2*i*(n*n)+(n*n), n, n);
        }}
    }}
}}

void solveBlock(real_t *mtxA_I, real_t *L, real_t *L_T,
                real_t *mtxA, const uint32_t dim, const real_t *mtxB, const uint32_t colsB,
                real_t *tmp)
{{
/*//     uint32_t ri;
//     real_t reg[] = {{.001}};
//     real_t dB[dim*colsB];
//     real_t dlsg[dim*colsB];
//     for (ri = 0; ri < dim; ri++){{
//         (mtxA+ri*dim+ri)[0] += reg[0];
//     }} */
    {prefix}_mtx_cholesky(L, mtxA, dim);
    {prefix}_mtx_transpose(L_T, L, dim, dim);
    {prefix}_mtx_fwd_subst(tmp, L, dim, mtxB, colsB);
    {prefix}_mtx_bwd_subst(mtxA_I, L_T, dim, tmp, colsB);
    
/*//     for (ri = 0; ri < dim; ri++){{
//         (mtxA+ri*dim+ri)[0] -= reg[0];
//     }}
// */     /* iterative refinement */ /*
//     {prefix}_mtx_multiply_mtx_mtx(dB, mtxA, mtxA_I, dim, dim, colsB);
//     {prefix}_mtx_scale_direct(dB, -1., dim, colsB);
//     
//     {prefix}_mtx_add_direct(dB, mtxB, dim, colsB);
// //     print_mtx(dB, dim, colsB );
//     fwd_subst(tmp, L, dim, dB, colsB);
//     bwd_subst(dlsg, L_T, dim, tmp, colsB);
// //     print_mtx(dlsg, dim, colsB );
//     {prefix}_mtx_add_direct(mtxA_I, dlsg, dim, colsB);
//     
//     {prefix}_mtx_multiply_mtx_mtx(dB, mtxA, mtxA_I, dim, dim, colsB);
//     {prefix}_mtx_scale_direct(dB, -1., dim, colsB);
//     {prefix}_mtx_add_direct(dB, mtxB, dim, colsB);
//     
//     fwd_subst(tmp, L, dim, dB, colsB);
//     bwd_subst(dlsg, L_T, dim, tmp, colsB);
// //     print_mtx(dlsg, dim, colsB );
//     {prefix}_mtx_add_direct(mtxA_I, dlsg, dim, colsB);
//     
//     {prefix}_mtx_multiply_mtx_mtx(dB, mtxA, mtxA_I, dim, dim, colsB);
//     {prefix}_mtx_scale_direct(dB, -1., dim, colsB);
//     {prefix}_mtx_add_direct(dB, mtxB, dim, colsB);
//     print_mtx(dB, dim, colsB ); */
}} 
/* TODO Maybe for the following three functions, one function is sufficient. */
void form_Yii(real_t sol[],
              const real_t A_B[], const uint32_t rowsA, const uint32_t colsA,
              const real_t last_PhiBlock_I[], const uint32_t rowsB, const uint32_t colsB,
              const real_t A_T_B_T[], const uint32_t rowsC, const uint32_t colsC,
              const real_t Qi[],
              real_t *tmp)
{{
    {prefix}_mtx_multiply_mtx_mtx(tmp, last_PhiBlock_I, A_T_B_T, rowsB, colsB, colsC);
    {prefix}_mtx_multiply_mtx_mtx(sol, A_B, tmp, rowsA, colsA, colsC);
    {prefix}_mtx_add_direct(sol, Qi, rowsA, rowsA);
}}

void form_Y_i_ip1(real_t sol[],
                  const real_t A_T_B_T[], const uint32_t rows,
                  const uint32_t cols,
                  const real_t QiSi_Block_I[])
{{    
    {prefix}_mtx_multiply_mtx_mtx(sol, QiSi_Block_I, A_T_B_T, cols, rows, cols);
    {prefix}_mtx_scale_direct(sol, -1, cols, cols);
}}

void form_Y11(real_t sol[],
              const real_t B[], const real_t *B_T,
              const uint32_t n, const uint32_t m,
              const real_t L_R0[], const real_t *L_R0_T,
              const real_t Q_tilde[],
              real_t *tmp1_mxn, real_t *tmp2_mxn)
{{
    {prefix}_mtx_multiply_mtx_mtx(sol, B, tmp2_mxn, n, m, n);
    {prefix}_mtx_add_direct(sol, Q_tilde, n, n);
}}

void setBlock(real_t mtx[], const uint32_t dim,
              const real_t blk[], const uint32_t s_row, const uint32_t s_col,
              const uint32_t row_fst, const uint32_t col_fst)
{{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s_col; i++){{
        for (j = 0; j < s_row; j++){{
            mtx[(i+row_fst)*dim+j+col_fst] = blk[i*s_row+j];
        }}
    }}
}}

void getBlock(real_t blk[],
              const real_t mtx[], const uint32_t dim,
              const uint32_t row_fst, const uint32_t col_fst,
              const uint32_t s_row, const uint32_t s_col)
{{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s_col; i++){{
        for (j = 0; j < s_row; j++){{
            blk[i*s_row+j] = mtx[(i+row_fst)*dim+j+col_fst];
        }}
    }}
}}
