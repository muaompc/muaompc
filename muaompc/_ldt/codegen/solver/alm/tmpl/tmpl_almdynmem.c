#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */

#include "cjson.h"

#include "{prefix}cvp.h"
#include "{prefix}fgmdynmem.h"
#include "{prefix}almdynmem.h"

/* Static functions declarations */
{prefix}_dynmem_error_t {prefix}_alm_parse_elements(struct {prefix}_alm *alm, cJSON *data);

/* Extern function definitions */

struct {prefix}_alm *{prefix}_alm_allocate_solver(void)
{{
    struct {prefix}_alm *alm = (struct {prefix}_alm*)malloc(sizeof(struct {prefix}_alm));
    if (NULL == alm) {{return {prefix}_alm_free_solver(alm);}}
    alm->conf = NULL;
    alm->fgm = NULL;
    alm->mu = NULL;
    alm->Linv = NULL;
    alm->u_opt = NULL;
    alm->u_ini = NULL;
    alm->l_opt = NULL;
    alm->l_ini = NULL;
    alm->tmp1_optvar_seqlen = NULL;
    alm->tmp2_optvar_seqlen = NULL;
    alm->tmp3_optvar_seqlen = NULL;
    alm->tmp4_optvar_seqlen = NULL;
    alm->tmp5_optvar_seqlen = NULL;
    alm->tmp6_optvar_seqlen = NULL;
    alm->tmp7_optvar_seqlen = NULL;
    alm->tmp1_lagmul_seqlen = NULL;
    alm->tmp2_lagmul_seqlen = NULL;
    alm->tmp3_lagmul_seqlen = NULL;
    alm->tmp4_lagmul_seqlen = NULL;
    alm->tmp5_lagmul_seqlen = NULL;
    alm->tmp6_lagmul_seqlen = NULL;
    alm->tmp7_lagmul_seqlen = NULL;
    alm->tmp1_lagmul_optvar_seqlen = NULL;

    alm->conf = (struct {prefix}_alm_conf*)malloc(sizeof(struct {prefix}_alm_conf));
    if (NULL == alm->conf) {{return {prefix}_alm_free_solver(alm);}}

    alm->mu = (real_t *)malloc(sizeof(real_t));
    if (NULL == alm->mu) {{return {prefix}_alm_free_solver(alm);}}
    alm->Linv = (real_t *)malloc(sizeof(real_t));
    if (NULL == alm->Linv) {{return {prefix}_alm_free_solver(alm);}}

    alm->fgm = {prefix}_fgm_allocate_solver();
    if (NULL == alm->fgm) {{return {prefix}_alm_free_solver(alm);}}

    return alm;
}}

{prefix}_dynmem_error_t {prefix}_alm_setup_solver(struct {prefix}_alm *alm, 
                struct {prefix}_cvp_prb *prb, char *fname)
{{
    uint32_t i;
    {prefix}_dynmem_error_t ret;
    cJSON *data;

    data = {prefix}_dynmem_get_data(fname);
    if (NULL == data) {{return {PREFIX}_DYNMEM_FAIL;}}
    ret = {prefix}_alm_parse_elements(alm, data);
    cJSON_Delete(data);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_fgm_setup_solver(alm->fgm, prb, fname);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    alm->V = prb->V->data;
    alm->v_lb = prb->v_lb->data;
    alm->v_ub = prb->v_ub->data;
    alm->optvar_seqlen = prb->u_lb->rows;
    alm->lagmul_seqlen = prb->v_lb->rows;
    alm->sizeof_optvar_seqlen = sizeof(real_t) * prb->u_lb->rows;
    alm->sizeof_lagmul_seqlen = sizeof(real_t) * prb->v_lb->rows;

    alm->conf->in_iter = 1;
    alm->conf->ex_iter = 1;
    alm->conf->warm_start = 0;
    alm->i_ex = &(alm->conf->ex_iter);
    alm->fgm->j_in = &(alm->conf->in_iter);

    alm->u_ini = alm->fgm->u_ini;
    alm->u_opt = alm->fgm->u_opt;
    alm->l_ini = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->l_ini) {{return {PREFIX}_DYNMEM_FAIL;}}
    for (i=0; i<alm->lagmul_seqlen; i++) {{
	    alm->l_ini[i] = (real_t)0.;
    }}
    alm->l_opt = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->l_opt) {{return {PREFIX}_DYNMEM_FAIL;}}


    alm->tmp1_optvar_seqlen = (real_t *)malloc(alm->sizeof_optvar_seqlen);
    if (NULL == alm->tmp1_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp2_optvar_seqlen = (real_t *)malloc(alm->sizeof_optvar_seqlen);
    if (NULL == alm->tmp2_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp3_optvar_seqlen = (real_t *)malloc(alm->sizeof_optvar_seqlen);
    if (NULL == alm->tmp3_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp4_optvar_seqlen = (real_t *)malloc(alm->sizeof_optvar_seqlen);
    if (NULL == alm->tmp4_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp5_optvar_seqlen = (real_t *)malloc(alm->sizeof_optvar_seqlen);
    if (NULL == alm->tmp5_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp6_optvar_seqlen = (real_t *)malloc(alm->sizeof_optvar_seqlen);
    if (NULL == alm->tmp6_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp7_optvar_seqlen = (real_t *)malloc(alm->sizeof_optvar_seqlen);
    if (NULL == alm->tmp7_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp1_lagmul_seqlen = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->tmp1_lagmul_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp2_lagmul_seqlen = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->tmp2_lagmul_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp3_lagmul_seqlen = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->tmp3_lagmul_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp4_lagmul_seqlen = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->tmp4_lagmul_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp5_lagmul_seqlen = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->tmp5_lagmul_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp6_lagmul_seqlen = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->tmp6_lagmul_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp7_lagmul_seqlen = (real_t *)malloc(alm->sizeof_lagmul_seqlen);
    if (NULL == alm->tmp7_lagmul_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    alm->tmp1_lagmul_optvar_seqlen = (real_t *)malloc(
                    alm->sizeof_lagmul_seqlen * alm->sizeof_optvar_seqlen);
    if (NULL == alm->tmp1_lagmul_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}

    return {PREFIX}_DYNMEM_OK;
}}

struct {prefix}_alm *{prefix}_alm_free_solver(struct {prefix}_alm *alm) {{
    if (alm) {{
        {prefix}_free(alm->tmp1_lagmul_optvar_seqlen);
        {prefix}_free(alm->tmp1_lagmul_seqlen);
        {prefix}_free(alm->tmp2_lagmul_seqlen);
        {prefix}_free(alm->tmp3_lagmul_seqlen);
        {prefix}_free(alm->tmp4_lagmul_seqlen);
        {prefix}_free(alm->tmp5_lagmul_seqlen);
        {prefix}_free(alm->tmp6_lagmul_seqlen);
        {prefix}_free(alm->tmp7_lagmul_seqlen);
        {prefix}_free(alm->tmp1_optvar_seqlen);
        {prefix}_free(alm->tmp2_optvar_seqlen);
        {prefix}_free(alm->tmp3_optvar_seqlen);
        {prefix}_free(alm->tmp4_optvar_seqlen);
        {prefix}_free(alm->tmp5_optvar_seqlen);
        {prefix}_free(alm->tmp6_optvar_seqlen);
        {prefix}_free(alm->tmp7_optvar_seqlen);
        {prefix}_free(alm->l_ini);
        {prefix}_free(alm->l_opt);
        {prefix}_free(alm->Linv);
        {prefix}_free(alm->mu);
        {prefix}_free(alm->conf);
        {prefix}_fgm_free_solver(alm->fgm);
        {prefix}_free(alm);
    }}
    return NULL;
}}

/* Static function definitions */

{prefix}_dynmem_error_t {prefix}_alm_parse_elements(struct {prefix}_alm *alm,
                cJSON *data)
{{
    cJSON *c, *lagmul, *veclen;

    c = cJSON_GetObjectItem(data, "mu");
    if (NULL == c) {{
        printf("ERROR: could not parse item %s \n", "mu");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    *(alm->mu) = (real_t)c->valuedouble;

    c = cJSON_GetObjectItem(data, "Linv");
    if (NULL == c) {{
        printf("ERROR: could not parse item %s \n", "Linv");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    *(alm->Linv) = (real_t)c->valuedouble;

    lagmul = cJSON_GetObjectItem(data, "lagmul");
    if (NULL == lagmul) {{
        printf("ERROR: could not parse item %s \n", "lagmul");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    veclen = cJSON_GetObjectItem(lagmul, "veclen");
    if (NULL == veclen) {{
        printf("ERROR: could not parse item %s \n", "veclen");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    alm->lagmul_veclen = (uint32_t)veclen->valueint;

    return {PREFIX}_DYNMEM_OK;
}}

