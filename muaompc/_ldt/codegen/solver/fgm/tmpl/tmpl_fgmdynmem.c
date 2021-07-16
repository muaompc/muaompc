#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */

#include "cjson.h"

#include "{prefix}dynmem.h"
#include "{prefix}cvp.h"
#include "{prefix}fgmdynmem.h"

/* Static functions declarations */
static {prefix}_dynmem_error_t {prefix}_fgm_parse_elements(struct {prefix}_fgm *fgm, cJSON *data);

/* Extern function definitions */

struct {prefix}_fgm *{prefix}_fgm_allocate_solver(void)
{{
    struct {prefix}_fgm *fgm = (struct {prefix}_fgm*)malloc(sizeof(struct {prefix}_fgm));
    if (NULL == fgm) {{return {prefix}_fgm_free_solver(fgm);}}
    fgm->conf = NULL;
    fgm->nu = NULL;
    fgm->u_opt = NULL;
    fgm->u_ini = NULL;
    fgm->tmp1_optvar_seqlen = NULL;
    fgm->tmp2_optvar_seqlen = NULL;
    fgm->tmp3_optvar_seqlen = NULL;
    fgm->tmp4_optvar_seqlen = NULL;
    fgm->tmp5_optvar_seqlen = NULL;
    fgm->tmp6_optvar_seqlen = NULL;

    fgm->conf = (struct {prefix}_fgm_conf*)malloc(sizeof(struct {prefix}_fgm_conf));
    if (NULL == fgm->conf) {{return {prefix}_fgm_free_solver(fgm);}}

    fgm->nu = (real_t *)malloc(sizeof(real_t));
    if (NULL == fgm->nu) {{return {prefix}_fgm_free_solver(fgm);}}

    return fgm;
}}

{prefix}_dynmem_error_t {prefix}_fgm_setup_solver(struct {prefix}_fgm *fgm, struct {prefix}_cvp_prb *prb, char *fname)
{{
    uint32_t i;
    {prefix}_dynmem_error_t ret;
    cJSON *data;

    data = {prefix}_dynmem_get_data(fname);
    if (NULL == data) {{return {PREFIX}_DYNMEM_FAIL;}}
    ret = {prefix}_fgm_parse_elements(fgm, data);
    cJSON_Delete(data);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    fgm->goL = prb->g->data;
    fgm->HoL = prb->H->data;
    fgm->u_lb = prb->u_lb->data;
    fgm->u_ub = prb->u_ub->data;
    fgm->optvar_seqlen = prb->u_lb->rows;
    fgm->sizeof_optvar_seqlen = sizeof(real_t) * prb->u_lb->rows;

    fgm->conf->in_iter = 1;
    fgm->conf->warm_start = 0;
    fgm->j_in = &(fgm->conf->in_iter);

    fgm->u_ini = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->u_ini) {{return {PREFIX}_DYNMEM_FAIL;}}
    for (i=0; i<fgm->optvar_seqlen; i++) {{
	    fgm->u_ini[i] = (real_t)0.;
    }}
    fgm->u_opt = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->u_opt) {{return {PREFIX}_DYNMEM_FAIL;}}

    fgm->tmp1_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp1_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    fgm->tmp2_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp2_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    fgm->tmp3_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp3_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    fgm->tmp4_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp4_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    fgm->tmp5_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp5_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}
    fgm->tmp6_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp6_optvar_seqlen) {{return {PREFIX}_DYNMEM_FAIL;}}

    return {PREFIX}_DYNMEM_OK;
}}

struct {prefix}_fgm *{prefix}_fgm_free_solver(struct {prefix}_fgm *fgm) {{
    if (fgm) {{
        {prefix}_free(fgm->tmp1_optvar_seqlen);
        {prefix}_free(fgm->tmp2_optvar_seqlen);
        {prefix}_free(fgm->tmp3_optvar_seqlen);
        {prefix}_free(fgm->tmp4_optvar_seqlen);
        {prefix}_free(fgm->tmp5_optvar_seqlen);
        {prefix}_free(fgm->tmp6_optvar_seqlen);
        {prefix}_free(fgm->u_ini);
        {prefix}_free(fgm->u_opt);
        {prefix}_free(fgm->nu);
        {prefix}_free(fgm->conf);
        {prefix}_free(fgm);
    }}
    return NULL;
}}

/* Static function definitions */

{prefix}_dynmem_error_t {prefix}_fgm_parse_elements(struct {prefix}_fgm *fgm, cJSON *data)
{{
    cJSON *nu, *optvar, *veclen;

    nu = cJSON_GetObjectItem(data, "nu");
    if (NULL == nu) {{
        printf("ERROR: could not parse item %s \n", "nu");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    *(fgm->nu) = (real_t)nu->valuedouble;

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
    fgm->optvar_veclen = (uint32_t)veclen->valueint;

    return {PREFIX}_DYNMEM_OK;
}}

