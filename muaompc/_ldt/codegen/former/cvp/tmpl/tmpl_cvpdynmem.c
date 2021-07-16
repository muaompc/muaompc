#include <stdio.h>  /* printf */
#include <stdlib.h>  /* malloc */

#include "cjson.h"

#include "{prefix}dynmem.h"
#include "{prefix}cvpdynmem.h"

/* Static functions declarations */
#if {PREFIX}_CVP_PRB_SOCC
static void {prefix}_init_socc_prb(struct {prefix}_cvp_prb_socc *prb_socc, struct {prefix}_cvp_socc *socc);
#endif
static {prefix}_dynmem_error_t {prefix}_get_json_socc(struct {prefix}_cvp *cvp, cJSON *data);
static {prefix}_dynmem_error_t {prefix}_get_json_pmetric(struct {prefix}_pmetric *pmetric, cJSON *data, char *jname, char *term_name);
static {prefix}_dynmem_error_t {prefix}_get_json_pmetric_(struct {prefix}_pmetric *pmetric, cJSON *data, char *term_name);
static {prefix}_dynmem_error_t {prefix}_get_json_pmetric_facs(struct {prefix}_pmetric *pmetric, cJSON *jpmetric);
static {prefix}_dynmem_error_t {prefix}_get_json_pmetric_fac_term(struct {prefix}_term *term, cJSON *jpmetric, uint32_t fac_pos);
static {prefix}_dynmem_error_t {prefix}_get_json_pmetric_sub_term(struct {prefix}_term *term, cJSON *jpmetric, char *sub_name);
static {prefix}_dynmem_error_t {prefix}_get_json_term(struct {prefix}_term *term, cJSON *data, char *jname, char *term_name);
static {prefix}_dynmem_error_t {prefix}_get_json_term_(struct {prefix}_term *term, cJSON *data, char *term_name);
static {prefix}_dynmem_error_t {prefix}_parse_elements(struct {prefix}_cvp *cvp, cJSON *data);
static {prefix}_dynmem_error_t {prefix}_get_json_term_items(struct {prefix}_term *term, cJSON *jobj);
static {prefix}_dynmem_error_t {prefix}_allocate_socc(struct {prefix}_cvp *cvp);
static {prefix}_dynmem_error_t {prefix}_allocate_cvp_socc(struct {prefix}_cvp_socc **socc);
static {prefix}_dynmem_error_t {prefix}_allocate_data(real_t **data, uint32_t elems);
static {prefix}_dynmem_error_t {prefix}_allocate_pmetric(struct {prefix}_pmetric **pmetric);
static void {prefix}_free_pmetric(struct {prefix}_pmetric* p);
static void {prefix}_free_term(struct {prefix}_term* term);
static void {prefix}_free_socc(struct {prefix}_cvp* cvp);
static void {prefix}_free_cvp_socc(struct {prefix}_cvp_socc* socc);

/* Extern function definitions */

/* FIXME: inconsistency, some allocate functions return dynmem_error_t whereas this returns cvp* */

struct {prefix}_cvp *{prefix}_cvp_allocate_former(void) {{
    uint32_t i;
    {prefix}_dynmem_error_t ret;
    struct {prefix}_cvp *cvp = (struct {prefix}_cvp*)malloc(sizeof(struct {prefix}_cvp));
    if (NULL == cvp) {{return {prefix}_cvp_free_former(cvp);}}

    for (i=0; i<{PREFIX}_PAR_NUM; i++) {{
        cvp->par[i] = NULL;
    }}
    for (i=0; i<{PREFIX}_PMETRIC_NUM; i++) {{
        cvp->pmetric[i] = NULL;
    }}

    /* parameters */
    for (i=0; i<{PREFIX}_PAR_NUM; i++) {{
        cvp->par[i] = (struct {prefix}_term*)malloc(sizeof(struct {prefix}_term));
        if (NULL == cvp->par[i]) {{return {prefix}_cvp_free_former(cvp);}}
    }}

    /* parametric */
    for (i=0; i<{PREFIX}_PMETRIC_NUM; i++) {{
        ret = {prefix}_allocate_pmetric(&(cvp->pmetric[i]));
        if ({PREFIX}_DYNMEM_OK != ret) {{return {prefix}_cvp_free_former(cvp);}}
    }}

    /* second-order cone constraints */
    cvp->socc_num = (uint32_t*)malloc(sizeof(uint32_t*));
    if (NULL == cvp->socc_num) {{return {prefix}_cvp_free_former(cvp);}}

    /* the evaluated problem itself */
    cvp->prb = (struct {prefix}_cvp_prb*)malloc(sizeof(struct {prefix}_cvp_prb));
    if (NULL == cvp->prb) {{return {prefix}_cvp_free_former(cvp);}}

    return cvp;
}}

struct {prefix}_cvp *{prefix}_cvp_free_former(struct {prefix}_cvp *cvp) {{
    uint32_t i;
    if (cvp) {{
        if (cvp->par) {{
            for (i=0; i<{PREFIX}_PAR_NUM; i++) {{
                {prefix}_free_term(cvp->par[i]);
            }}
        }}
        if (cvp->pmetric) {{
            for (i=0; i<{PREFIX}_PMETRIC_NUM; i++) {{
                {prefix}_free_pmetric(cvp->pmetric[i]);
            }}
        }}
        {prefix}_free_socc(cvp);
        {prefix}_free(cvp->prb);
        {prefix}_free(cvp);
    }}
    return NULL;
}}

{prefix}_dynmem_error_t {prefix}_cvp_setup_former(struct {prefix}_cvp *cvp, char *fname) {{
    {prefix}_dynmem_error_t ret;
    cJSON *data;
    uint32_t i, j, par_id;
    uint32_t k;
    data = {prefix}_dynmem_get_data(fname);
    if (NULL == data) {{return {PREFIX}_DYNMEM_FAIL;}}
    ret = {prefix}_parse_elements(cvp, data);
    cJSON_Delete(data);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    for (i=0; i<{PREFIX}_PMETRIC_NUM; i++) {{
        for (j=0; j<cvp->pmetric[i]->fac_num[0]; j++) {{
            par_id = cvp->pmetric[i]->par_id[j];
            cvp->pmetric[i]->par[j] = cvp->par[par_id];
        }}
    }}
    for (j=0; j<cvp->socc_num[0]; j++) {{
        for (i=0; i<{PREFIX}_SOCC_PMETRIC_NUM; i++) {{
            for (k=0; k<cvp->socc[j]->pmetric[i]->fac_num[0]; k++) {{
            par_id = cvp->socc[j]->pmetric[i]->par_id[k];
            cvp->socc[j]->pmetric[i]->par[k] = cvp->par[par_id];
            }}
        }}
    }}

{init_struct_prb}

#if {PREFIX}_CVP_PRB_SOCC
    cvp->prb->socc = (struct {prefix}_cvp_prb_socc**)
                              calloc(cvp->socc_num[0], sizeof(struct {prefix}_cvp_prb_socc*));
    for (j=0; j<cvp->socc_num[0]; j++) {{
    cvp->prb->socc[j] = (struct {prefix}_cvp_prb_socc*)
                              calloc(1, sizeof(struct {prefix}_cvp_prb_socc));
        {prefix}_init_socc_prb(cvp->prb->socc[j], cvp->socc[j]);
    }}
#endif
    return {PREFIX}_DYNMEM_OK;
}}

struct {prefix}_cvp_parameters *{prefix}_cvp_allocate_parameters(void) {{
    struct {prefix}_cvp_parameters *p = (struct {prefix}_cvp_parameters*)
                              malloc(sizeof(struct {prefix}_cvp_parameters));
    if (NULL == p) {{return NULL;}}
    return p;
}}

struct {prefix}_cvp_parameters *{prefix}_cvp_free_parameters(struct {prefix}_cvp_parameters *p) {{
    {prefix}_free(p);
    return NULL;
}}

{prefix}_dynmem_error_t {prefix}_cvp_setup_parameters(struct {prefix}_cvp_parameters *parameters,
                struct {prefix}_cvp *cvp) {{
{alloc_parameters}
    return {PREFIX}_DYNMEM_OK;
}}

/* Static function definitions */
#if {PREFIX}_CVP_PRB_SOCC
void {prefix}_init_socc_prb(struct {prefix}_cvp_prb_socc *prb_socc, struct {prefix}_cvp_socc *socc)
{{
    prb_socc->Wm = socc->pmetric[{PREFIX}_SOCC_WM]->val;
    prb_socc->wn = socc->pmetric[{PREFIX}_SOCC_WN]->val;
    prb_socc->wvT = socc->pmetric[{PREFIX}_SOCC_WVT]->val;
    prb_socc->ws = socc->pmetric[{PREFIX}_SOCC_WS]->val;

    return;
}}
#endif

{prefix}_dynmem_error_t {prefix}_parse_elements(struct {prefix}_cvp *cvp, cJSON *data)
{{
    {prefix}_dynmem_error_t ret;
{init_data}
    ret = {prefix}_get_json_pmetric(cvp->pmetric[{PREFIX}_G], data, "pmetric", "g");
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric(cvp->pmetric[{PREFIX}_H], data, "pmetric", "H");
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric(cvp->pmetric[{PREFIX}_U_LB], data, "pmetric", "u_lb");
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric(cvp->pmetric[{PREFIX}_U_UB], data, "pmetric", "u_ub");
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric(cvp->pmetric[{PREFIX}_V], data, "pmetric", "V");
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric(cvp->pmetric[{PREFIX}_V_LB], data, "pmetric", "v_lb");
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric(cvp->pmetric[{PREFIX}_V_UB], data, "pmetric", "v_ub");
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    ret = {prefix}_get_json_socc(cvp, data);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    return {PREFIX}_DYNMEM_OK;
}}

{prefix}_dynmem_error_t {prefix}_get_json_term(struct {prefix}_term *term, cJSON *data, char *jname, char *term_name)
{{
    {prefix}_dynmem_error_t ret;
    cJSON *jobj;

    jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {{
        printf("ERROR: could not get item %s \n", jname);
        return {PREFIX}_DYNMEM_FAIL;
    }}
    ret = {prefix}_get_json_term_(term, jobj, term_name);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    return {PREFIX}_DYNMEM_OK;
 }}

{prefix}_dynmem_error_t {prefix}_get_json_term_(struct {prefix}_term *term, cJSON *jobj,  char *term_name)
{{
    cJSON *jterm;

    jterm = cJSON_GetObjectItem(jobj, term_name);
    if (NULL == jterm) {{
        printf("ERROR: could not get item %s \n", term_name);
        return {PREFIX}_DYNMEM_FAIL;
    }}
    {prefix}_get_json_term_items(term, jterm);
    if (NULL == term) {{return {PREFIX}_DYNMEM_FAIL;}}

    return {PREFIX}_DYNMEM_OK;
 }}

{prefix}_dynmem_error_t {prefix}_get_json_socc(struct {prefix}_cvp *cvp,
		cJSON *data)
{{
    uint32_t i;
    {prefix}_dynmem_error_t ret;
    cJSON *jsocc, *jsoccnum, *jpmetric_list, *jpmetric;
    struct {prefix}_cvp_socc **socc;

    jsocc = cJSON_GetObjectItem(data, "socc");
    if (NULL == jsocc) {{
        printf("ERROR: could not get item %s \n", "socc");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    jsoccnum = cJSON_GetObjectItem(jsocc, "socc_num");
    if (NULL == jsoccnum) {{
        printf("ERROR: could not get item %s \n", "socc_num");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    jpmetric_list = cJSON_GetObjectItem(jsocc, "pmetric");
    if (NULL == jpmetric_list) {{
        printf("ERROR: could not get item %s \n", "pmetric");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    cvp->socc_num[0] = (uint32_t)jsoccnum->valueint;

    ret = {prefix}_allocate_socc(cvp);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    socc = cvp->socc;

    for (i=0; i<cvp->socc_num[0]; i++) {{
        jpmetric = cJSON_GetArrayItem(jpmetric_list, (int)i); 
        if (NULL == jpmetric) {{
            printf("ERROR: could not get %s %s array item in position %d \n", "socc", "pmetric", i);
            return {PREFIX}_DYNMEM_FAIL;
        }}
        ret = {prefix}_get_json_pmetric_(socc[i]->pmetric[{PREFIX}_SOCC_WN], jpmetric, "wn");
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
        ret = {prefix}_get_json_pmetric_(socc[i]->pmetric[{PREFIX}_SOCC_WS], jpmetric, "ws");
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
        ret = {prefix}_get_json_pmetric_(socc[i]->pmetric[{PREFIX}_SOCC_WM], jpmetric, "Wm");
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
        ret = {prefix}_get_json_pmetric_(socc[i]->pmetric[{PREFIX}_SOCC_WVT], jpmetric, "wvT");
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    }}

    return {PREFIX}_DYNMEM_OK;
}}

{prefix}_dynmem_error_t {prefix}_allocate_socc(struct {prefix}_cvp *cvp)
{{
    uint32_t i;
    {prefix}_dynmem_error_t ret;

    cvp->socc = (struct {prefix}_cvp_socc**)calloc(cvp->socc_num[0], sizeof(struct {prefix}_cvp_socc*));
    if (NULL == cvp->socc) {{return {PREFIX}_DYNMEM_FAIL;}}

    for (i=0; i<cvp->socc_num[0]; i++) {{
        ret = {prefix}_allocate_cvp_socc(&(cvp->socc[i]));
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    }}

    return {PREFIX}_DYNMEM_OK;
}}

{prefix}_dynmem_error_t {prefix}_get_json_pmetric(struct {prefix}_pmetric *pmetric,
		cJSON *data, char *jname, char *term_name)
{{
    {prefix}_dynmem_error_t ret;
    cJSON *jobj;

    jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {{
        printf("ERROR: could not get item %s \n", jname);
        return {PREFIX}_DYNMEM_FAIL;
    }}
    ret = {prefix}_get_json_pmetric_(pmetric, jobj, term_name);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    return {PREFIX}_DYNMEM_OK;
}}

{prefix}_dynmem_error_t {prefix}_get_json_pmetric_(struct {prefix}_pmetric *pmetric,
		cJSON *jobj, char *term_name)
{{
    {prefix}_dynmem_error_t ret;
    cJSON *jpmetric;

    jpmetric = cJSON_GetObjectItem(jobj, term_name);
    if (NULL == jpmetric) {{
        printf("ERROR: could not get item %s \n", term_name);
        return {PREFIX}_DYNMEM_FAIL;
    }}
	ret = {prefix}_get_json_pmetric_facs(pmetric, jpmetric);
	if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric_sub_term(pmetric->val, jpmetric, "val");
	if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric_sub_term(pmetric->aux, jpmetric, "aux");
	if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    ret = {prefix}_get_json_pmetric_sub_term(pmetric->fac0, jpmetric, "fac0");
	if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    return {PREFIX}_DYNMEM_OK;
}}

{prefix}_dynmem_error_t {prefix}_get_json_pmetric_facs(struct {prefix}_pmetric *pmetric,
		cJSON *jpmetric)
{{
    {prefix}_dynmem_error_t ret;
    cJSON *jfacnum, *jid_list;
    uint32_t i;

    jfacnum = cJSON_GetObjectItem(jpmetric, "fac_num");
    if (NULL == jfacnum) {{
        printf("ERROR: could not get item %s \n", "fac_num");
        return {PREFIX}_DYNMEM_FAIL;
    }}

    pmetric->fac_num[0] = (uint32_t)jfacnum->valueint;

    pmetric->par_id = (uint32_t*)calloc(pmetric->fac_num[0], sizeof(uint32_t));
    if (NULL == pmetric->par_id) {{return {PREFIX}_DYNMEM_FAIL;}}
    pmetric->fac = (struct {prefix}_term**)calloc(pmetric->fac_num[0], sizeof(struct {prefix}_term*));
    if (NULL == pmetric->fac) {{return {PREFIX}_DYNMEM_FAIL;}}
    pmetric->par = (struct {prefix}_term**)calloc(pmetric->fac_num[0], sizeof(struct {prefix}_term*));
    if (NULL == pmetric->par) {{return {PREFIX}_DYNMEM_FAIL;}}

    jid_list = cJSON_GetObjectItem(jpmetric, "par_id");
    for (i=0; i<pmetric->fac_num[0]; i++) {{
        pmetric->par_id[i] = (uint32_t)cJSON_GetArrayItem(jid_list, (int)i)->valueint;

        pmetric->fac[i] = (struct {prefix}_term*)calloc(pmetric->fac_num[0], sizeof(struct {prefix}_term));
        if (NULL == pmetric->fac[i]) {{return {PREFIX}_DYNMEM_FAIL;}}
        ret = {prefix}_get_json_pmetric_fac_term(pmetric->fac[i], jpmetric, i);
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    }}

    return {PREFIX}_DYNMEM_OK;
}}

{prefix}_dynmem_error_t {prefix}_get_json_pmetric_fac_term(struct {prefix}_term *term,
		cJSON *jpmetric, uint32_t fac_pos)
{{
    {prefix}_dynmem_error_t ret;
    cJSON *jfac_list, *jfac;

    jfac_list = cJSON_GetObjectItem(jpmetric, "fac");
    if (NULL == jfac_list) {{
        printf("ERROR: could not get item %s \n", "fac");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    jfac = cJSON_GetArrayItem(jfac_list, (int)fac_pos);
    if (NULL == jfac) {{
        printf("ERROR: could not get array item in position %d \n", fac_pos);
        return {PREFIX}_DYNMEM_FAIL;
    }}

    ret = {prefix}_get_json_term_items(term, jfac);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    return {PREFIX}_DYNMEM_OK;
 }}

{prefix}_dynmem_error_t {prefix}_get_json_pmetric_sub_term(struct {prefix}_term *term,
		cJSON *jpmetric, char *sub_name)
{{
    {prefix}_dynmem_error_t ret;
    cJSON *jterm;

    jterm = cJSON_GetObjectItem(jpmetric, sub_name);
    if (NULL == jterm) {{
        printf("ERROR: could not get item %s \n", sub_name);
        return {PREFIX}_DYNMEM_FAIL;
    }}
    ret = {prefix}_get_json_term_items(term, jterm);
    if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

    return {PREFIX}_DYNMEM_OK;
 }}

{prefix}_dynmem_error_t {prefix}_get_json_term_items(struct {prefix}_term *term, cJSON *jobj)
{{
    uint32_t elems;
    uint32_t i;
    cJSON *c, *jdata;

    {prefix}_dynmem_error_t ret;
    c = cJSON_GetObjectItem(jobj, "cols");
    if (NULL == c) {{
        printf("ERROR: could not get item %s \n", "cols");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    term->cols = (uint32_t)c->valueint;
    c = cJSON_GetObjectItem(jobj, "rows");
    if (NULL == c) {{
        printf("ERROR: could not get item %s \n", "rows");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    term->rows = (uint32_t)c->valueint;
    jdata = cJSON_GetObjectItem(jobj, "data");
    if (NULL == jdata) {{
        printf("ERROR: could not get item %s \n", "data");
        return {PREFIX}_DYNMEM_FAIL;
    }}
    elems = (term->rows * term->cols);
    if (elems != cJSON_GetArraySize(jdata)) {{
        printf("Size error, expected: %d*%d (rows*cols); actual: %d \n",
            term->rows, term->cols, cJSON_GetArraySize(jdata));
        {{return {PREFIX}_DYNMEM_FAIL;}}
    }} else {{
        ret = {prefix}_allocate_data(&(term->data), elems);
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
        for (i=0; i<elems; i++) {{
            term->data[i] = (real_t) cJSON_GetArrayItem(jdata, (int)i)->valuedouble;
        }}
    }}

    return {PREFIX}_DYNMEM_OK;
 }}

{prefix}_dynmem_error_t {prefix}_allocate_data(real_t **data, uint32_t elems)
{{
    data[0] = (real_t*) malloc(elems*sizeof(real_t));
    if (NULL == data[0]) {{return {PREFIX}_DYNMEM_FAIL;}}

    return {PREFIX}_DYNMEM_OK;
}}

void {prefix}_free_term(struct {prefix}_term* term)
{{
    if (term) {{
        {prefix}_free(term->data);
    }}
    {prefix}_free(term);
    return;
}}

void {prefix}_free_pmetric(struct {prefix}_pmetric* p)
{{
    uint32_t i;
    if (p) {{
        if (p->fac_num) {{
            if (p->fac) {{
                for (i=0; i < p->fac_num[0]; i++) {{
                    {prefix}_free_term(p->fac[i]);
                }}
                {prefix}_free(p->fac);
                {prefix}_free(p->par);
            }}
            {prefix}_free(p->fac_num);
        }}
        {prefix}_free(p->par_id);
        {prefix}_free_term(p->val);
        {prefix}_free_term(p->aux);
        {prefix}_free_term(p->fac0);
        {prefix}_free(p);
    }}
    return;
}}

void {prefix}_free_socc(struct {prefix}_cvp* cvp)
{{
    uint32_t i;
    if (cvp) {{
        if (cvp->socc_num) {{
            if (cvp->socc) {{
                for (i=0; i < cvp->socc_num[0]; i++) {{
                    {prefix}_free_cvp_socc(cvp->socc[i]);
                }}
            {prefix}_free(cvp->socc);
            }}
            {prefix}_free(cvp->socc_num);
        }}
    }}
    return;
}}

void {prefix}_free_cvp_socc(struct {prefix}_cvp_socc *socc) {{
    uint32_t i;
    if (socc) {{
        if (socc->pmetric) {{
            for (i=0; i<{PREFIX}_SOCC_PMETRIC_NUM; i++) {{
                {prefix}_free_pmetric(socc->pmetric[i]);
            }}
        }}
        {prefix}_free(socc);
    }}
    return;
}}

{prefix}_dynmem_error_t {prefix}_allocate_pmetric(struct {prefix}_pmetric **pmetric)
{{
        pmetric[0] = (struct {prefix}_pmetric*)malloc(sizeof(struct {prefix}_pmetric));
        if (NULL == pmetric[0]) {{return {PREFIX}_DYNMEM_FAIL;}}
        pmetric[0]->val = NULL;
        pmetric[0]->aux = NULL;
        pmetric[0]->fac0 = NULL;
        pmetric[0]->fac = NULL;
        pmetric[0]->par = NULL;
        pmetric[0]->par_id = NULL;

        pmetric[0]->fac_num = (uint32_t*)malloc(sizeof(uint32_t*));
        if (NULL == pmetric[0]->fac_num) {{return {PREFIX}_DYNMEM_FAIL;}}
        pmetric[0]->fac_num[0] = 0;
        pmetric[0]->val = (struct {prefix}_term*)malloc(sizeof(struct {prefix}_term));
        if (NULL == pmetric[0]->val) {{return {PREFIX}_DYNMEM_FAIL;}}
        pmetric[0]->val->data = NULL;
        pmetric[0]->aux = (struct {prefix}_term*)malloc(sizeof(struct {prefix}_term));
        if (NULL == pmetric[0]->aux) {{return {PREFIX}_DYNMEM_FAIL;}}
        pmetric[0]->aux->data = NULL;
        pmetric[0]->fac0 = (struct {prefix}_term*)malloc(sizeof(struct {prefix}_term));
        if (NULL == pmetric[0]->fac0) {{return {PREFIX}_DYNMEM_FAIL;}}
        pmetric[0]->fac0->data = NULL;

        return {PREFIX}_DYNMEM_OK;
}}

{prefix}_dynmem_error_t {prefix}_allocate_cvp_socc(struct {prefix}_cvp_socc **socc) {{
    {prefix}_dynmem_error_t ret;
    uint32_t j;

    socc[0] = (struct {prefix}_cvp_socc*)calloc(1, sizeof(struct {prefix}_cvp_socc));
    if (NULL == socc[0]) {{return {PREFIX}_DYNMEM_FAIL;}}
    for (j=0; j<{PREFIX}_SOCC_PMETRIC_NUM; j++) {{
        socc[0]->pmetric[j] = NULL;
    }}
    for (j=0; j<{PREFIX}_SOCC_PMETRIC_NUM; j++) {{
        ret = {prefix}_allocate_pmetric(&(socc[0]->pmetric[j]));
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
    }}

    return {PREFIX}_DYNMEM_OK;
}}
