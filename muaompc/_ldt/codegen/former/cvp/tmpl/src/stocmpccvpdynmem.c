#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */
#include "include/cjson.h"
#include "include/stocmpcdynmem.h"
#include "include/stocmpccvpdynmem.h"

/* Static functions declarations */
static stocmpc_dynmem_error_t stocmpc_get_json_pmetric(struct stocmpc_pmetric *pmetric, cJSON *data, char *jname, char *term_name);
static stocmpc_dynmem_error_t stocmpc_get_json_pmetric_facs(struct stocmpc_pmetric *pmetric, cJSON *jpmetric);
static stocmpc_dynmem_error_t stocmpc_get_json_pmetric_fac_term(struct stocmpc_term *term, cJSON *jpmetric, int fac_pos);
static stocmpc_dynmem_error_t stocmpc_get_json_pmetric_sub_term(struct stocmpc_term *term, cJSON *jpmetric, char *sub_name);
static stocmpc_dynmem_error_t stocmpc_get_json_term(struct stocmpc_term *term, cJSON *data, char *jname, char *term_name);
static stocmpc_dynmem_error_t stocmpc_parse_elements(struct stocmpc_cvp *cvp, cJSON *data);
static stocmpc_dynmem_error_t stocmpc_get_json_term_items(struct stocmpc_term *term, cJSON *jobj);
static stocmpc_dynmem_error_t stocmpc_alloc_data(real_t **data, int elems);
static void stocmpc_free_pmetric(struct stocmpc_pmetric* p);
static void stocmpc_free_term(struct stocmpc_term* term);

/* Extern function definitions */

struct stocmpc_cvp *stocmpc_cvp_allocate_former(void) {
    int i;
    struct stocmpc_cvp *cvp = (struct stocmpc_cvp*)malloc(sizeof(struct stocmpc_cvp));
    if (NULL == cvp) {return stocmpc_cvp_free_former(cvp);}

    for (i=0; i<STOCMPC_PAR_NUM; i++) {
        cvp->par[i] = NULL;
    }
    for (i=0; i<STOCMPC_CONSTANT_NUM; i++) {
        cvp->constant[i] = NULL;
    }
    for (i=0; i<STOCMPC_PMETRIC_NUM; i++) {
        cvp->pmetric[i] = NULL;
    }

    /* parameters */
    for (i=0; i<STOCMPC_PAR_NUM; i++) {
        cvp->par[i] = (struct stocmpc_term*)malloc(sizeof(struct stocmpc_term));
        if (NULL == cvp->par[i]) {return stocmpc_cvp_free_former(cvp);}
    }

    /* constants */
    for (i=0; i<STOCMPC_CONSTANT_NUM; i++) {
        cvp->constant[i] = (struct stocmpc_term*)malloc(sizeof(struct stocmpc_term));
        if (NULL == cvp->constant[i]) {return stocmpc_cvp_free_former(cvp);}
    }

    /* parametric */
    for (i=0; i<STOCMPC_PMETRIC_NUM; i++) {
        cvp->pmetric[i] = (struct stocmpc_pmetric*)malloc(sizeof(struct stocmpc_pmetric));
        if (NULL == cvp->pmetric[i]) {return stocmpc_cvp_free_former(cvp);}
        cvp->pmetric[i]->val = NULL;
        cvp->pmetric[i]->aux = NULL;
        cvp->pmetric[i]->fac0 = NULL;
        cvp->pmetric[i]->fac = NULL;
        cvp->pmetric[i]->par = NULL;
        cvp->pmetric[i]->par_id = NULL;

        cvp->pmetric[i]->fac_num = (uint32_t*)malloc(sizeof(uint32_t*));
        if (NULL == cvp->pmetric[i]->fac_num) {return stocmpc_cvp_free_former(cvp);}
        cvp->pmetric[i]->fac_num[0] = 0;
        cvp->pmetric[i]->val = (struct stocmpc_term*)malloc(sizeof(struct stocmpc_term));
        if (NULL == cvp->pmetric[i]->val) {return stocmpc_cvp_free_former(cvp);}
        cvp->pmetric[i]->val->data = NULL;
        cvp->pmetric[i]->aux = (struct stocmpc_term*)malloc(sizeof(struct stocmpc_term));
        if (NULL == cvp->pmetric[i]->aux) {return stocmpc_cvp_free_former(cvp);}
        cvp->pmetric[i]->aux->data = NULL;
        cvp->pmetric[i]->fac0 = (struct stocmpc_term*)malloc(sizeof(struct stocmpc_term));
        if (NULL == cvp->pmetric[i]->fac0) {return stocmpc_cvp_free_former(cvp);}
        cvp->pmetric[i]->fac0->data = NULL;
    }

    /* the evaluated problem itself */
    cvp->prb = (struct stocmpc_cvp_prb*)malloc(sizeof(struct stocmpc_cvp_prb));
    if (NULL == cvp->prb) {return stocmpc_cvp_free_former(cvp);}
    cvp->prb->ws1 = cvp->pmetric[STOCMPC_WS1]->val;
    cvp->prb->ws0 = cvp->pmetric[STOCMPC_WS0]->val;
    cvp->prb->wn0 = cvp->pmetric[STOCMPC_WN0]->val;
    cvp->prb->wn1 = cvp->pmetric[STOCMPC_WN1]->val;
    cvp->prb->g = cvp->pmetric[STOCMPC_G]->val;
    cvp->prb->u_lb = cvp->constant[STOCMPC_U_LB];
    cvp->prb->H = cvp->constant[STOCMPC_H];
    cvp->prb->Wm1 = cvp->constant[STOCMPC_WM1];
    cvp->prb->Wm0 = cvp->constant[STOCMPC_WM0];
    cvp->prb->wvT1 = cvp->constant[STOCMPC_WVT1];
    cvp->prb->wvT0 = cvp->constant[STOCMPC_WVT0];
    cvp->prb->u_ub = cvp->constant[STOCMPC_U_UB];


    return cvp;
}

struct stocmpc_cvp *stocmpc_cvp_free_former(struct stocmpc_cvp *cvp) {
    int i;
    if (cvp) {
        if (cvp->par) {
            for (i=0; i<STOCMPC_PAR_NUM; i++) {
                stocmpc_free_term(cvp->par[i]);
            }
        }
        if (cvp->constant) {
            for (i=0; i<STOCMPC_CONSTANT_NUM; i++) {
                stocmpc_free_term(cvp->constant[i]);
            }
        }
        if (cvp->pmetric) {
            for (i=0; i<STOCMPC_PMETRIC_NUM; i++) {
                stocmpc_free_pmetric(cvp->pmetric[i]);
            }
        }
        stocmpc_free(cvp->prb);
        stocmpc_free(cvp);
    }
    return NULL;
}

stocmpc_dynmem_error_t stocmpc_cvp_setup_former(struct stocmpc_cvp *cvp, char *fname) {
    stocmpc_dynmem_error_t ret;
    cJSON *data;
    int i, j, par_id;
    data = stocmpc_dynmem_get_data(fname);
    if (NULL == data) {return STOCMPC_DYNMEM_FAIL;}
    ret = stocmpc_parse_elements(cvp, data);
    cJSON_Delete(data);
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}

    for (i=0; i<STOCMPC_PMETRIC_NUM; i++) {
		for (j=0; j<cvp->pmetric[i]->fac_num[0]; j++) {
			par_id = cvp->pmetric[i]->par_id[j];
			cvp->pmetric[i]->par[j] = cvp->par[par_id];
		}
    }
    return STOCMPC_DYNMEM_OK;
}

struct stocmpc_cvp_parameters *stocmpc_cvp_allocate_parameters(void) {
    struct stocmpc_cvp_parameters *p = (struct stocmpc_cvp_parameters*)
                              malloc(sizeof(struct stocmpc_cvp_parameters));
    if (NULL == p) {return NULL;}
    return p;
}

struct stocmpc_cvp_parameters *stocmpc_cvp_free_parameters(struct stocmpc_cvp_parameters *p) {
    stocmpc_free(p);
    return NULL;
}

stocmpc_dynmem_error_t stocmpc_cvp_setup_parameters(struct stocmpc_cvp_parameters *parameters,
                struct stocmpc_cvp *cvp) {
    parameters->x_k = cvp->par[STOCMPC_X_K]->data;

    return STOCMPC_DYNMEM_OK;
}

/* Static function definitions */
stocmpc_dynmem_error_t stocmpc_parse_elements(struct stocmpc_cvp *cvp, cJSON *data)
{
    stocmpc_dynmem_error_t ret;
    ret = stocmpc_get_json_term(cvp->par[STOCMPC_X_K], data, "par", "x_k");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_pmetric(cvp->pmetric[STOCMPC_WS1], data, "pmetric", "ws1");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_pmetric(cvp->pmetric[STOCMPC_WS0], data, "pmetric", "ws0");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_pmetric(cvp->pmetric[STOCMPC_WN0], data, "pmetric", "wn0");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_pmetric(cvp->pmetric[STOCMPC_WN1], data, "pmetric", "wn1");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_pmetric(cvp->pmetric[STOCMPC_G], data, "pmetric", "g");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_term(cvp->constant[STOCMPC_U_LB], data, "constant", "u_lb");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_term(cvp->constant[STOCMPC_H], data, "constant", "H");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_term(cvp->constant[STOCMPC_WM1], data, "constant", "Wm1");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_term(cvp->constant[STOCMPC_WM0], data, "constant", "Wm0");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_term(cvp->constant[STOCMPC_WVT1], data, "constant", "wvT1");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_term(cvp->constant[STOCMPC_WVT0], data, "constant", "wvT0");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_term(cvp->constant[STOCMPC_U_UB], data, "constant", "u_ub");
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}

    return STOCMPC_DYNMEM_OK;
}

stocmpc_dynmem_error_t stocmpc_get_json_term(struct stocmpc_term *term, cJSON *data, char *jname, char *term_name)
{
    cJSON *jobj, *jterm;

    jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return STOCMPC_DYNMEM_FAIL;
    }
    jterm = cJSON_GetObjectItem(jobj, term_name);
    if (NULL == jterm) {
        printf("ERROR: could not get item %s \n", term_name);
        return STOCMPC_DYNMEM_FAIL;
    }
    stocmpc_get_json_term_items(term, jterm);
    if (NULL == term) {return STOCMPC_DYNMEM_FAIL;}

    return STOCMPC_DYNMEM_OK;
 }

stocmpc_dynmem_error_t stocmpc_get_json_pmetric(struct stocmpc_pmetric *pmetric,
		cJSON *data, char *jname, char *term_name)
{
    stocmpc_dynmem_error_t ret;
    cJSON *jobj, *jpmetric;

    jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return STOCMPC_DYNMEM_FAIL;
    }
    jpmetric = cJSON_GetObjectItem(jobj, term_name);
    if (NULL == jpmetric) {
        printf("ERROR: could not get item %s \n", term_name);
        return STOCMPC_DYNMEM_FAIL;
    }

	ret = stocmpc_get_json_pmetric_facs(pmetric, jpmetric);
	if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_pmetric_sub_term(pmetric->val, jpmetric, "val");
	if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_pmetric_sub_term(pmetric->aux, jpmetric, "aux");
	if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    ret = stocmpc_get_json_pmetric_sub_term(pmetric->fac0, jpmetric, "fac0");
	if (STOCMPC_DYNMEM_OK != ret) {return ret;}

    return STOCMPC_DYNMEM_OK;
}

stocmpc_dynmem_error_t stocmpc_get_json_pmetric_facs(struct stocmpc_pmetric *pmetric,
		cJSON *jpmetric)
{
    stocmpc_dynmem_error_t ret;
    cJSON *jfacnum, *jid_list;
    int i;

    jfacnum = cJSON_GetObjectItem(jpmetric, "fac_num");
    if (NULL == jfacnum) {
        printf("ERROR: could not get item %s \n", "fac_num");
        return STOCMPC_DYNMEM_FAIL;
    }

    pmetric->fac_num[0] = (uint32_t)jfacnum->valueint;

    pmetric->par_id = (uint32_t*)calloc(pmetric->fac_num[0], sizeof(uint32_t));
    if (NULL == pmetric->par_id) {return STOCMPC_DYNMEM_FAIL;}
    pmetric->fac = (struct stocmpc_term**)calloc(pmetric->fac_num[0], sizeof(struct stocmpc_term*));
    if (NULL == pmetric->fac) {return STOCMPC_DYNMEM_FAIL;}
    pmetric->par = (struct stocmpc_term**)calloc(pmetric->fac_num[0], sizeof(struct stocmpc_term*));
    if (NULL == pmetric->par) {return STOCMPC_DYNMEM_FAIL;}

    jid_list = cJSON_GetObjectItem(jpmetric, "par_id");
    for (i=0; i<pmetric->fac_num[0]; i++) {
        pmetric->par_id[i] = (uint32_t)cJSON_GetArrayItem(jid_list, i)->valueint;

        pmetric->fac[i] = (struct stocmpc_term*)calloc(pmetric->fac_num[0], sizeof(struct stocmpc_term));
        if (NULL == pmetric->fac[i]) {return STOCMPC_DYNMEM_FAIL;}
        ret = stocmpc_get_json_pmetric_fac_term(pmetric->fac[i], jpmetric, i);
        if (STOCMPC_DYNMEM_OK != ret) {return ret;}
    }

    return STOCMPC_DYNMEM_OK;
}

stocmpc_dynmem_error_t stocmpc_get_json_pmetric_fac_term(struct stocmpc_term *term,
		cJSON *jpmetric, int fac_pos)
{
    stocmpc_dynmem_error_t ret;
    cJSON *jfac_list, *jfac;

    jfac_list = cJSON_GetObjectItem(jpmetric, "fac");
    if (NULL == jfac_list) {
        printf("ERROR: could not get item %s \n", "fac");
        return STOCMPC_DYNMEM_FAIL;
    }
    jfac = cJSON_GetArrayItem(jfac_list, fac_pos);
    if (NULL == jfac) {
        printf("ERROR: could not get array item in position %d \n", fac_pos);
        return STOCMPC_DYNMEM_FAIL;
    }

    ret = stocmpc_get_json_term_items(term, jfac);
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}

    return STOCMPC_DYNMEM_OK;
 }

stocmpc_dynmem_error_t stocmpc_get_json_pmetric_sub_term(struct stocmpc_term *term,
		cJSON *jpmetric, char *sub_name)
{
    stocmpc_dynmem_error_t ret;
    cJSON *jterm;

    jterm = cJSON_GetObjectItem(jpmetric, sub_name);
    if (NULL == jterm) {
        printf("ERROR: could not get item %s \n", sub_name);
        return STOCMPC_DYNMEM_FAIL;
    }
    ret = stocmpc_get_json_term_items(term, jterm);
    if (STOCMPC_DYNMEM_OK != ret) {return ret;}

    return STOCMPC_DYNMEM_OK;
 }

stocmpc_dynmem_error_t stocmpc_get_json_term_items(struct stocmpc_term *term, cJSON *jobj)
{
    int elems;
    uint32_t i;
    cJSON *c, *jdata;

    stocmpc_dynmem_error_t ret;
    c = cJSON_GetObjectItem(jobj, "cols");
    if (NULL == c) {
        printf("ERROR: could not get item %s \n", "cols");
        return STOCMPC_DYNMEM_FAIL;
    }
    term->cols = (uint32_t)c->valueint;
    c = cJSON_GetObjectItem(jobj, "rows");
    if (NULL == c) {
        printf("ERROR: could not get item %s \n", "rows");
        return STOCMPC_DYNMEM_FAIL;
    }
    term->rows = (uint32_t)c->valueint;
    jdata = cJSON_GetObjectItem(jobj, "data");
    if (NULL == jdata) {
        printf("ERROR: could not get item %s \n", "data");
        return STOCMPC_DYNMEM_FAIL;
    }
    elems = (int)(term->rows * term->cols);
    if (elems != cJSON_GetArraySize(jdata)) {
        printf("Size error, expected: %d*%d (rows*cols); actual: %d \n",
            term->rows, term->cols, cJSON_GetArraySize(jdata));
        {return STOCMPC_DYNMEM_FAIL;}
    } else {
        ret = stocmpc_alloc_data(&(term->data), elems);
        if (STOCMPC_DYNMEM_OK != ret) {return ret;}
        for (i=0;i<elems;i++) {
            term->data[i] = (real_t) cJSON_GetArrayItem(jdata, i)->valuedouble;
        }
    }

    return STOCMPC_DYNMEM_OK;
 }

stocmpc_dynmem_error_t stocmpc_alloc_data(real_t **data, int elems)
{
    data[0] = (real_t*) malloc(elems*sizeof(real_t));
    if (NULL == data[0]) {return STOCMPC_DYNMEM_FAIL;}

    return STOCMPC_DYNMEM_OK;
}

void stocmpc_free_term(struct stocmpc_term* term)
{
    if (term) {
        stocmpc_free(term->data);
    }
    stocmpc_free(term);
    return;
}

void stocmpc_free_pmetric(struct stocmpc_pmetric* p)
{
    int i;
    if (p) {
        if (p->fac_num) {
            if (p->fac) {
                for (i=0; i < p->fac_num[0]; i++) {
                    stocmpc_free_term(p->fac[i]);
                }
                stocmpc_free(p->fac);
                stocmpc_free(p->par);
            }
            stocmpc_free(p->fac_num);
        }
        stocmpc_free(p->par_id);
        stocmpc_free_term(p->val);
        stocmpc_free_term(p->aux);
        stocmpc_free_term(p->fac0);
        stocmpc_free(p);
    }
    return;
}

