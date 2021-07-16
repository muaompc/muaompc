#include <stdlib.h>  /* malloc */

#include "{prefix}ctldynmem.h"


/* Extern function definitions */

struct {prefix}_ctl *{prefix}_ctl_allocate_ctl(void) 
{{
        struct {prefix}_ctl *ctl = (struct {prefix}_ctl*)malloc(sizeof(struct {prefix}_ctl));
        if (NULL == ctl) {{return {prefix}_ctl_free_ctl(ctl);}}
        ctl->former = NULL;
        ctl->parameters = NULL;
        ctl->solver = NULL;
        ctl->former = {prefix}_{former}_allocate_former();
        if (NULL == ctl->former) {{return {prefix}_ctl_free_ctl(ctl);}}
        ctl->parameters = {prefix}_{former}_allocate_parameters();
        if (NULL == ctl->parameters) {{return {prefix}_ctl_free_ctl(ctl);}}
        ctl->solver = {prefix}_{solver}_allocate_solver();
        if (NULL == ctl->solver) {{return {prefix}_ctl_free_ctl(ctl);}}

    return ctl;
}}

{prefix}_dynmem_error_t {prefix}_ctl_setup_ctl(struct {prefix}_ctl *ctl, char *fname)
{{
        {prefix}_dynmem_error_t ret;

        ret = {prefix}_{former}_setup_former(ctl->former, fname);
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
        ret = {prefix}_{former}_setup_parameters(ctl->parameters, ctl->former);
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}
        ret = {prefix}_{solver}_setup_solver(ctl->solver, ctl->former->prb, fname);
        if ({PREFIX}_DYNMEM_OK != ret) {{return ret;}}

        ctl->u_opt = ctl->solver->u_opt;

        return {PREFIX}_DYNMEM_OK;
}}

struct {prefix}_ctl *{prefix}_ctl_free_ctl(struct {prefix}_ctl *ctl) {{
    if (ctl) {{
        {prefix}_{former}_free_former(ctl->former);
        {prefix}_{former}_free_parameters(ctl->parameters);
        /* {prefix}_{solver}_free_solver(ctl->solver); TODO uncomment*/
        {prefix}_free(ctl);
    }}
    return NULL;
}}

