#include <stddef.h>  /* NULL */
#include "{prefix}cvp.h"
{alloc_terms}

{alloc_socc_terms}

void {prefix}{dataname}_cvp_setup_former(struct {prefix}_cvp *cvp)
{{
{init_terms}
    cvp->socc_num = &{pname}_socc_num;
    cvp->socc = {pname}_socc;
{init_socc_terms}

    cvp->prb = &{pname}_prb;
{init_prb}
    return;
}}

void {prefix}{dataname}_cvp_setup_parameters(struct {prefix}_cvp_parameters *parameters,
    struct {prefix}_cvp *cvp)
{{
{init_parameters}
    return;
}}
