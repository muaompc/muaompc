#define USE_DYNMEM  /* Delete this line to use static memory allocation */

#include <stdio.h>  /* printf, NULL */

#include "{prefix}ctl.h"  /* struct {prefix}_ctl */
#ifdef USE_DYNMEM
#include "{prefix}ctldynmem.h"  /* {prefix}_ctl_setup_ctl and friends */
#else  /* use static memory allocation */
#include "{prefix}{dataname}ctldata.h"  /* {prefix}{dataname}_ctl_setup_ctl */
#endif

int main(int argc, char *argv[]) {{
    struct {prefix}_ctl *ctl;  /* pointer to the an allocated structure */
#ifndef USE_DYNMEM
    struct {prefix}_ctl ctlst;  /* Structure for static memory allocation */
#endif

#ifdef USE_DYNMEM  /* dynamic memory allocation for structure ctl */
    {prefix}_dynmem_error_t err;
    char *fname = "{prefix}{dataname}.json";

    ctl = {prefix}_ctl_allocate_ctl();
    if (NULL == ctl) {{
	printf("ERROR: could not allocate ctl \n");
	return 0;
    }}

    err = {prefix}_ctl_setup_ctl(ctl, fname);
    if (err) {{
	printf("ERROR: could not setup ctl \n");
	return 0;
    }}
#else  /* use static memory allocation */
    ctl = &ctlst;
    {prefix}{dataname}_ctl_setup_ctl(ctl);
#endif

    /* Before calling this function, setup your parameters */
    /* ctl->parameters->xxx[0] = 0.; */ /* parameters in the problem file */
    ctl->solver->conf->in_iter = 1;
    /* ctl->solver->conf->ex_iter = 1; */ /* if state constraints present */ 
    ctl->solver->conf->warm_start = 1; 
    {prefix}_ctl_solve_problem(ctl);
    /* ctl->u_opt */ /* u_opt contains the computed control input */

#ifdef USE_DYNMEM  /* free the memory dynamically allocated for ctl */
    {prefix}_ctl_free_ctl(ctl);
#endif

    return 0;
}}
