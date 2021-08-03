#include <stdio.h>

#include "mpcctl.h"
#include "mpcctldynmem.h"

int main(int argc, char *argv[]) {
    struct mpc_ctl ctlst;  /* Structure for static memory allocation */
    struct mpc_ctl *ctl;  /* pointer to the an allocated structure */

    mpc_dynmem_error_t err;
    char *fname = "mpcmydat.json";

    ctl = mpc_ctl_allocate_ctl();
    if (NULL == ctl) {
	printf("ERROR: could not allocate ctl \n");
	return 0;
    }

    err = mpc_ctl_setup_ctl(ctl, fname);
    if (err) {
	printf("ERROR: could not setup ctl \n");
	return 0;
    }

    ctl->parameters->x_bar[0] = 0.1;
    ctl->parameters->x_bar[1] = -0.5; 
    ctl->solver->conf->in_iter = 1;
    ctl->solver->conf->warm_start = 1; 
    mpc_ctl_solve_problem(ctl);
    printf("ctl->u_opt[0] = %e \n", ctl->u_opt[0]);

    mpc_ctl_free_ctl(ctl);

    return 0;
}
