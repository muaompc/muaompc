#include <stdio.h>

#include "mpcctl.h"
#include "mpcmydatctldata.h"

int main(int argc, char *argv[]) {
    struct mpc_ctl ctlst;  /* Structure for static memory allocation */
    struct mpc_ctl *ctl;  /* pointer to the an allocated structure */

    ctl = &ctlst;
    mpcmydat_ctl_setup_ctl(ctl);

    ctl->parameters->x_bar[0] = 0.1;
    ctl->parameters->x_bar[1] = -0.5; 
    ctl->solver->conf->in_iter = 10;
    mpc_ctl_solve_problem(ctl);
    printf("ctl->u_opt[0] = %e \n", ctl->u_opt[0]);

    return 0;
}
