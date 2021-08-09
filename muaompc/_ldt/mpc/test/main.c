#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#include "mpcalmctl.h"
#include "mpcalmregmpcctldata.h"

int main(int argc, char *argv[])
{
    struct mpcalm_ctl *ctl, ctlst;

  if (8 != argc) {
    printf("[0., 0., 0., 0., 0.]\n");
  } else {
    ctl = &ctlst;
    mpcalmregmpc_ctl_setup_ctl(ctl);
    ctl->parameters->x_k[0] = (real_t)strtod(argv[1], NULL);
    ctl->parameters->x_k[1] = (real_t)strtod(argv[2], NULL);
    ctl->solver->u_ini[0] = (real_t)strtod(argv[3], NULL);
    ctl->solver->u_ini[1] = (real_t)strtod(argv[4], NULL);
    ctl->solver->l_ini[0] = (real_t)strtod(argv[5], NULL);
    ctl->solver->l_ini[1] = (real_t)strtod(argv[6], NULL);
    ctl->solver->l_ini[2] = (real_t)strtod(argv[7], NULL);
    ctl->solver->conf->in_iter = 2;
    ctl->solver->conf->ex_iter = 3;
    ctl->solver->conf->warm_start = 0;
    mpcalm_ctl_solve_problem(ctl);
    printf("[%e, %e, %e, %e, %e]\n", ctl->u_opt[0], ctl->u_opt[1],
        ctl->solver->l_opt[0], ctl->solver->l_opt[1], ctl->solver->l_opt[2]);
  }

  return 0;
}


