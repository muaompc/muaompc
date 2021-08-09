#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#include "stocmpcctl.h"
#include "stocmpccvp.h"
#include "stocmpcstocmpcctldata.h"

int main(int argc, char *argv[])
{
    struct stocmpc_ctl *ctl, ctlst;
    struct stocmpc_cvp_prb *prb;
    uint32_t i, j;

  if (5 != argc) {
    printf("[0., 0.]\n");
  } else {
    ctl = &ctlst;
    stocmpcstocmpc_ctl_setup_ctl(ctl);
    ctl->parameters->x_k[0] = (real_t)strtod(argv[1], NULL);
    ctl->parameters->x_k[1] = (real_t)strtod(argv[2], NULL);
    ctl->solver->z_ini[0] = (real_t)strtod(argv[3], NULL);
    ctl->solver->z_ini[1] = (real_t)strtod(argv[4], NULL);
    ctl->solver->conf->in_iter = 2;
    ctl->solver->conf->warm_start = 0;
    stocmpc_ctl_solve_problem(ctl);
    prb = ctl->former->prb;
    printf("[");
    for (j=0; j<prb->socc_num[0]; j++) {
        for (i=0; i<(prb->socc[j]->Wm->rows*prb->socc[j]->Wm->cols); i++) {
            printf("%e, ", prb->socc[j]->Wm->data[i]);
        }
        for (i=0; i<(prb->socc[j]->wn->rows*prb->socc[j]->wn->cols); i++) {
            printf("%e, ", prb->socc[j]->wn->data[i]);
        }
        for (i=0; i<(prb->socc[j]->wvT->rows*prb->socc[j]->wvT->cols); i++) {
            printf("%e, ", prb->socc[j]->wvT->data[i]);
        }
        for (i=0; i<(prb->socc[j]->ws->rows*prb->socc[j]->ws->cols); i++) {
            printf("%e, ", prb->socc[j]->ws->data[i]);
        }
    }
    printf("0.]\n");
  }

  return 0;
}
