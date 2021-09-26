#include "mpcctl.h"
#include "mpcmydatctldata.h"

struct mpc_ctl ctlst, *ctl;
unsigned long t1_cpu;
unsigned long dt_cpu;

void setup() {
    Serial.begin(9600);
    /* muaompc configuration */
    ctl = &ctlst;
    mpcmydat_ctl_setup_ctl(ctl);
    ctl->solver->conf->in_iter = 2;
    ctl->solver->conf->warm_start = 1; 
}

void loop() {
    /* Solve MPC problem and print the first element of input sequence */
    /* Initialize your parameter (i.e. the current state) */
    ctl->parameters->x_bar[0] = 10.;
    ctl->parameters->x_bar[1] = -10.; 
    t1_cpu = micros();
    mpc_ctl_solve_problem(ctl);  /* solve the MPC problem */
    dt_cpu = micros()-t1_cpu;
    Serial.print("u[0] = "); Serial.print(ctl->u_opt[0]); Serial.println();
    Serial.print("computation time in microseconds = "); Serial.print(dt_cpu);
    Serial.println();
    delay(1000);  /* Do not print too often */
}
