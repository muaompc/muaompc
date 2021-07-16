#include "{prefix}{dataname}{former}data.h"
#include "{prefix}{dataname}{solver}data.h"
#include "{prefix}ctl.h"

struct {prefix}_{former} {prefix}{dataname}_former;  /**< Problem former. */
struct {prefix}_{solver} {prefix}{dataname}_solver;  /**< Problem solver. */
struct {prefix}_{former}_parameters {prefix}{dataname}_parameters;  /**< Problem parameters. */

/* Extern function definitions */

void {prefix}{dataname}_ctl_setup_ctl(struct {prefix}_ctl *ctl)
{{
	ctl->former = &{prefix}{dataname}_former;
	ctl->solver = &{prefix}{dataname}_solver;
	ctl->parameters = &{prefix}{dataname}_parameters;

	{prefix}{dataname}_{former}_setup_former(ctl->former);
	{prefix}{dataname}_{solver}_setup_solver(ctl->solver, ctl->former->prb);
	{prefix}{dataname}_{former}_setup_parameters(ctl->parameters, ctl->former);

        ctl->u_opt = ctl->solver->u_opt;

        return;
}}
