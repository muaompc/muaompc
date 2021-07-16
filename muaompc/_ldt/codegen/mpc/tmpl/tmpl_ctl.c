#include "{prefix}ctl.h"

void {prefix}_ctl_form_problem(struct {prefix}_ctl *ctl)
{{
  {prefix}_{former}_form_problem(ctl->former);
  return;
}}

void {prefix}_ctl_solve_problem(struct {prefix}_ctl *ctl)
{{
  {prefix}_{former}_copy_parameters(ctl->former, ctl->parameters);
  {prefix}_{former}_form_problem(ctl->former);
  {prefix}_{solver}_solve_problem(ctl->solver);
  return;
}}

