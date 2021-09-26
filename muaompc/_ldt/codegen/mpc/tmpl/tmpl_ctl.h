#ifndef {PREFIX}CTL_H
#define {PREFIX}CTL_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}{solver}.h"
#include "{prefix}{former}.h"

struct {prefix}_ctl {{
struct {prefix}_{former}_parameters *parameters;  /**< Problem parameters. */
struct {prefix}_{former} *former;  /**< Problem former. */
struct {prefix}_{solver} *solver;  /**< Problem solver. */

real_t *u_opt;  /**< Optimal control input sequence. */
}};  /**< The main {PREFIX} structure. Contains algorithm, system and runtime data. */

extern void {prefix}_ctl_form_problem(struct {prefix}_ctl *ctl);
extern void {prefix}_ctl_solve_problem(struct {prefix}_ctl *ctl);

#ifdef __cplusplus
}}
#endif

#endif /* {PREFIX}CTL_H */
