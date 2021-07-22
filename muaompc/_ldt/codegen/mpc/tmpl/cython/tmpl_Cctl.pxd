cimport {prefix}C{former} as C{former}
cimport {prefix}C{solver} as C{solver}

cdef extern from "../include/{prefix}ctl.h":

    cdef struct {prefix}_ctl:
        C{former}.{prefix}_{former} *former
        C{solver}.{prefix}_{solver} *solver
        double *u_opt

    void {prefix}_ctl_form_problem({prefix}_ctl *ctl)
    void {prefix}_ctl_solve_problem({prefix}_ctl *ctl)

cdef extern from "../include/{prefix}ctldynmem.h":

    {prefix}_ctl *{prefix}_ctl_allocate_ctl()
    int {prefix}_ctl_setup_ctl({prefix}_ctl *ctl, char *fname)
    {prefix}_ctl *{prefix}_ctl_free_ctl({prefix}_ctl *ctl)
