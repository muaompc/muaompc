cimport {prefix}Ccvp as Ccvp

cdef extern from "../include/{prefix}alm.h":

    cdef struct {prefix}_alm_conf:
        int in_iter
        int ex_iter
        int warm_start

    cdef struct {prefix}_alm:
        int *i_ex
        double *u_ini
        double *u_opt
        double *l_ini
        double *l_opt
        int optvar_seqlen
        int lagmul_seqlen
        {prefix}_alm_conf *conf

    void {prefix}_alm_solve_problem({prefix}_alm *alm)

cdef extern from "../include/{prefix}almdynmem.h":

    {prefix}_alm *{prefix}_alm_allocate_solver()
    void {prefix}_alm_setup_solver({prefix}_alm *alm, Ccvp.{prefix}_cvp_prb *prb, char *fname)
