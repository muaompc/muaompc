cimport {prefix}Ccvp as Ccvp

cdef extern from "../include/{prefix}fgm.h":

    cdef struct {prefix}_fgm_conf:
        int in_iter
        int warm_start

    cdef struct {prefix}_fgm:
        int *j_in
        double *u_ini
        double *u_opt
        int optvar_seqlen
        {prefix}_fgm_conf *conf

    void {prefix}_fgm_solve_problem({prefix}_fgm *fgm)

cdef extern from "../include/{prefix}fgmdynmem.h":

    {prefix}_fgm *{prefix}_fgm_allocate_solver()
    void {prefix}_fgm_setup_solver({prefix}_fgm *fgm, Ccvp.{prefix}_cvp_prb *prb, char *fname)
