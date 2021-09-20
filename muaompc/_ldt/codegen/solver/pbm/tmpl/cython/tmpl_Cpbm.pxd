cimport {prefix}Ccvp as Ccvp

cdef extern from "../include/{prefix}pbm.h":

    cdef struct {prefix}_pbm_conf:
        int in_iter
        int warm_start

    cdef struct {prefix}_pbm:
        int *j_in
        double *z_ini
        double *z_opt
        int optvar_seqlen
        {prefix}_pbm_conf *conf

    void {prefix}_pbm_solve_problem({prefix}_pbm *pbm)

cdef extern from "../include/{prefix}pbmdynmem.h":

    {prefix}_pbm *{prefix}_pbm_allocate_solver()
    void {prefix}_pbm_setup_solver({prefix}_pbm *pbm, Ccvp.{prefix}_cvp_prb *prb, char *fname)
