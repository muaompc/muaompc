import json

from cython.view cimport array as cvarray

cimport {prefix}Ccvp as Ccvp
cimport {prefix}Calm as Calm

cdef class Solver:
    cdef Calm.{prefix}_alm *alm
    cdef Ccvp.{prefix}_cvp *cvp
    cdef cvarray u_opt
    cdef cvarray u_ini
    cdef cvarray l_opt
    cdef cvarray l_ini

    def __cinit__(self):
        self.alm = Calm.{prefix}_alm_allocate_solver()
        self.cvp = Ccvp.{prefix}_cvp_allocate_former()

    def __dealloc__(self):
        # TODO: free memory allocated by malloc
        pass

    property u_opt:
        def __get__(self):
          return self.u_opt

    property u_ini:
        def __get__(self):
          return self.u_ini
        def __set__(self, u_ini):
          cdef double [:]u_ini_cv = u_ini
          self.u_ini[:] = u_ini_cv

    property l_opt:
        def __get__(self):
          return self.l_opt

    property l_ini:
        def __get__(self):
          return self.l_ini
        def __set__(self, l_ini):
          cdef double [:]l_ini_cv = l_ini
          self.l_ini[:] = l_ini_cv

    cpdef setup_solver(self, data, fname='data.json'):
        cdef int optvar_seqlen
        cdef int lagmul_seqlen
        # TODO: check data has consistent sizes
        with open(fname, 'w') as f:
            json.dump(data, f)

        fname = fname.encode()  # python 3, char *
        Ccvp.{prefix}_cvp_setup_former(self.cvp, fname)
        cdef Ccvp.{prefix}_cvp_prb *prb
        prb = self.cvp.prb
        Calm.{prefix}_alm_setup_solver(self.alm, prb, fname)
        optvar_seqlen = self.alm.optvar_seqlen
        lagmul_seqlen = self.alm.lagmul_seqlen
        self.u_opt = <double[:optvar_seqlen]> self.alm.u_opt
        self.u_ini = <double[:optvar_seqlen]> self.alm.u_ini
        self.l_opt = <double[:lagmul_seqlen]> self.alm.l_opt
        self.l_ini = <double[:lagmul_seqlen]> self.alm.l_ini

    cpdef solve_problem(self, pardata):
        cdef double *u
        cdef Ccvp.{prefix}_cvp *ccvp = self.cvp
        Ccvp._py2c_pardata(ccvp, pardata)
        Ccvp.{prefix}_cvp_form_problem(self.cvp)
        Calm.{prefix}_alm_solve_problem(self.alm)
        return

    cpdef configure(self, int ex_iter, int in_iter, int warm_start):
        self.alm.conf.in_iter = in_iter
        self.alm.conf.ex_iter = ex_iter
        self.alm.conf.warm_start = warm_start

