import json

from cython.view cimport array as cvarray

cimport {prefix}Ccvp as Ccvp
cimport {prefix}Cpbm as Cpbm

cdef class Solver:
    cdef Cpbm.{prefix}_pbm *pbm
    cdef Ccvp.{prefix}_cvp *cvp
    cdef cvarray z_opt
    cdef cvarray z_ini

    def __cinit__(self):
        self.pbm = Cpbm.{prefix}_pbm_allocate_solver()
        self.cvp = Ccvp.{prefix}_cvp_allocate_former()

    def __dealloc__(self):
        # TODO: free memory allocated by malloc
        pass

    property z_opt:
        def __get__(self):
          return self.z_opt

    property z_ini:
        def __get__(self):
          return self.z_ini
        def __set__(self, z_ini):
          cdef double [:]z_ini_cv = z_ini
          self.z_ini[:] = z_ini_cv

    cpdef setup_solver(self, data, fname='data.json'):
        cdef int optvar_seqlen
        # TODO: check data has consistent sizes
        with open(fname, 'w') as f:
            json.dump(data, f)

        fname = fname.encode()  # python 3, char *
        Ccvp.{prefix}_cvp_setup_former(self.cvp, fname)
        cdef Ccvp.{prefix}_cvp_prb *prb
        prb = self.cvp.prb
        Cpbm.{prefix}_pbm_setup_solver(self.pbm, prb, fname)
        optvar_seqlen = self.pbm.optvar_seqlen
        self.z_opt = <double[:optvar_seqlen]> self.pbm.z_opt
        self.z_ini = <double[:optvar_seqlen]> self.pbm.z_ini

    cpdef solve_problem(self, pardata):
        cdef double *u
        cdef Ccvp.{prefix}_cvp *ccvp = self.cvp
        Ccvp._py2c_pardata(ccvp, pardata)
        Ccvp.{prefix}_cvp_form_problem(self.cvp)
        Cpbm.{prefix}_pbm_solve_problem(self.pbm)
        return

    cpdef configure(self, int in_iter, int warm_start):
        self.pbm.conf.in_iter = in_iter
        self.pbm.conf.warm_start = warm_start

