import json

cimport stocmpcCcvp as Ccvp
from cython.view cimport array

cdef class Term:
    cdef int rows
    cdef int cols
    cdef array data

    def __cinit__(self, int rows, int cols, array data):
        self.data = data
        self.rows = rows
        self.cols = cols

    property rows:
        def __get__(self):
            return self.rows

    property cols:
        def __get__(self):
            return self.cols

    property data:
        def __get__(self):
            return self.data

cdef class Prb:
    cdef object H
    def __cinit__(self, H):
        self.H = H
        #elf.g = g
        #elf.u_lb = u_lb
        #elf.u_ub = u_ub

    property H:
        def __get__(self):
            return self.H
#        def __set__(self, H):
#            self.H = H
#    property g:
#        def __get__(self):
#            return self.g
#    property u_lb:
#        def __get__(self):
#            return self.u_lb
#    property u_ub:
#        def __get__(self):
#            return self.u_ub


cdef class Former:
    cdef Ccvp.stocmpc_cvp *cvp
    cdef object prb
    # automatically generated, just as in Ccvp.pxd struct prb

    def __cinit__(self):
        self.cvp = Ccvp.stocmpc_cvp_allocate_former()

    def __dealloc__(self):
        # TODO: free memory allocated by malloc
        pass

    def initialize(self, data, fname='data.json'):
        # TODO: check data has consistent sizes
        with open(fname, 'w') as f:
            json.dump(data, f)

        fname = fname.encode()  # python 3, char *
        Ccvp.stocmpc_cvp_setup_former(self.cvp, fname)
        cdef Ccvp.stocmpc_cvp *ccvp = self.cvp
        cdef Ccvp.stocmpc_cvp_prb *cprb = ccvp.prb
        cdef Ccvp.stocmpc_term *H = cprb.H
        cdef Ccvp.stocmpc_term *g = cprb.g
        cdef Ccvp.stocmpc_term *u_lb = cprb.u_lb
        cdef Ccvp.stocmpc_term *u_ub = cprb.u_ub
        self.prb = Prb(
            Term(H.rows, H.cols, <double[:H.rows*H.cols]> H.data))
        #
        #   Term(g.rows, g.cols, <double[:g.rows*g.cols]> g.data),
        #   Term(u_lb.rows, u_lb.cols, <double[:u_lb.rows*u_lb.cols]> u_lb.data),
        #   Term(u_ub.rows, u_ub.cols, <double[:u_ub.rows*u_ub.cols]> u_ub.data)
        #   )

    def form_problem(self, pardata):
        cdef Ccvp.stocmpc_cvp *ccvp = self.cvp
        cdef Ccvp.stocmpc_cvp_prb *cprb = ccvp.prb
        Ccvp._py2c_pardata(ccvp, pardata)
        Ccvp.stocmpc_cvp_form_problem(ccvp)
        prb = Ccvp._c2py_qp(cprb)
        return prb

    property prb:
        def __get__(self):
            return self.prb
#       def __set__(self, prb):
#           self.prb = prb
