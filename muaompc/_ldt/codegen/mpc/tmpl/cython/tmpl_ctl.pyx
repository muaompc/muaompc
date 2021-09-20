"""{doc}
"""

from numpy import array
from cython.view cimport array as cvarray
ctypedef unsigned long long uint64_t

cimport {prefix}C{former} as C{former}
cimport {prefix}C{solver} as C{solver}
cimport {prefix}Cctl as Cctl

cdef class Prb(object):
    cdef Ccvp.{prefix}_cvp_prb *prb
{prb_class}
    def __cinit__(self, uint64_t prb):
        self.prb = <C{former}.{prefix}_{former}_prb *> prb
        self._update()

    def _update(self):
{prb_init}
{prb_property}
    def _set_property(self, uint64_t term):
        t = <C{former}.{prefix}_term *> term
        cdef cvarray cva = <double[:t.rows*t.cols]>t.data
        return (array(cva[:]), t.rows, t.cols)

cdef class Conf(object):
    cdef C{solver}.{prefix}_{solver}_conf *conf

    def __cinit__(self, uint64_t conf):
        self.conf = <C{solver}.{prefix}_{solver}_conf *> conf
{conf_property}

cdef class Parameters(object):
    cdef C{former}.{prefix}_{former} *former
{par_class}
    def __cinit__(self, uint64_t former):
        self.former = <C{former}.{prefix}_{former} *> former
{par_init}
{par_property}
    def _update(self):
        parnames = {parnames_list}
        for idx, parname in enumerate(parnames):
            self._set_parameter_pointer(idx, parname)

    def _set_parameter_property(self, idx):
        cdef C{former}.{prefix}_term *p = self.former.par[idx]
        cdef cvarray par = <double[:p.rows*p.cols]>p.data
        return array(par[:])

    def _set_parameter_pointer(self, idx, parname):
        cdef C{former}.{prefix}_term *p = self.former.par[idx]
        cdef cvarray par = <double[:p.rows*p.cols]>p.data
        ppar = getattr(self, parname)
        par[:] = ppar[:]

cdef class Ctl(object):
    cdef Cctl.{prefix}_ctl *ctl
    cdef Conf conf
    cdef Parameters parameters
    cdef Prb prb
    cdef cvarray u_opt

    def __cinit__(self):
        self.ctl = Cctl.{prefix}_ctl_allocate_ctl()
        if self.ctl is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self.ctl is not NULL:
            Cctl.{prefix}_ctl_free_ctl(self.ctl)

    property u_opt:
        def __get__(self):
            return array(self.u_opt)

    property conf:
        def __get__(self):
            return self.conf

    property parameters:
        def __get__(self):
            return self.parameters

    property prb:
        def __get__(self):
            return self.prb

    def __init__(self, fname):
        fname = fname.encode()  # python 3, char *
        Cctl.{prefix}_ctl_setup_ctl(self.ctl, fname)
        n = self.ctl.solver.optvar_seqlen
        self.u_opt = <double[:n]> self.ctl.u_opt
        cdef C{solver}.{prefix}_{solver}_conf *conf = self.ctl.solver.conf
        self.conf = Conf(<uint64_t> conf)
        cdef C{former}.{prefix}_{former} *former = self.ctl.former
        self.parameters = Parameters(<uint64_t> former)
        prb = <C{former}.{prefix}_{former}_prb *> former.prb
        self.prb = Prb(<uint64_t> prb)

    def solve_problem(self):
        self.parameters._update()
        Cctl.{prefix}_ctl_solve_problem(self.ctl)
        n = self.ctl.solver.optvar_seqlen
        self.u_opt = <double[:n]> self.ctl.u_opt
        self.prb._update()
