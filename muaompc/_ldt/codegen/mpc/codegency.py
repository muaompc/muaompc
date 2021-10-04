"""Generate Cython interface for C code for the controller.

This is the main interface exposed to the users.
It does not aim to expose all details of the solver nor former algorithms,
only the main configuration variables related to the MPC controller.
"""

from muaompc._ldt import codegen
import shutil
import os

from muaompc._ldt.codegen.codegen import codegen_doc
from muaompc._ldt.codegen.mpc.codegen import CCodeGenerator as MCG


class CythonCodeGenerator(MCG, object):

    def __init__(self, mcg, cvpcg):
        self.prefix = mcg.prefix
        self.path = mcg.path
        self.cvp = cvpcg.cvp
        self.former = mcg.former
        self.solver = mcg.solver
        self.base = dict(prefix=self.prefix,
                   PREFIX=self.prefix.upper(),
                   former=self.former,
                   solver=self.solver,
                   doc = codegen_doc
        )

    def generate_code(self):
        if self.solver == 'fgm':
            dod = _generate_dict_fgm(self.base)
        elif self.solver == 'alm':
            dod = _generate_dict_alm(self.base)
        elif self.solver == 'pbm':
            dod = _generate_dict_pbm(self.base)
        else:
            raise TypeError("Solver %s not recognized." % self.solver)
        dod['pyx'].update(_generate_dict_parameters(self.cvp.par))
        prbterms = ['H', 'g', 'u_lb', 'u_ub']
        if not self.cvp.is_inputbox_iec:  # mixed constraints
            prbterms += ['V', 'v_lb', 'v_ub']

        dod['pyx'].update(_generate_dict_prb(prbterms))
        self._generate_cython_code(dod)

    def _generate_cython_code(self, dod):
        # dod: a dict of dicts
        self._replace_dict(dod['setup'], self.prefix, 'ctlsetup.py', srcdir='cython')
        self._replace_dict(dod['pyx'], self.prefix, 'ctl.pyx', srcdir='cython')
        self._replace_dict(dod['pxd'], self.prefix, 'Cctl.pxd', srcdir='cython')

def _generate_dict_parameters(parnames):
    return dict(par_class=_get_class(parnames),
                par_property=_get_par_property(parnames),
                par_init=_get_par_init(parnames),
                parnames_list=str(list(parnames.keys()))
                )

def _get_class(names):
    fmt = ''
    tab = '    '
    for name in names:
        fmt += tab + 'cdef object %s\n' % (name)
    return fmt

_par_property = """
    property {parname}:
        def __get__(self):
            return self.{parname}
"""

def _get_par_property(parnames):
    fmt = ""
    for parname in parnames:
        fmt += _par_property.format(**dict(parname=parname))
    return fmt

def _get_par_init(parnames):
    fmt = ""
    tab = '    '
    tabs = 2*tab
    for idx, parname in enumerate(parnames):
        fmt += tabs + "self.%s = self._set_parameter_property(%d)\n" % (parname, idx)
    return fmt

def _generate_dict_prb(names):
    return dict(prb_class=_get_class(names),
                prb_property=_get_prb_property(names),
                prb_init=_get_prb_init(names),
                )

_prb_property = """
    property {name}:
        def __get__(self):
            data = self.{name}[0]
            rows = self.{name}[1]
            cols = self.{name}[2]
            return data.reshape((rows, cols))
"""

def _get_prb_property(names):
    fmt = ""
    for name in names:
        fmt += _prb_property.format(**dict(name=name))
    return fmt

def _get_prb_init(names):
    fmt = ""
    tab = '    '
    tabs = 2*tab
    for name in names:
        fmt += tabs + "self.%s = self._set_property(<uint64_t> self.prb.%s)\n" % (name, name)
    return fmt

def _generate_dict_fgm(base):
    dod = dict()
    dod['setup'] = dict(base,
                solver_dep="",
    )
    dod['pyx'] = _generate_dict_pyx_fgm(base)
    dod['pxd'] = dict(base)

    return dod

def _generate_dict_alm(base):
    # setup of ALM is identical to FGM, except that now FGM files are a dependency
    dod = _generate_dict_fgm(base)
    solver_dep = '"src/{prefix}{solver}dynmem.c","src/{prefix}{solver}.c",'
    dod['setup']['solver_dep'] = solver_dep.format(**dict(prefix=base['prefix'], 
                                                    solver='fgm'))
    dod['pyx'] = _generate_dict_pyx_alm(base)

    return dod

def _generate_dict_pbm(base):
    dod = _generate_dict_fgm(base)
    solver_dep = '"src/{prefix}pbmsolve.c",'
    dod['setup']['solver_dep'] = solver_dep.format(**dict(prefix=base['prefix']))
    return dod

def _generate_dict_pyx_fgm(base):
    items = ['warm_start', 'in_iter']
    return dict(base,
            conf_property=_gen_conf_property(items)
            )

def _generate_dict_pyx_alm(base):
    items = ['warm_start', 'in_iter', 'ex_iter']
    return dict(base,
            conf_property=_gen_conf_property(items)
            )

_property = """
    property {item}:
        def __get__(self):
            return self.conf.{item}

        def __set__(self, val):
            self.conf.{item} = val
"""

def _gen_conf_property(items):
    cp = ""
    for item in items:
        cp += _property.format(**dict(item=item))
    return cp
