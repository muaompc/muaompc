"""Generate Cython interface for C code for the controller."""
from muaompc._ldt import codegen
import shutil
import os

from muaompc._ldt.codegen.codegen import codegen_doc
from muaompc._ldt.codegen.mpc.codegen import CCodeGenerator as MCG


class CythonCodeGenerator(MCG, object):

    def __init__(self, mcg):
        self.prefix = mcg.prefix
        self.path = mcg.path
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
        else:
            raise TypeError
        self._generate_cython_code(dod)

    def _generate_cython_code(self, dod):
        # dod: a dict of dicts
        self._replace_dict(dod['setup'], self.prefix, 'ctlsetup.py', srcdir='cython')
        self._replace_dict(dod['pyx'], self.prefix, 'ctl.pyx', srcdir='cython')
        self._replace_dict(dod['pxd'], self.prefix, 'Cctl.pxd', srcdir='cython')


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

def _generate_dict_pyx_fgm(base):
    items = ['warm_start', 'in_iter']
    return dict(base,
            conf_class=_gen_conf_class(items),
            conf_cinit=_gen_conf_cinit(items),
            conf_property=_gen_conf_property(items)
            )

def _generate_dict_pyx_alm(base):
    items = ['warm_start', 'in_iter', 'ex_iter']
    return dict(base,
            conf_class=_gen_conf_class(items),
            conf_cinit=_gen_conf_cinit(items),
            conf_property=_gen_conf_property(items)
            )

def _gen_conf_class(items):
    cc = ""
    for item in items:
        cc += "    cdef int %s\n" % item
    return cc

def _gen_conf_cinit(items):
    ci = ""
    for item in items:
        ci += "        self.%s = self.conf.%s\n" % (item, item)
    return ci

_property = """
    property {item}:
        def __get__(self):
            return self.{item}

        def __set__(self, val):
            self.conf.{item} = val
            self.{item} = val
"""

def _gen_conf_property(items):
    cp = ""
    for item in items:
        cp += _property.format(**dict(item=item))
    return cp