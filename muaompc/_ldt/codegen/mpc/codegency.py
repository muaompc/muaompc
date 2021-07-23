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

    def generate_code(self):
        if self.solver == 'fgm':
            d = self._generate_dict_fgm()
        elif self.solver == 'alm':
            d = self._generate_dict_alm()
        else:
            raise TypeError

        self._generate_cython_code(d)

    def _generate_dict_fgm(self):
        return dict(prefix=self.prefix,
                   PREFIX=self.prefix.upper(),
                   former=self.former,
                   solver=self.solver,
                   solver_dep="",
                   doc = codegen_doc
        )

    def _generate_dict_alm(self):
        # code gen of ALM is identical to FGM, except that now FGM files are a depency
        d = self._generate_dict_fgm()
        solver_dep = '"src/{prefix}{solver}dynmem.c","src/{prefix}{solver}.c",'
        d['solver_dep'] = solver_dep.format(**dict(prefix=self.prefix, solver='fgm'))
        return d

    def _generate_cython_code(self, d):
        self._replace_dict(d, self.prefix, 'ctlsetup.py', srcdir='cython')
        self._replace_dict(d, self.prefix, 'ctl.pyx', srcdir='cython')
        self._replace_dict(d, self.prefix, 'Cctl.pxd', srcdir='cython')
