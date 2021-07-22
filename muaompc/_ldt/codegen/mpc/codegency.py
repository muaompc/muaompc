"""Generate Cython interface for C code for the controller."""
import shutil
import os

from muaompc._ldt.codegen.codegen import BaseCodeGenerator as BCG
from muaompc._ldt.codegen.mpc.codegen import CCodeGenerator as MCG


class CythonCodeGenerator(MCG, object):

    def __init__(self, mcg):
        self.prefix = mcg.prefix
        self.path = mcg.path
        self.former = mcg.former
        self.solver = mcg.solver

    def generate_code(self):
        self._generate_cython_code()

    def _generate_cython_code(self):
        d = dict(prefix=self.prefix,
                   PREFIX=self.prefix.upper(),
                   former=self.former,
                   solver=self.solver,
        )
        self._replace_dict(d, self.prefix, 'ctlsetup.py', srcdir='cython')
        self._replace_dict(d, self.prefix, 'ctl.pyx', srcdir='cython')
        self._replace_dict(d, self.prefix, 'Cctl.pxd', srcdir='cython')
