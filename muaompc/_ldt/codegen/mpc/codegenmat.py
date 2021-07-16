"""Generate Matlab interface for C code for an MPC controller.
"""

import os

from muaompc._ldt.codegen.mpc.codegen import CCodeGenerator as CCG


class MatlabCodeGenerator(CCG, object):

    """Generate Matlab code to interface with C code for an MPC controller

    This class overrides the following methods:
    generate_code: generates Matlab specific code on its own directory
    """

    def __init__(self, ccg):
        """Assign to self prefix, solver, former and path attributes from
        the given C code.

        This depends on C-code already being generated by C-code
        generation class. This creates the interface to that specific code.

        :param ccg: an instance of the C code generator class for which
        a Matlab interface is to be created
        :type ccg: muaompc._ldt.codegen.mpc.codegen.CCodeGenerator
        """
        self.prefix = ccg.prefix
        self.solver = ccg.solver
        self.former = ccg.former
        self.path = ccg.path
        self.classpath = os.path.join(self.path['dest'],
                'matlab', '@'+self.prefix+'ctl')

    def generate_code(self):
        self._make_destination_dir_tree()
        self._generate_matlab_code()

    def generate_matlab_make(self):
        # FIXME: Quick and dirty solution
        # this must be called last, it relies on all sources being generated
        fname = 'make.m'
        t = self._get_tmpl(fname, srcdir='matlab')
        c_src = self._get_c_src()
        s = t.format(prefix=self.prefix,
                c_src=c_src)
        pfname = self.prefix+fname
        self._write_file(s, pfname, destdir='matlab')

    def _get_c_src(self):
        src = []
        u = os.listdir(self.path['dest'])
        for file in u:
            if file.endswith('.c'):
                src.append(file)
        u_src = ''.join([' ../' + name for name in src])
        src = []
        m = os.listdir(os.path.join(self.path['dest'], 'matlab'))
        for file in m:
            if file.endswith('matlab.c'):
                src.append(file)
        m_src = ''.join([' ' + name for name in src])
        return u_src + m_src

    def _make_destination_dir_tree(self):
        try:
            os.mkdir(self.classpath)
        except OSError:
            pass

    def _generate_matlab_code(self):
        include_path = os.path.join('matlab', 'include')
        self._replace_prefix('ctlmatlab.h', srcdir='matlab', destdir=include_path)
        self._replace_psf('_matlab_create_ctl.c', srcdir='matlab', destdir='matlab')
        self._replace_psf('ctlmatlab.c', srcdir='matlab', destdir='matlab')
        self._replace_prefix('_matlab_form_problem.c', srcdir='matlab', destdir='matlab')
        self._replace_prefix('_matlab_solve_problem.c', srcdir='matlab', destdir='matlab')
        self._replace_prefix('ctl.m', srcdir='matlab', destdir=self.classpath)