"""Generate C code for an MPC problem."""

from muaompc._ldt.codegen.codegen import BaseCodeGenerator as BCG
from muaompc._ldt.codegen.codegen import BaseDataGenerator as BDG


class CCodeGenerator(BCG, object):

    def __init__(self, ccg, solver, former):
        """Generate C code for an MPC problem.

        :param str prefix: the text to be prepended to the names of the
        generated files and to the internal code global identifiers,
        e.g. names of functions, structures,
        defines, enums, etc. (like a namespace)
        :param str solver: the solver to be used to solve the MPC problem
        :param str former: the former to be used to form the MPC problem
        """
        self.prefix = ccg.prefix
        self.destdir = ccg.destdir
        self.prbname = ccg.prbname
        self.solver = solver
        self.former = former
        self.path = self._get_paths('codegen.mpc')

    def generate_code(self):
        self._make_destination_dir_tree()
        self._copy_static_src(['mpc.py', '__init__.py'])
        self._generate_setup_py()
        self._generate_c_code()

    def _generate_setup_py(self):
        self._replace_psf('setup.py', destdir='..')

    def _generate_c_code(self):
        self._replace_psf('ctl.h')
        self._replace_psf('ctldynmem.h')
        self._replace_psf('ctl.c')
        self._replace_psf('ctldynmem.c')

    def _replace_psf(self, fname, srcdir='.', destdir=None):
# replace prefix, solver, and former
        tmpl = self._get_tmpl(fname, srcdir)
        fmt = tmpl.format(prefix=self.prefix,
                          PREFIX=self.prefix.upper(),
                          solver=self.solver, former=self.former)
        pfname = self.prefix+fname
        self._write_file(fmt, pfname, destdir=destdir)


class CDataGenerator(BCG, BDG, object):

    def __init__(self, mpc, solver, former):
        self.prefix = mpc.prefix
        self.destdir = mpc.destdir
        self.prbname = mpc.prbname
        self.solver = solver
        self.former = former
        self.path = self._get_paths('codegen.mpc')
        self.dict_ = None  # Assignment made during data generation

    def generate_data(self, data, dataname):
        self.dict_ = dict(prefix=self.prefix,
                          PREFIX=self.prefix.upper(),
                          dataname=dataname,
                          DATANAME=dataname.upper(),
                          former=self.former,
                          solver=self.solver)
        #FIXME data is not being used
        self._write_static_data(data, dataname)

    def _write_static_data(self, data, dataname):
        self._generate_h_struct(data, dataname)
        self._generate_c_struct(data, dataname)
        self._generate_main(dataname)
        self._generate_Makefile(dataname)

    def _generate_h_struct(self, data, dataname):
        fname = 'ctldata'
        fname += '.h'
        prefix = self._get_full_fname_prefix(dataname)
        tmp = self.dict_
        self._replace_dict(tmp, prefix, fname)

    def _generate_c_struct(self, data, dataname):
        fname = 'ctldata'
        fname += '.c'
        prefix = self._get_full_fname_prefix(dataname)
        tmp = self.dict_
        self._replace_dict(tmp, prefix, fname)

    def _generate_main(self, dataname):
        fname = 'main'
        fname += '.c'
        prefix = self._get_full_fname_prefix(dataname)
        tmp = self.dict_
        self._replace_dict(tmp, prefix, fname)

    def _generate_Makefile(self, dataname):
        fname = 'Makefile.mk'
        prefix = self._get_full_fname_prefix(dataname)
        tmp = self.dict_
        self._replace_dict(tmp, prefix, fname)
