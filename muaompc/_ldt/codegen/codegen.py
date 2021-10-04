"""Generate C code for the structure of a given problem."""

import shutil
import os
import errno
import json
import importlib

import numpy as np

from muaompc.version import version


def os_makedirs(path, exist_ok=False):
# replicate functionality of Python3 os.makedirs
# the parameter exist_ok is not supported in Python2 makedirs
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

class BaseCodeGenerator(object):

    """Generate C code for the structure of a given problem.

    This class is intended to be subclassed to generate C-code,
    and interfaces to the C-code.
    """

    def __init__(self):
        pass

    def _get_paths(self, name):
        path = dict()
        module = importlib.import_module('muaompc._ldt.'+name)
        path['src'] = module.__path__[0]  # the path to the C template files
        path['tmpl'] = os.path.join(path['src'], 'tmpl')
        path['root'] = os.path.join(os.path.abspath(self.destdir),
                                    self.prefix+'_'+self.prbname)
        path['dest'] = os.path.join(path['root'], 'src')
        path['data'] = os.path.join(path['root'], 'data')

        return path

    def _make_destination_dir_tree(self):
        os_makedirs(self.path['dest'], exist_ok=True)
        os_makedirs(self.path['data'], exist_ok=True)

    def _copy_static_src(self, src):
        for name in src:
            fname = os.path.join(self.path['tmpl'], name)
            shutil.copy(fname, self.path['dest'])
        try:
            shutil.copytree(os.path.join(self.path['tmpl'], 'include'),
                            os.path.join(self.path['dest'], 'include'))
        except OSError:
            pass

    def _get_tmpl(self, fname, srcdir='.'):
        fname = os.path.join(self.path['src'], 'tmpl', srcdir, 'tmpl_'+fname)
        with open(fname, 'r') as f:
            tmpl = f.read()
        return tmpl

    def _write_file(self, fstring, fname, destdir=None):
        path = self.path['dest']
        if destdir is not None:
            path = os.path.join(self.path['dest'], destdir)
        else:
            if fname.endswith('.h'):
                path = os.path.join(self.path['dest'], 'include')
            if fname.endswith(('.py', '.pxd', '.pyx', 'pxi')):
                path = os.path.join(self.path['dest'], 'cython')

        with open(os.path.join(path, fname), 'w') as f:
            f.write(fstring)

    def _replace_dict(self, dict_, prefix, fname, srcdir='.', destdir=None):
        tmpl = self._get_tmpl(fname, srcdir)
        fmt = tmpl.format(**dict_)
        pfname = prefix+fname
        self._write_file(fmt, pfname, destdir)

    def _replace_prefix(self, fname, srcdir='.', destdir=None):
        if self.prefix=='hhmpc':  # TODO better condition.
            pbm_sparse = '#define %s_PBM_PRB_SPARSE 1\n' % (self.prefix.upper())
            pbm_cond = '#define %s_PBM_PRB_COND 0\n' % (self.prefix.upper())
        else:
            pbm_sparse = '#define %s_PBM_PRB_SPARSE 0\n' % (self.prefix.upper())
            pbm_cond = '#define %s_PBM_PRB_COND 1\n' % (self.prefix.upper())
        tmp = dict(prefix=self.prefix,
                   PREFIX=self.prefix.upper(),
                   pbm_sparse=pbm_sparse,
                   pbm_cond=pbm_cond)
        self._replace_dict(tmp, self.prefix, fname, srcdir=srcdir, destdir=destdir)

    def _get_enum_terms_c_tmpl(self, terms, ntabs=1):
        fmt = ''
        for k, name in enumerate(terms):
            fmt += self._get_term_enum_c_tmpl(name, k, ntabs)
        return fmt

    def _get_term_enum_c_tmpl(self, term_name, k, ntabs):
        tab = '    '
        tabs = ntabs*tab
        enum_i = tabs+'%s_%s = %d,\n' % (self.prefix, term_name, k)
        return enum_i.upper()

    def _np2Carray(self, arr):
        arr = str(list(arr)).replace('[', '{')
        arr = arr.replace(']', '}')
        return arr


class BaseDataGenerator(object):
    """Generate json data for the structure of a given problem.

    This class is intended to be used as is, to generate C code,
    or be subclassed to generate interfaces to the C code.

    """

    def __init__(self):
        pass

    def _get_full_fname_prefix(self, dataname):
        datadir = os.path.join(self.path['data'], dataname)
        os_makedirs(datadir, exist_ok=True)

        fname = os.path.join(datadir, self.prefix+dataname)
        return fname

    def _write_json_data(self, fname, data):
        with open(fname, 'w') as f:
            json.dump(data, f)

    def _get_par(self, num):
        par_t = dict()
        for par in self.smb.par:
            par_t[par.name] = self._get_seq_term(par, num)
        return par_t

    def _get_optvar_len(self, num):
        opt = self.smb.opt
        return dict(seqlen=opt.eval_seqlen(num), veclen=eval(opt.veclen, num))

    def _get_constant(self, constant_):
        constant = dict()
        for name, mtx in constant_.items():
            constant[name] = self._get_constant_mtx_term(mtx)
        return constant

    def _get_seq_term(self, term, num):
        rows = term.eval_seqlen(num)

        cols = 1
        term_lst = np.zeros(rows).tolist()  # values provided online
        return dict(data=term_lst, rows=rows, cols=cols)

    def _get_constant_mtx_term(self, mtx):
        rows, cols = mtx.shape
        mtx_lst = mtx.astype('float').flatten().tolist()  # C array
        return dict(data=mtx_lst, rows=rows, cols=cols)

codegen_doc = "This file was automatically generated by muaompc " + version  