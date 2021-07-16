"""Generate C code and data for the fast gradient method."""

import json

import numpy as np

from muaompc._ldt.codegen.codegen import BaseCodeGenerator as BCG
from muaompc._ldt.codegen.former.cvp.codegen import CVPDataGenerator as CVPDG
from muaompc._ldt.codegen.former.cvp.codegen import CDataGenerator as CVPCDG


class CCodeGenerator(BCG, object):
    """Generate C code for the fast gradient method."""

    def __init__(self, ccg):
        """Assign the argument to self, and set the C code paths.

        This depends on common C-code already being generated.

        :param ccg: constains the prefix and the destination directory of the
        generated common C-code.
        :type ccg: muaompc._ldt.codegen.common.codegen.CCodeGenerator
        """
        self.prefix = ccg.prefix
        self.destdir = ccg.destdir
        self.prbname = ccg.prbname
        self.path = self._get_paths('codegen.solver.fgm')

    def generate_code(self):
        self._make_destination_dir_tree()
        self._generate_c_code()

    def _generate_c_code(self):
        self._generate_c_header('fgm')
        self._generate_c_body('fgm')
        self._generate_getdata_h()
        self._generate_getdata_c()

    def _generate_c_header(self, fname):
        fname += '.h'
        self._replace_prefix(fname)

    def _generate_c_body(self, fname):
        fname += '.c'
        self._replace_prefix(fname)

    def _generate_getdata_h(self):
        fname = 'fgmdynmem.h'
        self._replace_prefix(fname)

    def _generate_getdata_c(self):
        fname = 'fgmdynmem.c'
        self._replace_prefix(fname)


class FGMCVPDataGenerator(CVPDG, object):
    """
    Generate json data for the fast gradient method using a
    condensed vectorial parameters problem formulation.
    """

    def __init__(self, cvpcg):
        """Assign the arguments to self.

        :param cvp: an instance of condensed vectorial parameters formulation
        :type cvp: muaompc._ldt.codegen.former.cvp.cvp.CondensedVectorialParameters

        """
        CVPDG.__init__(self, cvpcg)

    def _get_data(self, num):
        data = CVPDG._get_data(self, num)
        H = np.array(data['pmetric']['H']['fac0']['data'])
        lenH = np.sqrt(H.shape[0])
        if lenH.is_integer() is False:
            raise TypeError
        nu, L = self._compute_extra_step_constant(np.reshape(H, (lenH, lenH)))
        data = self._scale_cost(data, H, L)
        data['nu'] = nu
        return data

    def _compute_extra_step_constant(self, H):
        from math import sqrt
        maxl = max(np.linalg.eigvalsh(H))
        minl = min(np.linalg.eigvalsh(H))
        nu = (sqrt(maxl) - sqrt(minl)) / (sqrt(maxl) + sqrt(minl))
        return nu, maxl

    def _scale_cost(self, data, H, L):
        data['pmetric']['H']['fac0']['data'] = (H / L).tolist()
        for k, fac in enumerate(data['pmetric']['g']['fac']):
            data['pmetric']['g']['fac'][k]['data'] = [di / L for di in fac['data']]

        data['pmetric']['g']['fac0']['data'] = [
            di / L for di in data['pmetric']['g']['fac0']['data']]
        return data


class CDataGenerator(CCodeGenerator, FGMCVPDataGenerator, object):

    def __init__(self, ccg, cvpcg):
        """Assign the arguments to self.

        :param ccg: an instance of C code generator of for the FGM
        :type ccg: muaompc._ldt.codegen.solver.fgm.codegen.CCodeGenerator

        :param cvpcg: an instance of C code generator of a
        condensed vectorial parameters formulation
        :type cvpcg: muaompc._ldt.codegen.former.cvp.codegen.CCodeGenerator

        """
        self.cvp = cvpcg.cvp
        self.smb = cvpcg.cvp.smb
        self.prefix = cvpcg.prefix
        self.prbname = ccg.prbname
        self.destdir = ccg.destdir
        self.path = self._get_paths('codegen.solver.fgm')
        self.cvpcg = cvpcg
        self.cvpcdg = CVPCDG(cvpcg)
        self.dict_ = None  # Assignment made during data generation

    def generate_data(self, num, dataname):
        self.dict_ = dict(prefix=self.prefix,
                          PREFIX=self.prefix.upper(),
                          dataname=dataname,
                          DATANAME=dataname.upper())
        self.cvpcdg.dict_ = self.dict_
        fgmcvpdg = FGMCVPDataGenerator(self.cvpcg)
        data = fgmcvpdg._get_data(num)
        self.cvpcdg._write_static_data(data, dataname)
        self._write_static_data(data, dataname)

    def _write_static_data(self, data, dataname):
        self._generate_h_struct(data, dataname)
        self._generate_c_struct(data, dataname)

    def _generate_h_struct(self, data, dataname):
        fname = 'fgmdata'
        fname += '.h'
        prefix = self._get_full_fname_prefix(dataname)
        tmp = self.dict_
        self._replace_dict(tmp, prefix, fname)

    def _generate_c_struct(self, data, dataname):
        fname = 'fgmdata'
        fname += '.c'
        prefix = self._get_full_fname_prefix(dataname)
        zeros_optvar_seqlen = self._np2Carray(np.zeros(data['optvar']['seqlen']))
        tmp = dict(nu=data['nu'],
                   zeros_optvar_seqlen=zeros_optvar_seqlen,
                   optvar_veclen=data['optvar']['veclen'])
        tmp.update(self.dict_)
        self._replace_dict(tmp, prefix, fname)
