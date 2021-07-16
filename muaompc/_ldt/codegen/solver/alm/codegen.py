"""Generate C code and data for the augmented Lagrangian method."""

import numpy as np

from muaompc._ldt.codegen.codegen import BaseCodeGenerator as BCG
from muaompc._ldt.codegen.former.cvp.codegen import CVPDataGenerator as CVPDG
from muaompc._ldt.codegen.former.cvp.codegen import CDataGenerator as CVPCDG
from muaompc._ldt.codegen.solver.fgm.codegen import CCodeGenerator as FGMCG
from muaompc._ldt.codegen.solver.fgm.codegen import (FGMCVPDataGenerator
                                                     as FGMCVPDG)
from muaompc._ldt.codegen.solver.fgm.codegen import CDataGenerator as FGMCVPCDG


class CCodeGenerator(BCG, object):
    """Generate C code for the augmented Lagrangian method."""

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
        self.path = self._get_paths('codegen.solver.alm')
        self.fgm = FGMCG(ccg)

    def generate_code(self):
        self.fgm.generate_code()
        self._make_destination_dir_tree()
        self._generate_c_code()

    def _generate_c_code(self):
        self._generate_c_header('alm')
        self._generate_c_body('alm')
        self._generate_getdata_h()
        self._generate_getdata_c()

    def _generate_c_header(self, fname):
        fname += '.h'
        self._replace_prefix(fname)

    def _generate_c_body(self, fname):
        fname += '.c'
        self._replace_prefix(fname)

    def _generate_getdata_h(self):
        fname = 'almdynmem.h'
        self._replace_prefix(fname)

    def _generate_getdata_c(self):
        fname = 'almdynmem.c'
        self._replace_prefix(fname)


class ALMCVPDataGenerator(CVPDG, object):
    """
    Generate json data for the augmented Lagrangian method using a
    condensed vectorial parameters problem formulation.
    """

    def __init__(self, cvpcg):
        """Assign the arguments to self.

        :param cvp: an instance of condensed vectorial parameters formulation
        :type cvp:
        muaompc._ldt.codegen.former.cvp.cvp.CondensedVectorialParameters

        """
        CVPDG.__init__(self, cvpcg)
        self.fgm = FGMCVPDG(cvpcg)

    def _get_data(self, num):
        data = CVPDG._get_data(self, num)
        H, V = self._get_matrices_HV(num)
        try:
            num['mu']
        except KeyError:
            data['mu'] = self._find_penalty_parameter(H, V)
        else:
            data['mu'] = num['mu']
        nu, L = self._compute_extra_step_constant(H, data['mu'], V)
        data_scaled = self.fgm._scale_cost(data, H.flatten(), L)
        data_scaled['nu'] = nu
        data_scaled['Linv'] = 1. / L
        return data_scaled

    def _get_matrices_HV(self, num): 
        data = CVPDG._get_data(self, num)
        Ht = data['pmetric']['H']['fac0']
        H = np.reshape(np.array(Ht['data']), (Ht['rows'], Ht['cols']))
        Vt = data['pmetric']['V']['fac0']
        V = np.reshape(np.array(Vt['data']), (Vt['rows'], Vt['cols']))
        return H, V

    def _compute_extra_step_constant(self, H, mu, V):
        from math import sqrt
        L = max(np.linalg.eigvalsh(H + mu * np.dot(V.T, V)))
        phi = min(np.linalg.eigvalsh(H))
        nu = ((sqrt(L) - sqrt(phi)) / (sqrt(L) + sqrt(phi)))
        return (nu, L)

    def _find_penalty_parameter(self, H, V):
        """Returns a penalty parameter that makes the condition number of the
        internal problem (FGM) not much larger than the Hessian's
        condition number.
        """
        msg = '_find_penalty_parameter '

        lower_limit, upper_limit, max_iter = (2, 10, 100)
        mu = 1.
        cnH = self._compute_alg_condnum(H, 0., V)
        cn_fact = self._compute_alg_condnum(H, mu, V) / cnH

        for k in range(max_iter):
            if cn_fact > upper_limit:
                mu /= 2
            elif cn_fact < lower_limit:
                mu *= 2
            else:
                break
            cn_fact = self._compute_alg_condnum(H, mu, V) / cnH
            if k == max_iter - 1:
                msg += 'WARNING: reached maximum number of iterations.'

        return mu

    def _compute_alg_condnum(self, H, mu, V):
        dummy, L = self._compute_extra_step_constant(H, mu, V)
        return L / (min(np.linalg.eigvals(H)))


class CDataGenerator(CCodeGenerator, ALMCVPDataGenerator, object):

    def __init__(self, ccg, cvpcg):
        """Assign the arguments to self.

        :param ccg: an instance of C code generator for the ALM
        :type ccg: muaompc._ldt.codegen.solver.alm.codegen.CCodeGenerator

        :param cvpcg: an instance of C code generator of a
        condensed vectorial parameters formulation
        :type cvpcg: muaompc._ldt.codegen.former.cvp.codegen.CCodeGenerator

        """
        self.cvp = cvpcg.cvp
        self.smb = cvpcg.cvp.smb
        self.prbname = ccg.prbname
        self.prefix = cvpcg.prefix
        self.destdir = ccg.destdir
        self.path = self._get_paths('codegen.solver.alm')
        self.cvpcg = cvpcg
        self.cvpcdg = CVPCDG(cvpcg)
        self.fgmcvpcdg = FGMCVPCDG(ccg, cvpcg)
        self.dict_ = None  # Assignment made during data generation

    def generate_data(self, num, dataname):
        self.dict_ = dict(prefix=self.prefix,
                          PREFIX=self.prefix.upper(),
                          dataname=dataname,
                          DATANAME=dataname.upper())
        self.cvpcdg.dict_ = self.dict_
        self.fgmcvpcdg.dict_ = self.dict_

        almcvpdg = ALMCVPDataGenerator(self.cvpcg)
        data = almcvpdg._get_data(num)
        self.cvpcdg._write_static_data(data, dataname)
        self.fgmcvpcdg._write_static_data(data, dataname)
        self._write_static_data(data, dataname)

    def _write_static_data(self, data, dataname):
        self._generate_h_struct(data, dataname)
        self._generate_c_struct(data, dataname)

    def _generate_h_struct(self, data, dataname):
        fname = 'almdata'
        fname += '.h'
        prefix = self._get_full_fname_prefix(dataname)
        tmp = self.dict_
        self._replace_dict(tmp, prefix, fname)

    def _generate_c_struct(self, data, dataname):
        fname = 'almdata'
        fname += '.c'
        prefix = self._get_full_fname_prefix(dataname)
        zeros_optvar_seqlen = self._np2Carray(
            np.zeros(data['optvar']['seqlen']))
        zeros_lagmul_seqlen = self._np2Carray(
            np.zeros(data['lagmul']['seqlen']))
        zeros_lagmul_optvar_seqlen = self._np2Carray(np.zeros(
            data['lagmul']['seqlen']*data['optvar']['seqlen']))

        tmp = dict(mu=data['mu'],
                   Linv=data['Linv'],
                   zeros_optvar_seqlen=zeros_optvar_seqlen,
                   zeros_lagmul_seqlen=zeros_lagmul_seqlen,
                   zeros_lagmul_optvar_seqlen=zeros_lagmul_optvar_seqlen,
                   lagmul_veclen=data['lagmul']['veclen'])
        tmp.update(self.dict_)
        self._replace_dict(tmp, prefix, fname)
