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
        self.path = self._get_paths('codegen.solver.pbm')

    def generate_code(self):
        self._make_destination_dir_tree()
        self._generate_c_code()

    def _generate_c_code(self):
        self._replace_prefix('pbmsolve.h')
        self._replace_prefix('pbmsolve.c')
        self._generate_c_header('pbm')
        self._generate_c_body('pbm')
        self._generate_getdata_h()
        self._generate_getdata_c()

    def _generate_c_header(self, fname):
        fname += '.h'
        self._replace_prefix(fname)

    def _generate_c_body(self, fname):
        fname += '.c'
        self._replace_prefix(fname)

    def _generate_getdata_h(self):
        fname = 'pbmdynmem.h'
        self._replace_prefix(fname)

    def _generate_getdata_c(self):
        fname = 'pbmdynmem.c'
        self._replace_prefix(fname)


class PBMCVPDataGenerator(CVPDG, object):
    """
    Generate json data for the primal barrier method using a
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
        else:
            lenH = int(lenH)
        nu, L = self._compute_extra_step_constant(np.reshape(H, (lenH, lenH)))
        data = self._scale_cost(data, H, L)
        data['nu'] = nu
        data = self._get_pbm_constants(data)
        return data

    def _get_pbm_constants(self, data):
        data['kappa'] = 90.
        data['roh'] = 0.1
        data['optvar']['horizon'] = 1
        data['optvar']['state_veclen'] = 0
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


class CDataGenerator(CCodeGenerator, PBMCVPDataGenerator, object):

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
        self.path = self._get_paths('codegen.solver.pbm')
        self.cvpcg = cvpcg
        self.cvpcdg = CVPCDG(cvpcg)
        self.dict_ = None  # Assignment made during data generation

        self.is_sparse_data = 0  # TODO Set is_sparse_data.
        self.num = None

    def generate_data(self, num, dataname):
        self.num = num
        self.dict_ = dict(prefix=self.prefix,
                          PREFIX=self.prefix.upper(),
                          dataname=dataname,
                          DATANAME=dataname.upper())
        self.cvpcdg.dict_ = self.dict_
        pbmcvpdg = PBMCVPDataGenerator(self.cvpcg)
        data = pbmcvpdg._get_data(num)
        self.cvpcdg._write_static_data(data, dataname)
        self._write_static_data(data, dataname)

    def _write_static_data(self, data, dataname):
        self._generate_h_struct(data, dataname)
        self._generate_c_struct(data, dataname)

    def _generate_h_struct(self, data, dataname):
        fname = 'pbmdata'
        fname += '.h'
        prefix = self._get_full_fname_prefix(dataname)
        tmp = self.dict_
        self._replace_dict(tmp, prefix, fname)

    def _generate_c_struct(self, data, dataname):
        fname = 'pbmdata'
        fname += '.c'
        prefix = self._get_full_fname_prefix(dataname)

        v_array = np.array(data['pmetric']['V']['fac0']['data'])
        number_affine = 2 * data['pmetric']['V']['fac0']['rows']
        number_qc = 0  # TODO Add nb_qc
        number_soft = 0
        number_ineq = number_affine + number_qc + data['socc']['socc_num']
        if number_affine == 0:
            zeros_affine = {0}  # no affine constraints TODO
        else:
            zeros_affine = self._np2Carray(np.zeros(number_affine))
        p_data = self._np2Carray(np.hstack((-v_array, v_array)))

        if self.is_sparse_data == 1:
            # from a condensed problem, get the dimensions for a sparse problem
            # in particular, note that the sparse optvar is z = [u_0, x_1, ..., u_(N-1) x_N]
            state_veclen = data['par']['x_k']['rows']
            control_veclen = data['optvar']['veclen']
            pbm_horizon = int(data['optvar']['seqlen'] // data['optvar']['veclen'])
            pbm_optvar_veclen = state_veclen + control_veclen
            pbm_optvar_seqlen = pbm_optvar_veclen * pbm_horizon
            dual_seqlen = pbm_horizon*state_veclen
            zeros_dual_seqlen = self._np2Carray(np.zeros(dual_seqlen))
            zeros_state_veclen = self._np2Carray(np.zeros(state_veclen))
            zeros_state_veclen_x_optvar_veclen = self._np2Carray(
                np.zeros(state_veclen*pbm_optvar_veclen))
            zeros_state_veclen_x_control_veclen = self._np2Carray(
                np.zeros(state_veclen*control_veclen))
            zeros_state_veclen_x_state_veclen = self._np2Carray(
                np.zeros(state_veclen*state_veclen))
            zeros_2_x_horizon_x_state_veclen_x_state_veclen = self._np2Carray(
                np.zeros((2*pbm_horizon)*state_veclen*state_veclen))  # eigentlich reicht (2*pbm_horizon-1)
            zeros_optvar_seqlen_x_dual_seqlen = self._np2Carray(
                np.zeros(pbm_optvar_seqlen*dual_seqlen))

            a_data = self._np2Carray(np.array(self.num['A']).flatten())
            b_data = self._np2Carray(np.array(self.num['B']).flatten())

            # Equality constraints matrix
            c_array = np.zeros([pbm_horizon*state_veclen,
                                pbm_horizon*(state_veclen+control_veclen)])
            c_array[0:state_veclen, 0:control_veclen+state_veclen] =\
                    np.hstack([-np.array(self.num['B']),
                               np.eye(state_veclen, state_veclen)])
            for i in range(1, pbm_horizon):
                c_array[i*state_veclen:(i+1)*state_veclen, control_veclen+(i-1)*(control_veclen+state_veclen):control_veclen+i*(control_veclen+state_veclen)+state_veclen] = np.hstack([-np.array(self.num['A']), -np.array(self.num['B']), np.eye(state_veclen, state_veclen)])
            c_data = self._np2Carray(c_array.flatten())

            pbm_sparse = '#define HHMPC_CVP_PRB_SPARSE 1\n'

        else:
            pbm_horizon = data['optvar']['horizon']
            pbm_optvar_veclen = data['optvar']['seqlen']
            pbm_optvar_seqlen = data['optvar']['seqlen']
            state_veclen = data['optvar']['state_veclen']
            control_veclen = data['optvar']['veclen']
            dual_seqlen = 0
            zeros_dual_seqlen = {0}  # ISO C forbids zero-size array
            zeros_state_veclen = {0}
            zeros_state_veclen_x_state_veclen = {0}
            zeros_state_veclen_x_optvar_veclen = {0}
            zeros_state_veclen_x_control_veclen = {0}
            zeros_2_x_horizon_x_state_veclen_x_state_veclen = {0}
            zeros_optvar_seqlen_x_dual_seqlen = {0}

            a_data = {0}
            b_data = {0}
            c_data = {0}

            pbm_sparse = ''
            pbm_sparse += '#define HHMPC_CVP_PRB_SPARSE 0\n'

        zeros_optvar_seqlen = self._np2Carray(np.zeros(pbm_optvar_seqlen))
        zeros_optvar_seqlen_x_optvar_seqlen = self._np2Carray(
            np.zeros(pbm_optvar_seqlen*pbm_optvar_seqlen))
        if number_ineq == 0:
            zeros_ineq = {0}
            zeros_ineq_x_ineq = {0}
            zeros_optvar_seqlen_x_number_ineq = {0}
        else:
            zeros_ineq = self._np2Carray(np.zeros(number_ineq))
            zeros_ineq_x_ineq = self._np2Carray(
                np.zeros(number_ineq*number_ineq))
            zeros_optvar_seqlen_x_number_ineq = self._np2Carray(
                np.zeros(pbm_optvar_seqlen*number_ineq))
        if number_soft == 0:
            zeros_soft = {0}
            zeros_soft_x_soft = {0}
            zeros_optvar_seqlen_x_soft = {0}
        else:
            zeros_soft = self._np2Carray(np.zeros(number_soft))
            zeros_soft_x_soft = self._np2Carray(np.zeros(number_soft*number_soft))
            zeros_optvar_seqlen_x_soft = self._np2Carray(np.zeros(pbm_optvar_seqlen*number_soft))


        pbm_zeros_optvar_veclen_x_optvar_veclen = self._np2Carray(
                np.zeros(pbm_optvar_veclen*pbm_optvar_veclen))
        pbm_zeros_horizon_x_optvar_veclen_x_optvar_veclen = self._np2Carray(
                np.zeros(pbm_horizon*pbm_optvar_veclen*pbm_optvar_veclen))

        tmp = dict(nu=data['nu'],
                   zeros_optvar_seqlen=zeros_optvar_seqlen,
                   zeros_dual_seqlen=zeros_dual_seqlen,
                   zeros_state_veclen=zeros_state_veclen,
                   zeros_affine=zeros_affine,
                   zeros_ineq=zeros_ineq,
                   zeros_soft=zeros_soft,
                   zeros_soft_x_soft=zeros_soft_x_soft,
                   zeros_optvar_seqlen_x_soft=zeros_optvar_seqlen_x_soft,
                   zeros_ineq_x_ineq=zeros_ineq_x_ineq,
                   zeros_optvar_seqlen_x_optvar_seqlen=zeros_optvar_seqlen_x_optvar_seqlen,
                   zeros_optvar_seqlen_x_number_ineq=zeros_optvar_seqlen_x_number_ineq,
                   zeros_state_veclen_x_state_veclen=zeros_state_veclen_x_state_veclen,
                   zeros_state_veclen_x_optvar_veclen=zeros_state_veclen_x_optvar_veclen,
                   zeros_state_veclen_x_control_veclen=zeros_state_veclen_x_control_veclen,
                   pbm_zeros_optvar_veclen_x_optvar_veclen=pbm_zeros_optvar_veclen_x_optvar_veclen,
                   pbm_zeros_horizon_x_optvar_veclen_x_optvar_veclen=pbm_zeros_horizon_x_optvar_veclen_x_optvar_veclen,
                   zeros_2_x_horizon_x_state_veclen_x_state_veclen=zeros_2_x_horizon_x_state_veclen_x_state_veclen,
                   zeros_optvar_seqlen_x_dual_seqlen=zeros_optvar_seqlen_x_dual_seqlen,
                   pbm_p_data=p_data,
                   pbm_a_data=a_data,
                   pbm_b_data=b_data,
                   pbm_c_data=c_data,
                   pbm_horizon=pbm_horizon,
                   state_veclen=state_veclen,
                   control_veclen=control_veclen,
                   pbm_optvar_veclen=pbm_optvar_veclen,
                   optvar_seqlen=pbm_optvar_seqlen,
                   dual_seqlen=dual_seqlen,
                   number_affine=number_affine,
                   number_socc=data['socc']['socc_num'],
                   number_ineq=number_ineq,
                   pbm_sparse=pbm_sparse)

        tmp.update(self.dict_)
        self._replace_dict(tmp, prefix, fname)
