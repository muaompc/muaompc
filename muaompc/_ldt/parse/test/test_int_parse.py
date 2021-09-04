"""Integration test for parsing of problem and data files."""

import sys

import unittest

import numpy as np
import pandas as pd
import json
from numpy.core.numeric import allclose
from scipy.linalg import block_diag
from numpy.linalg import matrix_power

from muaompc._ldt.parse.prbstruct import DataError
from muaompc import ldt

class TestParseDimensions(unittest.TestCase):
    def test_dat_bad_A(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        num = ldt._get_data(mpc, 'fixtures/test_bad_A.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_bad_Q(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        self.assertRaises(ValueError, ldt._get_data, mpc, 'fixtures/test_bad_Q.dat', safe_mode=True)

    def test_dat_R_not_sq(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        num = ldt._get_data( mpc, 'fixtures/test_R_not_sq.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_dat_R_too_big(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        num = ldt._get_data( mpc, 'fixtures/test_R_too_big.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_dat_Q_too_small(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        num = ldt._get_data( mpc, 'fixtures/test_Q_too_small.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_dat_n_too_big(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        num = ldt._get_data( mpc, 'fixtures/test_n_too_big.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_dat_m_too_big(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        num = ldt._get_data( mpc, 'fixtures/test_m_too_big.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_dat_u_bound_size(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        num = ldt._get_data( mpc, 'fixtures/test_u_bound_size.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_dat_state_bound_short(self):
        mpc = ldt.setup_mpc_problem('fixtures/test_state_constr.prb')
        num = ldt._get_data( mpc, 'fixtures/test_state_bound_short.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_dat_state_bound_long(self):
        mpc = ldt.setup_mpc_problem('fixtures/test_state_constr.prb')
        num = ldt._get_data( mpc, 'fixtures/test_state_bound_long.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')

    def test_dat_bad_Kx(self):
        mpc = ldt.setup_mpc_problem('fixtures/test_state_constr.prb')
        num = ldt._get_data( mpc, 'fixtures/test_bad_Kx.dat', safe_mode=True)
        self.assertRaises(DataError, mpc.ddg.generate_data, num, 'xyz')


class TestParseData(unittest.TestCase):
    def setUp(self):
        self.fname = 'xyz'
        self.prefix = 'mpc'
        self.fixtures = 'fixtures/'

    def tearDown(self):
        pass

    def test_dat_input_constr_ok(self):
        prbname = 'test'
        prbpath = self.fixtures + prbname + '.prb'
        mpc = ldt.setup_mpc_problem(prbpath)
        num = ldt._get_data(mpc, 'fixtures/test_input_constr_ok.dat', safe_mode=True)
        mpc.ddg.generate_data(num, self.fname)
        (H, G) = self._load_H_G(prbname)
        (HH, GG) = self._form_H_G(num)
        scale = max(np.linalg.eigvals(HH))
        np.testing.assert_allclose(HH, H*scale) 
        np.testing.assert_allclose(GG, G*scale)

    def test_dat_state_constr_ok(self):
        prbname = 'test_state_constr'
        prbpath = self.fixtures + prbname + '.prb'
        mpc = ldt.setup_mpc_problem(prbpath)
        num = ldt._get_data(mpc, 'fixtures/test_state_constr_ok.dat', safe_mode=True)
        mpc.ddg.generate_data(num, self.fname)
        (H, G) = self._load_H_G(prbname)
        (HH, GG) = self._form_H_G(num)
        np.testing.assert_allclose(HH, H) 
        np.testing.assert_allclose(GG, G)

    def test_dat_path_following(self):
        sys.path.append('fixtures')
        prbname = 'test_path_following'
        prbpath = self.fixtures + prbname + '.prb'
        mpc = ldt.setup_mpc_problem(prbpath)
        num = ldt._get_data(mpc, 'test_path_following.py', safe_mode=False)
        mpc.ddg.generate_data(num, self.fname)
        (H, G) = self._load_H_G(prbname, fac=1)  # get G: g(z_k) = G*z_k
        (HH, GG) = self._form_H_G_pf(num)
        np.testing.assert_allclose(HH, H) 
        np.testing.assert_allclose(GG, G)

    def _form_H_G_pf(self, num):
        # Hessian of path following
        N = num['N']
        QQ = self._form_big_Q(num['Q'], num['P'], N)
        RR = self._form_big_R(num['R'], N)
        AAx = self._form_big_A(num['Ax'], N)
        AAz = self._form_big_A(num['Az'], N)
        BBx = self._form_big_B(num['Ax'], num['Bx'], N)
        BBz = self._form_big_B(num['Az'], num['Bz'], N)
        BB = BBx - BBz

        H = BB.T @ QQ @ BB + RR
        # we are only interested in testing G: g(z_k) = G*z_k
        G = BB.T @ QQ @ AAz * -1  # |xi - zi|^2_Q: i.e. -1 * zi 
        return (H, G)

    def _form_H_G(self, num):
        N = num['N']
        AA = self._form_big_A(num['A'], N)
        BB = self._form_big_B(num['A'], num['B'], N)
        QQ = self._form_big_Q(num['Q'], num['P'], N)
        RR = self._form_big_R(num['R'], N)

        H = BB.T @ QQ @ BB + RR
        G = BB.T @ QQ @ AA
        return (H, G)

    def _form_big_A(self, A, N):
        n = A.shape[0]
        AA = np.zeros([n*N, n])
        for i in range(1, N+1):
            AA[(i-1)*n:i*n,:] = np.linalg.matrix_power(A, i)
        return AA

    def _form_big_B(self, A, B, N):
        (n, m) = B.shape
        BB = np.zeros([n*N, m*N])
        BB[:n, :m] = B[:]
        for i in range(1, N):
            AB = BB[(i-1)*n:i*n, 0:m]
            AAB = A @ AB
            BB[i*n:(i+1)*n, 0:m] = AAB[:]
            BBi = BB[(i-1)*n:i*n, 0:-m]
            BB[i*n:(i+1)*n, m:] = BBi

        return BB

    def _form_big_Q(self, Q, P, N):
        qq = [Q]*(N-1)
        qq.append(P)
        QQ = block_diag(*qq)
        return QQ

    def _form_big_R(self, R, N):
        RR = block_diag(*[R]*N)
        return RR

    def _load_H_G(self, prbname, fac=0):
        # load the Hessian matrix as generated by muaompc
        path = '%s_%s/data/%s/%s%s.json' % (self.prefix, prbname, self.fname,
        self.prefix, self.fname)
        with open(path) as f:
            d = json.load(f)

        Hp = d['pmetric']['H']['fac0']
        # get one of the parametric matrices Gp: g(pk) = g0 + G0*fac0 + G1*fac1 + ...
        Gp = d['pmetric']['g']['fac'][fac]

        H = np.array(Hp['data'])
        H = H.reshape((Hp['rows'], Hp['cols']))
        G = np.array(Gp['data'])
        G = G.reshape((Gp['rows'], Gp['cols']))
        # muaompc returns the cost function's Hessian matrix divided by a factor
        try:
            scale = d['Linv']*2  # state constraints problems (ALM+FGM)
        except KeyError:
            scale = 1.  # H has been divided by max eigenval (FGM)
        return (H/scale, G/scale)


if __name__ == '__main__':
    unittest.main()
