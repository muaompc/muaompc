"""Integration test for parsing of problem and data files."""

import unittest

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

    def test_dat_input_constr_ok(self):
        mpc = ldt.setup_mpc_problem('fixtures/test.prb')
        num = ldt._get_data(mpc, 'fixtures/test_input_constr_ok.dat', safe_mode=True)
        mpc.ddg.generate_data(num, 'xyz')

    def test_dat_state_constr_ok(self):
        mpc = ldt.setup_mpc_problem('fixtures/test_state_constr.prb')
        num = ldt._get_data( mpc, 'fixtures/test_state_constr_ok.dat', safe_mode=True)
        mpc.ddg.generate_data(num, 'xyz')

if __name__ == '__main__':
    unittest.main()
