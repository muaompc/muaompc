import os
import unittest
import shutil
from setuptools import setup
from sys import path

import numpy as np
from numpy.testing import assert_allclose
from Cython.Build import cythonize

from muaompc._ldt.parse.prbsym import ProblemSym
from muaompc._ldt.parse.prbstruct import ProblemStruct
from muaompc._ldt.parse import prbdsl
from muaompc._ldt.codegen.codegen import BaseDataGenerator
from muaompc._ldt.codegen.former.cvp.test.regmpcdata import data
import muaompc._ldt.codegen.former.cvp.test.regmpcdata as regmpcdata
from muaompc._ldt.codegen.former.cvp.cvp import CondensedVectorialParameters
from muaompc._ldt.codegen.common import codegen as commoncodegen
from muaompc._ldt.codegen.former.cvp import codegen as cvpcodegen
from muaompc._ldt.codegen.former.cvp import codegency as cvpcodegency
from muaompc._ldt.codegen.solver.alm import codegen
from muaompc._ldt.codegen.solver.alm import codegency

class TestConstantHessian(unittest.TestCase):

    def setUp(self):
        self.num = data
        self.destdir = 'testdestdir'
        self.prbnameorig = 'regorigmpc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, self.prbnameorig+'.prb')

        prsori = prbdsl.parse_file(prb_path)
        symori = ProblemSym(prsori)
        self.stcori = ProblemStruct(symori)
        self.cvpori = CondensedVectorialParameters(self.stcori)
        self.paramori = dict(x_k=np.array([2, 3]))

        self.basedir = shutil.os.path.abspath(shutil.os.curdir)
        shutil.rmtree(self.destdir, ignore_errors=True)

    def tearDown(self):
        shutil.os.chdir(self.basedir)
        shutil.rmtree(self.destdir, ignore_errors=True)

#   @unittest.skip('This is how you skip')
    def test_stc(self):
        prefix = 'mpcstc'
        H = 2*np.array([[ 1.00157741, 0.00133198],
            [0.00133198, 1.00146001]])
        g = 2*np.array([[0.33038257],
               [0.3054006]])
        V = np.array([[0.,  0.],
                    [0.10501746,  0.],
                    [0.01899904,  0.02095512]])
        u_ub = np.array([[100], [100]])
        v_ub = np.array([[-22.5], [-21.02973811], [ -2.61268764]])
        u_opt = np.array([-100., -12.749567])
        l_opt = np.array([-1184.000000,
                   28085.497984, 1800.198024])

        ccg = commoncodegen.CCodeGenerator(self.prbnameorig, prefix, self.destdir)
        ccg.generate_code('float64_t')
        cdg = cvpcodegen.CCodeGenerator(self.cvpori, ccg)
        cdg.generate_code()
        cycdg = cvpcodegency.CythonCodeGenerator(cdg)
        cycdg.generate_code()
        slvccg = codegen.CCodeGenerator(ccg)
        tdtg = codegen.CDataGenerator(slvccg, cdg)
        tdtg.generate_data(self.num, prefix)
        slv = codegen.CCodeGenerator(ccg)
        slv.generate_code()
        dtg = codegen.ALMCVPDataGenerator(cdg)
        data = dtg._get_data(self.num)
        cycdg = codegency.CythonCodeGenerator(slv)
        cycdg.generate_code()

        testdir = os.path.join(self.basedir, '%s/%s_%s/' % (
            self.destdir, prefix, self.prbnameorig))
        shutil.os.chdir(testdir)
        path.append(testdir + '/src/cython/')
        cvpsetup = __import__(prefix+'cvpsetup')
        almsetup = __import__(prefix+'almsetup')
        setup(ext_modules=cythonize(cvpsetup.ext_modules+almsetup.ext_modules),
                script_args=['build_ext', '--inplace'])
        path.append(testdir)
        import mpcstcalm
        shutil.os.chdir('../../')
        dtg = codegen.ALMCVPDataGenerator(cdg)
        data = dtg._get_data(self.num)
        alm = mpcstcalm.Solver()
        alm.setup_solver(data)
        pardata = dict()
        ex_iter = 3
        in_iter = 2
        alm.configure(ex_iter, in_iter, False)
        alm.u_ini = np.array([1000., -1000.])
        alm.l_ini = np.array([-20e3, 20e3, 1e3])
        pardata['x_k'] = np.matrix(self.paramori['x_k']).T
        alm.solve_problem(pardata)
        assert_allclose(alm.u_opt, u_opt, rtol=1e-6,
                err_msg='Solution')
        assert_allclose(alm.l_opt, l_opt, rtol=1e-6,
                err_msg='Solution')
        pardata['x_k'] = np.matrix([[1., 2.]]).T
        alm.solve_problem(pardata)
        assert_allclose(alm.u_opt, np.array([-100., -5.55513153]), rtol=1e-6,
                err_msg='Solution')
        assert_allclose(alm.l_opt, np.array([-8096., 21549.88502795, 1021.37388622]), rtol=1e-6,
                err_msg='Solution')

        alm.configure(ex_iter, in_iter, True)
        pardata['x_k'] = np.matrix(self.paramori['x_k']).T
        alm.solve_problem(pardata)
        assert_allclose(alm.u_opt, u_opt, rtol=1e-6,
                err_msg='Solution')
        assert_allclose(alm.l_opt, l_opt, rtol=1e-6,
                err_msg='Solution')
        pardata['x_k'] = np.matrix([[1., 2.]]).T

        alm.solve_problem(pardata)
        assert_allclose(alm.u_opt, np.array([-100., -0.14098062]), rtol=1e-6,
                err_msg='Solution')
        assert_allclose(alm.l_opt, np.array([ 38453.49798376, 3350.08305197, 0.]), rtol=1e-6,
                err_msg='Solution')

    def _vecseq2array(self, vecseq):
        outer_rows = len(vecseq)
        n = vecseq[0].shape[0]
        arr = np.zeros((n*outer_rows, 1))

        for i, mtx in enumerate(vecseq):
                arr[n*i:n*(i+1), 0] = mtx

        return arr


if __name__ == '__main__':
    unittest.main()
