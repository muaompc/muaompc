# test for codegen.py and codegency.py
import os
import unittest
import shutil
from setuptools import setup
from sys import path
from importlib import import_module

import numpy as np
from numpy.testing import assert_allclose
from Cython.Build import cythonize

from muaompc._ldt.parse.prbsym import ProblemSym
from muaompc._ldt.parse.prbstruct import ProblemStruct
from muaompc._ldt.parse import prbdsl
from muaompc._ldt.codegen.former.cvp.test.regmpcdata import data
from muaompc._ldt.codegen.former.cvp.cvp import CondensedVectorialParameters
from muaompc._ldt.codegen.common import codegen as commoncodegen
from muaompc._ldt.codegen.former.cvp import codegen
from muaompc._ldt.codegen.former.cvp import codegency


class TestConstantHessian(unittest.TestCase):

    def setUp(self):
        self.num = data
        self.destdir = 'testdestdir'
        self.prefix = 'mpc'
        self.base_dir = os.path.dirname(__file__)

        self.paramref = dict(x_k=np.array([2, 3]),
                             xr=[np.array([5, 7]), np.array(
                                 [5, 7]), np.array([11, 13])],
                             ur=[np.array([17]), np.array([19])])

        self.paramori = dict(x_k=np.array([2, 3]))

        self.prbnamestoc = 'stocmpc'
        prb_path = os.path.join(self.base_dir, self.prbnamestoc+'.prb')
        prsstoc = prbdsl.parse_file(prb_path)
        symstoc = ProblemSym(prsstoc)
        self.stcstoc = ProblemStruct(symstoc)
        self.cvpstoc = CondensedVectorialParameters(self.stcstoc)
        self.prbnameoutref = 'regoutrefmpc'

        self.prbnameref = 'regrefmpc'
        prb_path = os.path.join(self.base_dir, self.prbnameref+'.prb')
        prsref = prbdsl.parse_file(prb_path)
        symref = ProblemSym(prsref)
        self.stcref = ProblemStruct(symref)
        self.cvpref = CondensedVectorialParameters(self.stcref)
        shutil.rmtree(self.destdir, ignore_errors=True)
        self.basedir = shutil.os.path.abspath(shutil.os.curdir)

    def tearDown(self):
        del self.cvpref
        shutil.rmtree(self.destdir, ignore_errors=True)
        shutil.os.chdir(self.basedir)

    @unittest.skip('This is how you skip')
    def test_socc(self):
        prefix = 'stocmpc'
        H = 2*np.array([[1.00157741, 0.00133198],
                        [0.00133198, 1.00146001]])
        g = 2*np.array([[0.33038257],
                        [0.3054006]])
        V = np.array([[0.,  0.],
                      [0.10501746,  0.],
                      [0.01899904,  0.02095512]])
        Wm0 = np.array([[9.67483607e-05,   0.00000000e+00],
                        [1.90325164e-02,   0.00000000e+00]])
        wn0 = np.array([[2.02854877],
                        [2.71451225]])
        wvT0 = np.array([[-9.67483607e-05,   0.00000000e+00]])
        ws0 = np.array([[-2.37754877]])
        Wm1 = np.array([[2.77866701e-04,   9.67483607e-05],
                        [1.72213330e-02,   1.90325164e-02]])
        wn1 = np.array([[2.05438077],
                        [2.45619226]])
        wvT1 = np.array([[-2.77866701e-04,  -9.67483607e-05]])
        ws1 = np.array([[-2.40338077]])

        u_ub = np.array([[100], [100]])
        v_ub = np.array([[-22.5], [-21.02973811], [-2.61268764]])

        ccg = commoncodegen.CCodeGenerator(
            self.prbnamestoc, prefix, self.destdir)
        ccg.generate_code('float64_t')
        cdg = codegen.CCodeGenerator(self.cvpstoc, ccg)
        cdg.generate_code()
        cycdg = codegency.CythonCodeGenerator(cdg)
        cycdg.generate_code()
        tdtg = codegen.CDataGenerator(cdg)
        tdtg.generate_data(self.num, prefix)
        shutil.os.chdir('%s/%s_%s/' % (
            self.destdir, self.prbnamestoc, prefix))
        path.append('./src/cython/')
        cvpsetup = __import__(prefix+'cvpsetup')
        setup(ext_modules=cythonize(cvpsetup.ext_modules),
              script_args=['build_ext', '--inplace'])
        import stocmpccvp
        shutil.os.chdir('../../')
        qpx = stocmpccvp.Former()
        dtg = codegen.CVPDataGenerator(cdg)
        data = dtg._get_data(self.num)
        qpx.initialize(data)
        pardata = dict()
        pardata['x_k'] = np.matrix(self.paramref['x_k']).T

        qp = qpx.form_problem(pardata)

        assert_allclose(qp['H'], H, rtol=1e-6,
                        err_msg='Hessian matrix')

        assert_allclose(qp['g'], g, rtol=1e-6,
                        err_msg='gradient vector')

        assert_allclose(qp['u_ub'], u_ub, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['socc'][0]['Wm'], Wm0, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['socc'][0]['wvT'], wvT0, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['socc'][0]['wn'], wn0, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['socc'][0]['ws'], ws0, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['socc'][1]['Wm'], Wm1, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['socc'][1]['wvT'], wvT1, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['socc'][1]['wn'], wn1, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['socc'][1]['ws'], ws1, rtol=1e-6,
                        err_msg='Constraint matrix')

    #unittest.skip('This is how you skip')

    def test_H_mtx(self):
        # the original H was computed by muaompc 0.3.0, which uses half of the cost
        H = 2*np.array([[1.00157741, 0.00133198],
                        [0.00133198, 1.00146001]])
        prb = self.cvpref.form_condensed_prb(self.num)
        assert_allclose(prb['pmetric']['H']['fac0'], H, rtol=1e-6,
                        err_msg='Hessian matrix')

    #unittest.skip('This is how you skip')
    def test_inc(self):
        prefix = 'mpcinc'
        prbname = 'regrefmpc'
        H = 2*np.array([[1.00157741, 0.00133198],
                        [0.00133198, 1.00146001]])
        g = 2*np.array([[-18.28409196],
                        [-20.31826494]])

        cvpcy = self._get_cython_former(prbname, prefix)

        pardata = dict()
        pardata['xr'] = self._vecseq2array(self.paramref['xr'])
        pardata['ur'] = self._vecseq2array(self.paramref['ur'])
        pardata['x_k'] = np.matrix(self.paramref['x_k']).T

        qp = cvpcy.form_problem(pardata)

        assert_allclose(qp['H'], H, rtol=1e-6,
                        err_msg='Hessian matrix')

        assert_allclose(qp['g'], g, rtol=1e-6,
                        err_msg='gradient vector')

    #unittest.skip('This is how you skip')
    def test_stc(self):
        prefix = 'mpcstc'
        prbname = 'regorigmpc'
        H = 2*np.array([[1.00157741, 0.00133198],
                        [0.00133198, 1.00146001]])
        g = 2*np.array([[0.33038257],
                        [0.3054006]])
        V = np.array([[0.,  0.],
            [0.10501746,  0.],
            [0.01899904,  0.02095512]])
        u_ub = np.array([[100], [100]])
        v_ub = np.array([[-22.5],
            [-21.02973811], [-2.61268764]])

        cvpcy = self._get_cython_former(prbname, prefix)

        pardata = dict()
        pardata['x_k'] = np.matrix(self.paramori['x_k']).T

        qp = cvpcy.form_problem(pardata)

        assert_allclose(qp['H'], H, rtol=1e-6,
                        err_msg='Hessian matrix')

        assert_allclose(qp['g'], g, rtol=1e-6,
                        err_msg='gradient vector')

        assert_allclose(qp['V'], V, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['u_ub'], u_ub, rtol=1e-6,
                        err_msg='Constraint matrix')
        assert_allclose(qp['v_ub'], v_ub, rtol=1e-6,
                        err_msg='Constraint matrix')

    @unittest.skip('This is how you skip')
    def test_outref(self):
        from pudb import set_trace
        set_trace()
        prefix = 'mpcoutref'
        prbname = 'regoutrefmpc'
        H = 2*np.array([[1.00157741, 0.00133198],
                        [0.00133198, 1.00146001]])
        g = 2*np.array([[-18.28409196],
                        [-20.31826494]])

        self.num.update(Cx=np.eye(2))
        cvpcy = self._get_cython_former(prbname, prefix)

        pardata = dict()
        pardata['yr'] = self._vecseq2array(self.paramref['xr'])
        pardata['x_k'] = np.matrix(self.paramref['x_k']).T

        qp = cvpcy.form_problem(pardata)

        assert_allclose(qp['H'], H, rtol=1e-6,
                        err_msg='Hessian matrix')

        assert_allclose(qp['g'], g, rtol=1e-6,
                        err_msg='gradient vector')

    def _vecseq2array(self, vecseq):
        outer_rows = len(vecseq)
        n = vecseq[0].shape[0]
        arr = np.zeros((n*outer_rows, 1))

        for i, mtx in enumerate(vecseq):
            arr[n*i:n*(i+1), 0] = mtx

        return arr

    def _get_problem_cvp(self, prbname):
        prb_path = os.path.join(self.base_dir, prbname+'.prb')
        prsref = prbdsl.parse_file(prb_path)
        symref = ProblemSym(prsref)
        stc = ProblemStruct(symref)
        cvp = CondensedVectorialParameters(stc)
        return cvp

    def _generate_code_and_data(self, prbname, prefix):
        cvp = self._get_problem_cvp(prbname)
        print(self.destdir, prefix+prbname)
        ccg = commoncodegen.CCodeGenerator(prbname, prefix, self.destdir)
        ccg.generate_code('float64_t')
        cdg = codegen.CCodeGenerator(cvp, ccg)
        cdg.generate_code()
        cycdg = codegency.CythonCodeGenerator(cdg)
        cycdg.generate_code()
        tdtg = codegen.CDataGenerator(cdg)
        tdtg.generate_data(self.num, prefix)
        dtg = codegen.CVPDataGenerator(cdg)
        data = dtg._get_data(self.num)
        return data

    def _compile_cython_interface(self, prbname, prefix):
        path.append('./src/cython/')
        cvpsetup = import_module(prefix+'cvpsetup')
        setup(ext_modules=cythonize(cvpsetup.ext_modules),
              script_args=['build_ext', '--inplace'])
        return

    def _get_cython_former(self, prbname, prefix):
        data = self._generate_code_and_data(prbname, prefix)
        testdir = os.path.join(self.basedir, '%s/%s_%s/' % (
            self.destdir, prefix, prbname))
        print(testdir)
        shutil.os.chdir(testdir)
        self._compile_cython_interface(prbname, prefix)

        #print('\nCurrent path: ', shutil.os.path.abspath(shutil.os.curdir))
        path.append(testdir)
        modname = prefix + 'cvp'
        mpccvp = import_module(modname)
        shutil.os.chdir('../../')
        qpx = mpccvp.Former()
        qpx.initialize(data)
        return qpx

    @unittest.skip('This is how you skip')
    def test_mex(self):
        prefix = 'mpcinc'
        H = 2*np.array([[1.00157741, 0.00133198],
                        [0.00133198, 1.00146001]])
        g = 2*np.array([[-18.28409196],
                        [-20.31826494]])

        ccg = commoncodegen.CCodeGenerator(prefix, '.')
        ccg.generate_code('float64_t')
        cdg = codegen.CCodeGenerator(self.cvpori, ccg)
        cdg.generate_code()
        mexcdg = codegenmex.MatlabCodeGenerator(cdg)
        mexcdg.generate_code()


if __name__ == '__main__':
    unittest.main()
