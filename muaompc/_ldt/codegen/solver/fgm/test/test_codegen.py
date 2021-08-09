#!/usr/bin/python3
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
from muaompc._ldt.codegen.solver.fgm import codegen
from muaompc._ldt.codegen.solver.fgm import codegency
#rom muaompc._ldt.codegen.solver.fgm import codegenmex

class TestConstantHessian(unittest.TestCase):

    def setUp(self):
        self.num = data
        self.destdir = 'testdestdir'
        self.prbnameref = 'regrefmpc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, self.prbnameref+'.prb')

        prsref = prbdsl.parse_file(prb_path)
        symref = ProblemSym(prsref)
        self.stcref = ProblemStruct(symref)
        self.cvpref = CondensedVectorialParameters(self.stcref)
        self.paramref = dict(x_k=np.array([2, 3]),
                xr=[np.array([5, 7]), np.array([5, 7]), np.array([11, 13])],
                ur=[np.array([17]), np.array([19])])
        #rb_path = os.path.join(base_dir, 'regorigmpc.prb')
        prsori = prbdsl.parse_file(prb_path)
        symori = ProblemSym(prsori)
        #et_trace()
        self.stcori = ProblemStruct(symori)
        self.cvpori = CondensedVectorialParameters(self.stcori)
        self.paramori = dict(x_k=np.array([2, 3]))

        shutil.rmtree(self.destdir, ignore_errors=True)
        self.basedir = shutil.os.path.abspath(shutil.os.curdir)

    def tearDown(self):
        del self.cvpref
        shutil.os.chdir(self.basedir)
        shutil.rmtree(self.destdir, ignore_errors=True)

    #unittest.skip('This is how you skip')
    def test_inc(self):
        prefix = 'mpcinc'
        Linv = 0.4985780638467639
        H = 2*Linv*np.array([[ 1.00157741, 0.00133198],
            [0.00133198, 1.00146001]])
        g = 2*Linv*np.array([[-18.28409196],
               [-20.31826494]])
        u_opt = np.array([ 18.22316,  20.26981])

        ccg = commoncodegen.CCodeGenerator(self.prbnameref, prefix, self.destdir)
        ccg.generate_code('float64_t')
        cdg = cvpcodegen.CCodeGenerator(self.cvpref, ccg)
        cdg.generate_code()
        cycdg = cvpcodegency.CythonCodeGenerator(cdg)
        cycdg.generate_code()
        slvccg = codegen.CCodeGenerator(ccg)
        tdtg = codegen.CDataGenerator(slvccg, cdg)
        tdtg.generate_data(self.num, prefix)  # test static data gen
        slv = codegen.CCodeGenerator(ccg)
        slv.generate_code()
        dtg = codegen.FGMCVPDataGenerator(cdg)
        data = dtg._get_data(self.num)
        cycdg = codegency.CythonCodeGenerator(slv)
        cycdg.generate_code()
        testdir = os.path.join(self.basedir, '%s/%s_%s/' % (
            self.destdir, prefix, self.prbnameref))
        shutil.os.chdir(testdir)
        path.append(testdir + '/src/cython/')
        cvpsetup = __import__(prefix+'cvpsetup')
        fgmsetup = __import__(prefix+'fgmsetup')
        setup(ext_modules=cythonize(fgmsetup.ext_modules+cvpsetup.ext_modules),
                script_args=['build_ext', '--inplace'])
        path.append(testdir)
        import mpcincfgm
        shutil.os.chdir('../../')
        dtg = codegen.FGMCVPDataGenerator(cdg)
        data = dtg._get_data(self.num)
        fgm = mpcincfgm.Solver()
        fgm.setup_solver(data)
        pardata = dict()
        pardata['xr'] = self._vecseq2array(self.paramref['xr'])
        pardata['ur'] = self._vecseq2array(self.paramref['ur'])
        pardata['x_k'] = np.matrix(self.paramref['x_k']).T
        fgm.configure(2, False)
        fgm.u_ini = np.array([-1000., 1000.])
        fgm.solve_problem(pardata)
        assert_allclose(fgm.u_opt, u_opt, rtol=1e-6,
                err_msg='Solution')
        pardata['x_k'] = np.matrix([[1., 2.]]).T
        fgm.solve_problem(pardata)
        assert_allclose(fgm.u_opt, np.array([18.35037369, 20.39036194]), rtol=1e-9,
                err_msg='Solution')

        fgm.configure(2, True)
        pardata['x_k'] = np.matrix(self.paramref['x_k']).T
        fgm.solve_problem(pardata)
        assert_allclose(fgm.u_opt, u_opt, rtol=1e-6,
                err_msg='Solution')
        pardata['x_k'] = np.matrix([[1., 2.]]).T
        fgm.solve_problem(pardata)
        assert_allclose(fgm.u_opt, np.array([18.35561064, 20.38488912]), rtol=1e-9,
                err_msg='Solution')

    def _vecseq2array(self, vecseq):
        outer_rows = len(vecseq)
        n = vecseq[0].shape[0]
        arr = np.zeros((n*outer_rows, 1))

        for i, mtx in enumerate(vecseq):
                arr[n*i:n*(i+1), 0] = mtx

        return arr

    @unittest.skip('This is how you skip')
    def test_mexinc(self):
        prefix = 'mpcinc'

        ccg = commoncodegen.CCodeGenerator(prefix, self.destdir)
        ccg.generate_code('float64_t')
        cdg = cvpcodegen.CCodeGenerator(self.cvpref, ccg)
        cdg.generate_code()
        cycdg = cvpcodegency.CythonCodeGenerator(cdg)
        cycdg.generate_code()
        tdtg = codegen.CDataGenerator(cdg)
        tdtg.generate_static_data(self.num, prefix)  # test static data gen
        slv = codegen.CCodeGenerator(ccg)
        slv.generate_code()
        dtg = codegen.FGMCVPDataGenerator(cdg)
        data = dtg._get_data(self.num)
        #excdg = codegenmex.MatlabCodeGenerator(slv)
        #excdg.generate_code()

if __name__ == '__main__':
    unittest.main()
