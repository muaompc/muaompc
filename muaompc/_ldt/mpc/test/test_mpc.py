import os
import sys
import unittest
import shutil
import subprocess
from subprocess import call
from importlib import import_module

import numpy as np
from numpy.testing import assert_allclose

from muaompc.ldt import setup_mpc_problem, generate_mpc_data
from muaompc._ldt.codegen.former.cvp.test.regmpcdata import data
import muaompc._ldt.codegen.former.cvp.test.regmpcdata as regmpcdata

class TestMPC(unittest.TestCase):

    def setUp(self):
        self.num = data
        self.destdir = 'testdestdir'
        self.python = sys.executable
        self.basedir = shutil.os.path.abspath(shutil.os.curdir)

        self.paramref = dict(x_k=np.array([2, 3]),
        xr=[np.array([5, 7]), np.array([5, 7]), np.array([11, 13])],
        ur=[np.array([17]), np.array([19])])
        self.paramori = dict(x_k=np.array([2, 3]))
        shutil.rmtree(self.destdir, ignore_errors=True)

    def tearDown(self):
        shutil.rmtree(self.destdir, ignore_errors=True)
        shutil.os.chdir(self.basedir)

    def test_ctl_input_constraints(self):
        # Test that the condensed problem is accessible through 
        # the controller interface
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        data_name = 'regmpc'
        mod_name = 'mpcctlfgm'  # module name
        dat_path = os.path.join(base_dir, data_name+'.dat')
        mpc = setup_mpc_problem(prb_path, mod_name, destdir=self.destdir)
        generate_mpc_data(mpc, dat_path)
        testdir = os.path.join(self.basedir, '%s/%s_%s/' % (
            self.destdir, mod_name, 'regrefmpc'))
        shutil.os.chdir(testdir)
        call([self.python, mod_name + 'setup.py', 'install'])
        from mpcctlfgm import mpcctlfgmctl as ctl
        c = ctl.Ctl('data/%s/%s%s.json' % (data_name, mod_name, data_name))
        assert_allclose([[0.99872905, 0.00132819], [0.00132819, 0.99861199]], c.prb.H, rtol=1e-6)
        assert_allclose([[100.],[100]], c.prb.u_ub)
        assert_allclose([[-100.],[-100]], c.prb.u_lb)
        assert_allclose([[0.],[0]], c.prb.g)
        c.parameters.x_k[:] = self.paramref['x_k']
        c.parameters.xr[:] = self._vecseq2array(self.paramref['xr'])[:,0]
        c.parameters.ur[:] = self._vecseq2array(self.paramref['ur'])[:,0]
        c.solve_problem()
        # before solving the problem, the parametric terms (like g) are updated
        assert_allclose([[-18.23209433],[-20.26048239]], c.prb.g, rtol=1e-6)
        shutil.os.chdir('../..')

    #unittest.skip('This is how you skip')
    def test_fgm_solver(self):
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        setup_mpc_problem(prb_path, 'mpcfgm', destdir=self.destdir)
        u_opt = np.array([18.22316, 20.26981])
        testdir = os.path.join(self.basedir, '%s/%s_%s/' % (
            self.destdir, 'mpcfgm', 'regrefmpc'))
        shutil.os.chdir(testdir)

        call([self.python, 'mpcfgmsetup.py', 'install'])
        shutil.os.chdir('../..')
        from mpcfgm.mpc import setup_mpc_ctl
        ctl = setup_mpc_ctl(data)
        pardata = dict()
        pardata['xr'] = self._vecseq2array(self.paramref['xr'])
        pardata['ur'] = self._vecseq2array(self.paramref['ur'])
        pardata['x_k'] = np.matrix(self.paramref['x_k']).T
        ctl.configure(2, False)
        ctl.u_ini = np.array([-1000., 1000.])
        ctl.solve_problem(pardata)
        assert_allclose(ctl.u_opt, u_opt, rtol=1e-6,
                err_msg='Solution')

    #unittest.skip('This is how you skip')
    def test_alm_solver(self):
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regorigmpc.prb')
        setup_mpc_problem(prb_path, 'mpcalm', destdir=self.destdir)
        u_opt = np.array([-100., -12.749567])
        l_opt = np.array([-1184.000000,
                   28085.497984, 1800.198024])
        testdir = os.path.join(self.basedir, '%s/%s_%s/' % (
            self.destdir, 'mpcalm', 'regorigmpc'))
        shutil.os.chdir(testdir)
        print("PYTHON: ", self.python)
        call([self.python, 'mpcalmsetup.py', 'install'])
        shutil.os.chdir('../..')
        from mpcalm.mpc import setup_mpc_ctl
        ctl = setup_mpc_ctl(data)
        pardata = dict()
        pardata['x_k'] = np.matrix(self.paramori['x_k']).T
        ctl.configure(3, 2, False)
        ctl.u_ini = np.array([1000., -1000.])
        ctl.l_ini = np.array([-20e3, 20e3, 1e3])
        ctl.solve_problem(pardata)
        assert_allclose(ctl.u_opt, u_opt, rtol=1e-6,
                err_msg='Solution')
        assert_allclose(ctl.l_opt, l_opt, rtol=1e-6,
                err_msg='Solution')

    @unittest.skip('This is how you skip')
    def test_fgm_solver_static_data(self):
        import muaompc
        muaompc_dir = os.path.abspath(os.path.dirname(muaompc.__file__))

        prefix = 'mpcfgm'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        data_name = 'regmpc'
        dat_path = os.path.join(base_dir, data_name+'.dat')
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=self.destdir)
        generate_mpc_data(mpc, dat_path)
        #shutil.os.chdir('%s/%s_%s/src' % (self.destdir, 'regrefmpc', 'mpcfgm'))
        #shutil.os.chdir('../../..')
        shutil.copy(os.path.join(muaompc_dir + '/_ldt/mpc/test/mainfgm.c'), 'main.c')
        shutil.os.remove(prefix+data_name+'main.c')
        call(['make', '-f', prefix+data_name+'Makefile.mk', 'clean'])
        call(['make', '-f', prefix+data_name+'Makefile.mk'])
        args = list(self.paramref['x_k']) + list(self.paramref['xr'][0]) + list(
                self.paramref['xr'][1]) + list(self.paramref['xr'][2]) + list(
                        self.paramref['ur'][0]) + list(
                self.paramref['ur'][1]) + [-1000., +1000.]
        p = subprocess.Popen(['./main'] + [str(arg) for arg in args],
                stdout=subprocess.PIPE)
        shutil.os.chdir('../../../..')
        sol = p.stdout.readline()
        if isinstance(sol, str):  # probably Python2
            ctl_sol = np.matrix(sol)
        elif isinstance(sol, bytes):  # probably Python3
            ctl_sol = np.matrix(sol.decode())
        ctl_u_opt = np.array(ctl_sol[0, 0:2]).flatten()

        u_opt = np.array([18.22316, 20.26981])
        assert_allclose(ctl_u_opt, u_opt, rtol=1e-6,
                err_msg='Solution')

    def test_alm_solver_static_data(self):
        import muaompc
        muaompc_dir = os.path.abspath(os.path.dirname(muaompc.__file__))

        prefix = 'mpcalm'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regorigmpc.prb')
        data_name = 'regmpc'
        dat_path = os.path.join(base_dir, data_name+'.dat')
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=self.destdir)
        generate_mpc_data(mpc, dat_path)
        shutil.os.chdir(shutil.os.path.abspath(shutil.os.curdir))
        shutil.os.chdir('%s/%s_%s/src' % (self.destdir, 'mpcalm', 'regorigmpc'))
        #shutil.copy(os.path.join(muaompc_dir + '/_ldt/codegen/common/test/math/static_data.h'), '.')
        shutil.os.chdir('../../..')
        testdir = os.path.join(self.basedir, '%s/%s_%s/data/%s' % (
            self.destdir, 'mpcalm', 'regorigmpc', 'regmpc'))
        shutil.os.chdir(testdir)
        #hutil.copy('/home/zmenende/workspace/muaompc/muaompc/_ldt/mpc/test/main.c', '.')
        shutil.copy(os.path.join(muaompc_dir + '/_ldt/mpc/test/main.c'), '.')
        shutil.os.remove(prefix+data_name+'main.c')
        call(['make', '-f', prefix+data_name+'Makefile.mk', 'clean'])
        call(['make', '-f', prefix+data_name+'Makefile.mk'])
        args = list(self.paramori['x_k']) + (
                [1000., -1000.] + [-20e3, 20e3, 1e3])
        p = subprocess.Popen(['./main'] + [str(arg) for arg in args],
                stdout=subprocess.PIPE)
        shutil.os.chdir('../../../..')
        sol = p.stdout.readline()
        if isinstance(sol, str):  # probably Python2
            ctl_sol = np.matrix(sol)
        elif isinstance(sol, bytes):  # probably Python3
            ctl_sol = np.matrix(sol.decode())
        ctl_u_opt = np.array(ctl_sol[0, 0:2]).flatten()
        ctl_l_opt = np.array(ctl_sol[0, 2:5]).flatten()

        u_opt = np.array([-100., -12.749567])
        l_opt = np.array([-1184.000000,
                   28085.497984, 1800.198024])
        assert_allclose(ctl_u_opt, u_opt, rtol=1e-6,
                err_msg='Solution')
        assert_allclose(ctl_l_opt, l_opt, rtol=1e-6,
                err_msg='Solution')

    @unittest.skip('This is how you skip')
    def test_pbm_solver_cond_data(self):
        import muaompc
        muaompc_dir = os.path.abspath(os.path.dirname(muaompc.__file__))
        prefix = 'condhhmpc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regorigmpc.prb')  #gen cvp for alm
        import muaompc._ldt.codegen.former.cvp.test.pce as pce
        data_name = 'pce'
        #at_path = os.path.join(base_dir, data_name)
        #at_path = 'pce'
        destdir = self.destdir
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=destdir, 
                solver='pbm')
        generate_mpc_data(mpc, pce, safe_mode=False)
        shutil.os.chdir('%s/%s_%s/src' % (destdir, 'regorigmpc', 'condhhmpc'))
        # shutil.copy(os.path.join(
        #     muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm_cond/static_data.h'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm_cond/condhhmpcsocp.h'), '.')
        shutil.os.chdir('../../..')
        shutil.os.chdir('%s/%s_%s/data/%s' % (
            destdir, 'regorigmpc', 'condhhmpc', 'pce'))
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm_cond/main.c'), '.')
        shutil.os.remove(prefix+data_name+'main.c')
        call(['make', '-f', prefix+data_name+'Makefile.mk', 'clean'])
        call(['make', '-f', prefix+data_name+'Makefile.mk'])
        u_ref = np.array([-0.226552284045, 0.102580066361, -0.142418641446,
                          0.055273441036, -0.061111695401, 0.028558669224,
                          -0.025119132617, 0.013521211350, -0.000750050738,
                          0.010381000188])
        u_ref = np.array([-0.226519736693, 0.102530448633, -0.142376565305,
                          0.055243190007, -0.061088108164, 0.028803444911,
                          -0.025489050303, 0.013783357019, -0.000891999469,
                          0.010396940561])
        # vorrübergehend, weil H und g anders sind und dadurch auch kappa
        u_ref = np.array([-0.2280065 ,  0.1046547 , -0.1447085 ,  0.05744459, -0.06321969,
        0.03063431, -0.02772489,  0.01514494, -0.00699698,  0.00845786])

        p = subprocess.Popen(['./main'],# + [str(arg) for arg in args],
                stdout=subprocess.PIPE)
        # os.remove('ipm_data.c')
        # os.remove('socp_data.c')
        # os.remove('hhmpcsocp.c')
        # os.remove('hhmpcsocp.h')
        shutil.os.chdir('../../../..')
        sol = p.stdout.readline()
        all_sol = np.matrix(sol.decode())
        u_sol = np.array(all_sol[0, 0:10]).flatten()
        valid = np.array(all_sol[0, 10]).flatten()
        assert_allclose(u_sol, u_ref, rtol=1e-6,
                err_msg='Solve cond with pbm failed.')
        assert_allclose(valid, -1, rtol=1e-6,
                err_msg='Not valid.')

    @unittest.skip('This is how you skip')
    def test_pbm_solver_sparse_data(self):
        import muaompc
        muaompc_dir = os.path.abspath(os.path.dirname(muaompc.__file__))
        prefix = 'hhmpc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        data_name = 'regmpc'
        dat_path = os.path.join(base_dir, 'regmpc.dat')
        destdir = self.destdir
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=destdir,
                solver='pbm')
        generate_mpc_data(mpc, dat_path)
        shutil.os.chdir('%s/%s_%s/src' % (
            destdir, 'regrefmpc', 'hhmpc'))
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/static_data.h'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/hhmpcsocp.h'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/hhmpcsocp.c'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/socp_data.c'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/ipm_data.c'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/hhmpc_qpdatasmall10.h'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/pbm_dyn_data.json'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm_cond/mpcqpN10.json'), '.')
        shutil.os.chdir('../../..')
        shutil.os.chdir('%s/%s_%s/data/%s' % (
            destdir, 'regrefmpc', 'hhmpc', 'regmpc'))
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/main.c'), '.')
        shutil.os.remove(prefix+data_name+'main.c')
        call(['make', '-f', prefix+data_name+'Makefile.mk', 'clean'])
        call(['make', '-f', prefix+data_name+'Makefile.mk'])
        shutil.os.chdir('../../../..')
        shutil.os.chdir('%s/%s_%s/src' % (
            destdir, 'regrefmpc', 'hhmpc'))
        os.remove('ipm_data.c')
        os.remove('socp_data.c')
        os.remove('hhmpcsocp.c')
        os.remove('hhmpcsocp.h')
        shutil.os.chdir('../../..')
        shutil.os.chdir('%s/%s_%s/data/%s' % (
            destdir, 'regrefmpc', 'hhmpc', 'regmpc'))
        u_ref = np.array([-0.226552774757, 0.102873692334, -0.141520295300,
                          0.054100414765, -0.059100414932, 0.027108691716,
                          -0.022752888235, 0.013103179156, 0.006799832462,
                          0.013707830772])
        p = subprocess.Popen(['./main'],# + [str(arg) for arg in args],
                stdout=subprocess.PIPE)

        shutil.os.chdir('../../../..')
        sol = p.stdout.readline()
        all_sol = np.matrix(sol.decode())
        u_sol = np.array(all_sol[0, 0:10]).flatten()
        assert_allclose(u_sol, u_ref, rtol=1e-6,
                err_msg='Solve sparse with pbm failed.')

    @unittest.skip('This is how you skip')
    def test_common_math(self):
        import muaompc
        muaompc_dir = os.path.abspath(os.path.dirname(muaompc.__file__))
        prefix = 'hhmpc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        dat_path = os.path.join(base_dir, 'regmpc.dat')
        destdir = self.destdir
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=destdir)
        generate_mpc_data(mpc, dat_path)
        shutil.os.chdir('%s/%s_%s/src' % (
            destdir, 'regrefmpc', 'hhmpc'))
        shutil.copy(os.path.join(muaompc_dir + '/_ldt/codegen/common/test/math/main.c'), '.')
        shutil.copy(os.path.join(muaompc_dir + '/_ldt/codegen/common/test/math/static_data.h'), '.')
        shutil.copy(os.path.join(muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/hhmpcsocp.h'), '.')
        call(['make', '-f', prefix+'Makefile.mk', 'clean'])
        call(['make', '-f', prefix+'Makefile.mk'])
        shutil.os.remove('main.c')
        mtx = np.array([4., 2., 0., 2., 5., 2., 0., 2., 5.])
        mtx_dim = 3
        test_af = np.array([2., 0., 3., 1.])
        test_ab = np.array([4., 1., 0., 2.])
        test_b = np.array([4., 8., 4., 8.])
        # Test cholesky factorization
        args = list(mtx) + [mtx_dim] + (
                list(test_af) + list(test_ab) + list(test_b))
        chol_ref = np.array([2., 0., 0., 1., 2., 0., 0., 1., 2.])
        fwd_ref = np.array([2., 4., -2., -4.])
        bwd_ref = np.array([.5, 1., 2., 4.])
        p = subprocess.Popen(['./main'] + [str(arg) for arg in args],
                stdout=subprocess.PIPE)
        shutil.os.chdir('../../..')
        sol = p.stdout.readline()
        all_sol = np.matrix(sol.decode())
        chol_sol = np.array(all_sol[0, 0:9]).flatten()
        fwd_sol = np.array(all_sol[0, 9:13]).flatten()
        bwd_sol = np.array(all_sol[0, 13:17]).flatten()
        assert_allclose(chol_sol, chol_ref, rtol=1e-6,
                err_msg='Cholesky factorization failed.')
        assert_allclose(fwd_sol, fwd_ref, rtol=1e-6,
                err_msg='Forward substitution failed.')
        assert_allclose(bwd_sol, bwd_ref, rtol=1e-6,
                err_msg='Backward substitution failed.')

    @unittest.skip('This is how you skip')
    def test_common_solve(self):
        import muaompc
        muaompc_dir = os.path.abspath(os.path.dirname(muaompc.__file__))
        prefix = 'hhmpc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        dat_path = os.path.join(base_dir, 'regmpc.dat')
        destdir = self.destdir
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=destdir,
                solver='pbm')
        generate_mpc_data(mpc, dat_path)
        shutil.os.chdir('%s/%s_%s/src' % (
            destdir, 'regrefmpc', 'hhmpc'))
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/common/test/solve/main.c'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/common/test/solve/static_data.h'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm/hhmpcsocp.h'), '.')
        call(['make', '-f', prefix+'Makefile.mk', 'clean'])
        call(['make', '-f', prefix+'Makefile.mk'])
        shutil.os.remove('main.c')
        delta_z_ref = np.array([
            -13.997489150487263, 5.406972595368199, 4.776895838826961,
            10.733684585131803, 16.414840767556548, -156.578190818640479,
            -2.590516555119194, 3.871722379936190, -4.652813415677183,
            -9.071001815573380, -2.558390140372020, -48.856406163215773,
            16.464922653098494, 2.628392409139868, -0.939205567832684,
            -4.564250749299972, -0.856148616921228, -4.322564467882859,
            -1.350495556682932, 2.325126881373518, 9.124936092007225])
        delta_v_ref = ([
            629.359948524434230, -1066.284511494853632, 19.396311938666354,
            397.292038589655306, 5.268576650132512, 608.496770168496710,
            -731.217580993566003, 81.330459700416299, -6.689914187784114,
            -34.036261511762845, 3.590187702111360, 2.905130931383225,
            4.343104636331974, -11.941383543336487, -24.835093488191738])
        bz_ref = np.array([
            40.181558625426, 14.967635791964, -35.676453837091,
            93.756568935670, 24.153456576189, 403.981952777434,
            39.304838161894, 86.165267131823, 82.116526824200,
            28.660204330909, 31.319973202918, 109.377552081020,
            5.251469355554, -9.201168023570, 31.422291265929,
            25.938053761702, 3.590187702111, 2.905130931382,
            4.343104636333, -11.941383543336, -24.835093488191])
        bv_ref = np.array([
            40.425894872539, -23.196319452050, 20.169720013097, 7.756882323647,
            7.288940714645, -4.608134887490, 41.150471454975, 34.569259859659,
            3626.798702405204, 2.887658373030, -5.504725877673, 2.607627009342,
            -12.975952244762, 341.639130549762, -1.923378499420])
        p = subprocess.Popen(['./main'],  # + [str(arg) for arg in args],
                stdout=subprocess.PIPE)
        shutil.os.chdir('../../..')
        sol = p.stdout.readline()
        all_sol = np.matrix(sol.decode())
        delta_z_sol = np.array(all_sol[0, 0:21]).flatten()
        delta_v_sol = np.array(all_sol[0, 21:36]).flatten()
        bv_sol = np.array(all_sol[0, 36:51]).flatten()
        bz_sol = np.array(all_sol[0, 51:72]).flatten()
        # print(delta_v_sol)
        # print(delta_z_sol)
        # print(all_sol)
        assert_allclose(delta_z_sol, delta_z_ref, rtol=1e-6,
                err_msg='Calculation of delta_z failed.')
        assert_allclose(delta_v_sol, delta_v_ref, rtol=1e-6,
                err_msg='Calculation of delta_v failed.')
        assert_allclose(bv_sol, bv_ref, rtol=1e-6,
                err_msg='multiply_C_z failed.')
        assert_allclose(bz_sol, bz_ref, rtol=1e-6,
                err_msg='multiply_C_T_v failed.')

    @unittest.skip('This is how you skip')
    def test_socc_solver_static_data(self):
        H = 2*np.array([[ 1.00157741, 0.00133198],
            [0.00133198, 1.00146001]])
        g = 2*np.array([[0.33038257],
               [0.3054006]])
        V = np.array([[0.,  0.],
                    [0.10501746,  0.],
                    [0.01899904,  0.02095512]])
        Wm0 = np.array([[  9.67483607e-05,   0.00000000e+00],
                   [  1.90325164e-02,   0.00000000e+00]]).flatten()
        Wm1 = np.array([[  2.77866701e-04,   9.67483607e-05],
                   [  1.72213330e-02,   1.90325164e-02]]).flatten()
        wvT0 = np.array([[ -9.67483607e-05,   0.00000000e+00]]).flatten()
        wvT1 = np.array([[ -2.77866701e-04,  -9.67483607e-05]]).flatten()
        ws0 = np.array([[19.93649]]).flatten()  # np.array([[-2.37754877]]).flatten()
        ws1 = np.array([[20.19481]]).flatten()  # np.array([[-2.40338077]]).flatten()
        wn0 = -10.*np.array([[ 2.02854877],
                   [ 2.71451225]]).flatten()
        wn1 = -10.*np.array([[ 2.05438077],
                   [ 2.45619226]]).flatten()

        u_ub = np.array([[100], [100]])
        v_ub = np.array([[-22.5], [-21.02973811], [ -2.61268764]])

        import muaompc
        muaompc_dir = os.path.abspath(os.path.dirname(muaompc.__file__))

        prefix = 'stocmpc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'stocmpc.prb')
        data_name = 'stocmpc'
        dat_path = os.path.join(base_dir, data_name+'.dat')
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=self.destdir)
        generate_mpc_data(mpc, dat_path)
        shutil.os.chdir('%s/%s_%s/src' % (self.destdir, data_name, prefix))
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm_socc/static_data.h'), '.')
        shutil.copy(os.path.join(
            muaompc_dir + '/_ldt/codegen/solver/pbm/test/pbm_socc/stocmpcsocp.h'), '.')
        shutil.os.chdir('../../..')
        shutil.os.chdir('%s/%s_%s/data/%s' % (
            self.destdir, 'stocmpc', 'stocmpc', 'stocmpc'))
        #hutil.copy('/home/zmenende/workspace/muaompc/muaompc/_ldt/mpc/test/main.c', '.')
        shutil.copy(os.path.join(muaompc_dir + '/_ldt/mpc/test/mainstocmpc.c'), '.')
        shutil.os.remove(prefix+data_name+'main.c')
        call(['make', '-f', prefix+data_name+'Makefile.mk', 'clean'])
        call(['make', '-f', prefix+data_name+'Makefile.mk'])
        args = list(-10.*self.paramori['x_k']) + (
                [1000., -1000.])
        p = subprocess.Popen(['./main'] + [str(arg) for arg in args],
                stdout=subprocess.PIPE)
        shutil.os.chdir('../../../..')
        sol = p.stdout.readline()
        if isinstance(sol, str):  # probably Python2
            ctl_sol = np.matrix(sol)
        elif isinstance(sol, bytes):  # probably Python3
            ctl_sol = np.matrix(sol.decode())
        ctl_Wm0 = np.array(ctl_sol[0, 0:4]).flatten()
        ctl_wn0 = np.array(ctl_sol[0, 4:6]).flatten()
        ctl_wvT0 = np.array(ctl_sol[0, 6:8]).flatten()
        ctl_ws0 = np.array(ctl_sol[0, 8:9]).flatten()
        ctl_Wm1 = np.array(ctl_sol[0, 9:13]).flatten()
        ctl_wn1 = np.array(ctl_sol[0, 13:15]).flatten()
        ctl_wvT1 = np.array(ctl_sol[0, 15:17]).flatten()
        ctl_ws1 = np.array(ctl_sol[0, 17:18]).flatten()

        assert_allclose(ctl_Wm0, Wm0, rtol=1e-6,
                err_msg='Constraint matrix')
        assert_allclose(ctl_wvT0, wvT0, rtol=1e-6,
                err_msg='Constraint matrix')
        assert_allclose(ctl_wn0, wn0, rtol=1e-6,
                    err_msg='Constraint matrix')
        assert_allclose(ctl_ws0, ws0, rtol=1e-6,
                    err_msg='Constraint matrix')
        assert_allclose(ctl_Wm1, Wm1, rtol=1e-6,
                err_msg='Constraint matrix')
        assert_allclose(ctl_wvT1, wvT1, rtol=1e-6,
               err_msg='Constraint matrix')
        assert_allclose(ctl_wn1, wn1, rtol=1e-6,
                err_msg='Constraint matrix')
        assert_allclose(ctl_ws1, ws1, rtol=1e-6,
                err_msg='Constraint matrix')

    @unittest.skip('This is how you skip')
    def test_destdir_alm(self):
        prefix = 'mpcalm'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regorigmpc.prb')
        dat_path = os.path.join(base_dir, 'regmpc.dat')
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=self.destdir)
        generate_mpc_data(mpc, dat_path)
        #ys.path.append(os.path.join(self.destdir, 'regorigmpc'+'_'+prefix))
        #pc = import_module(prefix+'.mpc')
        #pc.generate_mpc_data(dat_path)

    # @unittest.skip('This is how you skip')
    def test_destdir_pbm(self):
        prefix = 'hhmpc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        dat_path = os.path.join(base_dir, 'regmpc.dat')
        destdir = self.destdir
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=destdir)
        generate_mpc_data(mpc, dat_path)

    @unittest.skip('This is how you skip')
    def test_destdir_fgm(self):
        prefix = 'mpcfgm'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        dat_path = os.path.join(base_dir, 'regmpc.dat')
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=self.destdir)
        generate_mpc_data(mpc, dat_path)

    @unittest.skip('This is how you skip')
    def test_matlab_codegen(self):
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regorigmpc.prb')
        setup_mpc_problem(prb_path, 'mpcmat')

    def test_datagen_microcontroller(self):
        prefix = 'muc'
        base_dir = os.path.dirname(regmpcdata.__file__)
        prb_path = os.path.join(base_dir, 'regrefmpc.prb')
        dat_path = os.path.join(base_dir, 'regmpc.dat')
        mpc = setup_mpc_problem(prb_path, prefix=prefix, destdir=self.destdir)
        generate_mpc_data(mpc, dat_path, muc=True)


    def _vecseq2array(self, vecseq):
        outer_rows = len(vecseq)
        n = vecseq[0].shape[0]
        arr = np.zeros((n*outer_rows, 1))

        for i, mtx in enumerate(vecseq):
                arr[n*i:n*(i+1), 0] = mtx

        return arr


if __name__ == '__main__':
    unittest.main()
