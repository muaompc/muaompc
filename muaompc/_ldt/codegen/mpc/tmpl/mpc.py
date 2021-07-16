import os
import pickle

from muaompc import ldt

def setup_mpc_ctl(num):
    mpc = _get_mpc()
    prefix = mpc.prefix
    __import__(prefix+'.'+prefix+'C'+mpc.former)
    s = __import__(prefix+'.'+prefix+mpc.solver)
    solver = s.__getattribute__(prefix+mpc.solver)
    data = mpc.ddg._get_data(num)
    ctl = solver.Solver()
    ctl.setup_solver(data)
    return ctl

def generate_mpc_data(fname):
    mpc = _get_mpc()
    ldt.generate_mpc_data(mpc, fname, safe_mode=False)

def _get_mpc():
    base_dir = os.path.dirname(__file__)
    prb_path = os.path.join(base_dir, 'mpc.pickle')
    with open(prb_path, 'rb') as f:
        mpc = pickle.load(f)
    return mpc
