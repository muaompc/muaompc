r"""
This module creates a model predictive control (MPC) controller for a
linear discrete-time (*ldt*) system with input and
(optionally) state constraints, and a convex cost function.
The MPC problem is reformulated as a condensed convex mathematical program,
which can be solved using off-the-shelf optimization algorithms,
or the ones included in muaompc.

"""

import os
import shutil
import pickle
import types

import numpy as np

from muaompc._ldt.parse import prbdsl
from muaompc._ldt.parse.prbsym import ProblemSym
from muaompc._ldt.parse.prbstruct import ProblemStruct
from muaompc._ldt.codegen.codegen import BaseDataGenerator
from muaompc._ldt.codegen.former.cvp.cvp import CondensedVectorialParameters as CVP
from muaompc._ldt.codegen.common import codegen as commoncodegen
from muaompc._ldt.codegen.common import codegenmat as commoncodegenmat
from muaompc._ldt.codegen.former.cvp import codegen as cvpcodegen
from muaompc._ldt.codegen.former.cvp import codegency as cvpcodegency
from muaompc._ldt.codegen.former.cvp import codegenmat as cvpcodegenmat
from muaompc._ldt.codegen.solver.pbm import codegen as pbmcodegen
from muaompc._ldt.codegen.solver.fgm import codegen as fgmcodegen
from muaompc._ldt.codegen.solver.fgm import codegency as fgmcodegency
from muaompc._ldt.codegen.solver.fgm import codegenmat as fgmcodegenmat
from muaompc._ldt.codegen.solver.alm import codegen as almcodegen
from muaompc._ldt.codegen.solver.alm import codegency as almcodegency
from muaompc._ldt.codegen.solver.alm import codegenmat as almcodegenmat
from muaompc._ldt.codegen.mpc import codegen as mpccodegen
from muaompc._ldt.codegen.mpc import codegency as mpccodegency
from muaompc._ldt.codegen.mpc import codegenmat as mpccodegenmat
from muaompc._ldt.mpc.mpc import MPC

ROWS, COLS = (0, 1)

class SetupError(Exception): pass


def setup_mpc_problem(fname, prefix='mpc', destdir='.', verbose=False, 
        numeric='float64', solver='alm+fgm'):
    """Form a problem description given an MPC problem high-level formulation.

    :param str fname: the name of the file with a valid MPC problem
    :param str prefix: the text to be prepended to the names of the generated
    files and to the internal code global identifiers,
    e.g. names of functions, structures,
    defines, enums, etc. (like a namespace)
    :param str destdir: the path to the directory where the generated code
    will be dumped. Default is the current working directory.
    :param bool verbose: if True, print out information about the code 
    generation process. 
    :param str numeric: specifies the numeric type for arithmetic operations.
    Supported values are:
    'float64' for 64-bit double precision floating point arithmetic,
    'float32' for 32-bit single precision floating point arithmetic.
    :param str solver: specifies the optimization algorithm to be used
    to solve the problem.
    Supported values are:
    'alm+fgm' a combination of augmented Lagrangian method (ALM) with fast gradient
    method (FGM). 
    If only input constraints are present, the FGM is used. Otherwise,
    ALM with the FGM as inner solver.
    """

    prs = prbdsl.parse_file(fname, verbose)
    sym = ProblemSym(prs)
    stt = ProblemStruct(sym)
    cvp = CVP(stt)
    former = 'cvp'
    prbname = _get_name(fname)
    ccg = commoncodegen.CCodeGenerator(prbname, prefix, destdir)
    real_t = _get_real_t(numeric)
    ccg.generate_code(real_t)
    matccg = commoncodegenmat.MatlabCodeGenerator(ccg)
    matccg.generate_code()
    cvpcg = cvpcodegen.CCodeGenerator(cvp, ccg)
    cvpcg.generate_code()
    cycdg = cvpcodegency.CythonCodeGenerator(cvpcg)
    cycdg.generate_code()
    matcdg = cvpcodegenmat.MatlabCodeGenerator(cvpcg)
    matcdg.generate_code()

    if (not cvp.is_socc_iec) and (solver == 'alm+fgm'):
        if cvp.is_inputbox_iec:
            cdg = fgmcodegen.CCodeGenerator(ccg)
            cycdg = fgmcodegency.CythonCodeGenerator(cdg)
            matcdg = fgmcodegenmat.MatlabCodeGenerator(cdg)
            ddg = fgmcodegen.FGMCVPDataGenerator(cvpcg)
            sdg = fgmcodegen.CDataGenerator(cdg, cvpcg)
            solver = 'fgm'
            mcg = mpccodegen.CCodeGenerator(ccg, solver=solver, former=former)
        else:
            cdg = almcodegen.CCodeGenerator(ccg)
            cycdg = almcodegency.CythonCodeGenerator(cdg)
            matcdg = almcodegenmat.MatlabCodeGenerator(cdg)
            ddg = almcodegen.ALMCVPDataGenerator(cvpcg)
            sdg = almcodegen.CDataGenerator(cdg, cvpcg)
            solver = 'alm'
            mcg = mpccodegen.CCodeGenerator(ccg, solver=solver, former=former)
    if cvp.is_socc_iec or (solver == 'pbm'):
        cdg = pbmcodegen.CCodeGenerator(ccg)  # Use pbmcodegen.
        cycdg = fgmcodegency.CythonCodeGenerator(cdg)
        matcdg = fgmcodegenmat.MatlabCodeGenerator(cdg)
        ddg = pbmcodegen.PBMCVPDataGenerator(cvpcg)
        sdg = pbmcodegen.CDataGenerator(cdg, cvpcg)
        solver = 'pbm'
        mcg = mpccodegen.CCodeGenerator(ccg, solver=solver, former=former)
    cdg.generate_code()
    cycdg.generate_code()
    matcdg.generate_code()
    mcg.generate_code()
    cymcg = mpccodegency.CythonCodeGenerator(mcg, cvpcg)
    cymcg.generate_code()
    matmcg = mpccodegenmat.MatlabCodeGenerator(mcg)
    matmcg.generate_code()
    matmcg.generate_matlab_make()  # this is always called last
    mpc = MPC(prbname, ccg, ddg, sdg, solver, former)
    base_dir = mpc.path['dest']
    prb_path = os.path.join(base_dir, 'mpc.pickle')
    with open(prb_path, 'wb') as f:
        pickle.dump(mpc, f)

    return mpc

def generate_mpc_data(mpc, fname, safe_mode=True, muc=False):
    """Generate data corresponding to the given MPC problem from a given a data file

    :param mpc: an MPC problem object created using the setup_mpc_problem function.
    :param str fname: the name of the file with a valid data file for 
    the given MPC problem object.
    :param bool safe_mode: only allow the simplest parsing method.
    This is useful for server applications, where there should be tighter control
    of what it is being parsed.
    :param bool muc: generate code for a microcontroller like Arduino.
    In general, this only means that the code for problem and data are put on a
    single directory.
    """

# fname can be an absolute path
    num = _get_data(mpc, fname, safe_mode)
    mcg = mpccodegen.CDataGenerator(mpc,
            solver=mpc.solver, former=mpc.former)
    dataname = _get_name(fname)
    mcg.generate_data(num, dataname)
    if not muc:
        mpc.ddg.generate_data(num, dataname)
    mpc.sdg.generate_data(num, dataname)
    if muc:
        _move_to_mucdir(mpc, dataname)
    return num

def find_penalty_parameters(prbfname, datfname, safe_mode=True, factors=[1], 
        verbose=False):
    mpc = setup_mpc_problem(prbfname)
    num = generate_mpc_data(mpc, datfname, safe_mode)
    if mpc.solver != 'alm':
        print('Solver %s does not have a penalty parameter' % mpc.solver)
        return

    mus, kappas = _calc_int_cn(mpc, num, factors, verbose)
    return dict(params=mus, condnums=kappas)

def _calc_int_cn(mpc, num, factors, verbose):
    H, V = mpc.ddg._get_matrices_HV(num)
    mu = mpc.ddg._find_penalty_parameter(H, V)
    penalties = [mu*factor for factor in factors]

    samples = len(penalties)
    mu_vec = np.zeros((samples))
    kappa_vec = np.zeros((samples))
    
    for k, mu in enumerate(penalties):
        kappa_vec[k] = mpc.ddg._compute_alg_condnum(H, mu, V) 
        mu_vec[k] = mu
    if verbose: 
        for k, mu in enumerate(penalties):
            print('For penalty:', mu_vec[k],
                  'the condition number of internal problem is:', 
                  kappa_vec[k])
        print('Strong convexity parameter:', 
              min(np.linalg.eigvals(H)))
    return (mu_vec, kappa_vec)

def _get_name(fname):
    if isinstance(fname, types.ModuleType):
        name_ext = fname.__name__
        name = name_ext.split('.')[-1]  # xxx.yyy.name
    else:
        name_ext = os.path.basename(fname)
        name = name_ext.split('.')[0]  # name.xxx.yyy -> name
    return name

def _get_data(mpc, fname, safe_mode):
    smb = mpc.sdg.smb
    cvp = mpc.sdg.cvp
#TODO: cvp might not always exist

    data = _get_data_from_file_name(fname, safe_mode)
    required_attr = smb.lst.fix

#TODO: there should be a way to know which attributes are expected to be
# matrices and which are scalars
    for attr_name in required_attr:
        attr = _get_attribute(data, attr_name)
        if attr is None:
            msg = ('Required attribute '+attr_name+' was not found.')
            raise SetupError(msg)

#TODO: there should be a list for reserved names, like mu
    if not cvp.is_inputbox_iec:  # state constraint, optionally look for mu
        mu = _get_attribute(data, 'mu')
        if (mu is not None):
            if not (isinstance(mu, float) or isinstance(mu, int)):
                print(type(mu))
                msg = ('Wrong type for penalty parameter mu. '
                        'Accepted types are float, int or None.')
                raise SetupError(msg)
            elif mu <= 0:
                msg = ('Penalty parameter must be greater than zero')
                raise SetupError(msg)
            else:
                data['mu'] = mu
    return data

def _get_attribute(d, key):
    try:
        attr = d[key]
    except(KeyError,):
        attr = None
    return attr

def _get_data_from_file_name(fname, safe_mode):
    if safe_mode:
        data = _parse_data_file(fname)
    else:
        data = _get_data_from_file(fname)

    return data

def _get_data_from_file(fname):
    if isinstance(fname, types.ModuleType):
        data = _get_data_dict_from_py(fname)
    elif isinstance(fname, str):
        data = _get_data_from_str(fname)
    else:
        msg = ('Wrong fname type. ' + str(fname) +
                ' is of type ' + str(type(fname)) +
                '. Check the documentation for accepted types.')
        raise SetupError(msg)
    return data

def _get_data_from_str(fname):
    name, dot, extension = fname.partition('.')
    if extension == 'mat':
        data = _get_data_from_mat(name)
    elif (extension == 'py') or (extension == ''):
        data = _get_data_from_py(name)
    else: 
        data = _parse_data_file(fname)
    
    return data

def _get_data_from_mat(name):
    try:
        from scipy import io
    except(ImportError,):
        msg = ('Attempting to open the data module from a '
              'MATLAB mat file but scipy.io could not be imported. '
              'Check that Scipy is correctly installed.')
        raise SetupError(msg)
    else:
        try:
            data = io.loadmat(name)
        except(IOError,):
            msg = ('Could not find MATLAB data module.')   
            raise SetupError(msg)

    return data

def _get_data_dict_from_py(module):
    try:
        data = module.data
    except(AttributeError,):
        msg = ('Could not find data dictionary in Python module.'
                'Check that data module contains a dictionary called "data"') 
        raise SetupError(msg)
    if not isinstance(data, dict):
        msg = ('Found "data" attribute in data module, but it is not a dictionary.'
                'Check that data module contains a dictionary called "data"') 
        raise SetupError(msg)
    return data

def _get_data_from_py(name):
    try:
        datamod = __import__(name)
    except(ImportError,):
        msg = ('Could not find Python data module.')   
        raise SetupError(msg)
        
    return _get_data_dict_from_py(datamod) 

def _parse_data_file(fname):
    with open(fname) as f:
        lines = f.readlines()

    d = dict()
    for linenum, line in enumerate(lines):
        try:
            (assign, comment) = line.split('#')
        except ValueError:  # comment not found
            assign = line
        try:
            (name_ws, data) = assign.split('=')
        except ValueError:  # assignment not found
            if len(assign.strip()) != 0:
                print("Warning: ignoring line", linenum+1)
                print(line)
        else:
            if data.strip() == '':
                msg = "Invalid assignment in line"+str(linenum+1)+"\n"
                msg += line
                raise SetupError(msg)
            name = name_ws.strip()
            try:
                d[name] = int(data)
            except ValueError:  # not an int
                pass
            else:
                continue

            try:
                d[name] = float(data)
            except ValueError:  # not a float
                pass
            else:
                continue

            values = data.strip('[]')  # safe input for np.matrix

            try: 
                npmtx = np.matrix(values)
            except ValueError as err:
                msg = '%s Offending matrix: %s' % (err, values)
                raise ValueError(msg)

            d[name] = np.array(npmtx)

    return d

def _move_to_mucdir(mpc, dataname):
    mucdirpath = _move_to_single_dir(mpc, dataname)
    mpc.path['muc'] = mucdirpath
    _remove_dynmem(mucdirpath)
    return

def _move_to_single_dir(mpc, dataname):
    singledirpath = os.path.join(mpc.destdir, 
                '%s_%s_%s' % (mpc.prbname, dataname, mpc.prefix))
    try:
        os.mkdir(singledirpath)
    except FileExistsError:
        pass
    _copy_files_to_dir('.c', mpc.path['dest'], singledirpath)
    _copy_files_to_dir('.h', os.path.join(mpc.path['dest'], 'include'), 
                            singledirpath)
    _copy_files_to_dir('.c', os.path.join(mpc.path['data'], dataname), 
                            singledirpath)
    _copy_files_to_dir('.h', os.path.join(mpc.path['data'], dataname), 
                            singledirpath)
    return singledirpath

def _copy_files_to_dir(ext, src, dest):
    for fname in os.listdir(src):
        if fname.endswith(ext):
            fpath = os.path.join(src, fname)
            shutil.copy(fpath, dest) 
    return

def _remove_dynmem(dirpath):
    for fname in os.listdir(dirpath):
        if (fname.endswith('dynmem.c') or fname.endswith('dynmem.h') or
            fname.endswith('main.c')):
            os.remove(os.path.join(dirpath, fname))
    return

def _get_real_t(numeric):
    if (numeric == 'float64') or (numeric == 'float32'):
        return numeric + '_t'
    else:
        msg = '"%s" is not a valid numeric value.' 
        raise SetupError(msg)

