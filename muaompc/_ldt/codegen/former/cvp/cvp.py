"""The class representing a Condensed problem with Vectorial Parameters."""

from copy import deepcopy

import numpy as np
from scipy.linalg import block_diag


def get_pmetric_none():
    return dict(fac0=np.zeros((0,0)))

class CondensedVectorialParameters(object):

    """Describe a condensed problem based on a structured problem description.

    The condensed problem is given by the cost function:
    1/2 (u.T * H * u) + (u.T * g(p_k))
    with u the optimization variable, and k=1:n_p
    and the double sided constraints:
    u_lb(p_k) <= u <= u_ub(p_k),
    v_lb(p_k) <= V * u <= v_ub(p_k)
    u.T * P_i * u <= a_i(p_k), i=1,...,n_q
    | Wm_j * u + wn_j(p_k) | <= wmT_j * u + ws_j(p_k), j=1,...,n_s

    Each parametric vector f(p_k) is described as follows:
    f(p_1, ..., p_{n_p}) = fac0 + sum(fac_k * p_k, k=1,...,n_p).
    P_i, is a matrix and a_i is a scalar.
    Wm_j is a matrix, wn_j is a column vector, wnT_j is a row vector,
    ws_j is a scalar.
    H, V and P_i, i=1,...,n_q, Wm_j, j=1,...,n_s, can be zero matrices.

    At the moment, only SOCPs without quadratic constraints are supported,
    i.e. P_i=0, for all i in {1, ..., n_q}.
    """

    def __init__(self, stc):

        """Identify parametric terms in given structured problem.

        The terms depending on zero or more parameters are assigned to
        the dict self.pmetric. Similarly, its keys are the names
        of the parametric vectors (e.g. 'g').
        If a terms depends on zero parameters, it is considered constant.
        The items in self.pmetric, however, are each a dictionary with the
        keys p[k].name identifying the factor that multiplies parameter k.
        The names of the parameters are given in the MPC problem definition.
        For example, if 'g' depends on the parameter 'x_k', then parameter
        information is accessed by self.pmetric['g']['x_k'].

        pmetric (and constant) terms depend on the structure of the problem.
        Problems with only input box constraint have different (less) keys
        in both dictionaries compared to more general problems.

        This function also assign to self.par a dict with the parameter names
        as keys, and the contents of smb.par as values,
        where smb is an instance of muaompc._ldt.parse.prbsym.SymbolElements.

        :param stc: an instance of the class representing the problem structure
        :type stc: muaompc._ldt.parse.prbstruct.ProblemStruct

        """

        self.stc = stc
        self.smb = stc.smb
        self.par = self._par_list2dict(stc.smb.par)
        self.min_ = stc.min_
        self.eqcs = stc.eqcs
        self.iecs = stc.iecs
        self.soccs = stc.soccs
        self.is_inputbox_iec = self._get_is_inputbox_iec()
        self.is_socc_iec = self._get_is_socc_iec()
        self.is_sparse_data = self._get_is_sparse_data()
        self.pmetric = self._get_pmetric()

    def form_condensed_prb(self, num):
        """Return numeric values of the pmetric terms.

        It replaces the equality aux variable found in eqcs into the
        cost function min_ and the inequality constraints iecs.
        These are assigned to the corresponding keys in
        the returned dict pmetric.

        In contrast to self.pmetric, the returned pmetric dict
        has the additional key 'fac0' for each entry, e.g.
        pmetric['g']['fac0'].

        :param num: contains numeric values of fix symbols
        :param type: dict

        :return: a dict with the pmetric terms of a QP, 
        description of SOCCs, and
        the length of Lagrange multipliers.
        """
        eqcsnum = self._propagate_multistage(num)
        pmetric = deepcopy(self.pmetric)
        pmetric = self._set_costfunc(pmetric, eqcsnum, num)
        pmetric, lagmul_len = self._set_iecs(pmetric, eqcsnum, num)
        socc = self._set_soccs(eqcsnum, num)
        return dict(pmetric=pmetric, lagmul_len=lagmul_len,
                socc=socc)

    def _get_flags(self):
        flags = dict()
        flags['inputbox'] = self.is_inputbox_iec
        flags['affine'] = not self.is_inputbox_iec
        flags['socc'] = self.is_socc_iec
        flags['sparse'] = self.is_sparse_data
        for key in flags:
            flags[key] = int(flags[key])
        return flags

    def _par_list2dict(self, par):
        pard = dict()
        for par_i in par:
            pard[par_i.name] = par_i
        return pard

    def _get_pmetric(self):
        g = get_pmetric_none()
        H = get_pmetric_none()
        u_lb = get_pmetric_none()
        u_ub = get_pmetric_none()
        V = get_pmetric_none()
        v_lb = get_pmetric_none()
        v_ub = get_pmetric_none()
        return dict(H=H, g=g, 
                u_lb=u_lb, u_ub=u_ub, 
                V=V, v_lb=v_lb, v_ub=v_ub)

    def _get_prb_terms(self):
        prb_terms = ['H', 'g', 'u_lb', 'u_ub']
        if not self.is_inputbox_iec:
            prb_terms += ['V', 'v_lb', 'v_ub']
        return prb_terms

    def _propagate_multistage(self, num):
        structnum = []
        eqcsnum = self.eqcs.sym2num(num)
        for k, struct in enumerate(self.eqcs.struct):
            structnum.append(eqcsnum[k])
            rngend = 1 + eval(struct['rng'].end, num)
            rows = struct['inifunc']._get_rows(num)
            for i in range(rngend):
                structnum[k] = self._mul_fac_struct(structnum[k], i, rows)
        return structnum

    def _mul_fac_struct(self, structnum, i, rows):
        ini = (i+1)*rows
        end = (i+2)*rows
        for j, aux in enumerate(structnum['auxs']):
            aux_ = aux[ini:end,:]
            for key, mtxfacs in structnum.items():
                if type(mtxfacs) is list:
                    for k, mtxfac in enumerate(mtxfacs):
                        if (key == 'auxs' and k==j):  # avoid recursive MAC operation on aux
                            continue
                        structnum[key][k][ini:end,:] += np.dot(aux_, mtxfac)
                else:
                    structnum[key][ini:end,:] += np.dot(aux_, mtxfacs)
        return structnum

    def _replace_eqcs(self, eqcsnum, funcnum):
        for j, aux in enumerate(funcnum['auxs']):
            if not np.allclose(aux, np.zeros(aux.shape)):
                for key, mtxfacs in funcnum.items():
                    if (key == 'quadmtx'):
                        continue  # quadmtx handled in a different function
                    if type(mtxfacs) is list:
                        for k, mtxfac in enumerate(mtxfacs):
                            if (key == 'auxs' and k==j):
                                continue  # avoid iterative aux[j] with itself
                            funcnum[key][k] += np.dot(aux, eqcsnum[j][key][k])
                    else:
                        funcnum[key] += np.dot(aux, eqcsnum[j][key])
        return funcnum

    def _set_costfunc(self, pmetric, eqcsnum, num):
        structnum = self.min_.sym2num(num)
        varlen = eval(self.smb.opt.seqlen, num)
        H_const = np.zeros((varlen, varlen))
        H_param = [np.zeros((varlen, parlen)) for parlen in self._get_parlen(num)]
        g_const = np.zeros((varlen, 1))
        g_param = [np.zeros((varlen, parlen)) for parlen in self._get_parlen(num)]
        for k, funcnum in enumerate(structnum):
            rows, cols = funcnum['quadmtx'].shape
            vecinseq = rows // cols
            quadmtx = block_diag(*([funcnum['quadmtx'][cols*i:cols*(i+1), :]
                                    for i in range(vecinseq)]))
            funcnum = self._replace_eqcs(eqcsnum, funcnum)
            varmtx = funcnum['variable']
            varTquadmtx = np.dot(varmtx.T, quadmtx)
            H_const += np.dot(varTquadmtx, varmtx)
            g_const += np.dot(varTquadmtx, funcnum['fix'])
            for k, param in enumerate(funcnum['parameters']):
                g_param[k] += np.dot(varTquadmtx, param)

        H_const *= 2.  # we want 1/2 u.T * H_const * u
        g_const *= 2.  # this is due to the symmetric non-quadratic terms
        for k, dummy in enumerate(g_param):
            g_param[k] *= 2.
            pmetric['g'][self.smb.par[k].name] = g_param[k]
        pmetric['g']['fac0'] = g_const
        pmetric['H']['fac0'] = H_const  # H assumed constant
        return pmetric

    def _get_parlen(self, num):
        par_len = []
        for par in self.smb.par:
            seqlen = par.eval_seqlen(num)
            par_len.append(seqlen)

        return par_len

    def _set_iecs(self, pmetric, eqcsnum, num):
        varlen = eval(self.smb.opt.seqlen, num)
        V_const = np.zeros((0, varlen))
        v_lb = np.zeros((0, 1))
        v_ub = np.zeros((0, 1))
        v_par = [np.zeros((0, parlen)) for parlen in self._get_parlen(num)]
        structnum = self.iecs.sym2num(num)
        lagmul_veclens = []
        for k, iec in enumerate(structnum):
            funcnum = iec['centres']
            is_inputbox = self._check_inputbox_iec(
                self.iecs.struct['centres'][k].struct)
            if is_inputbox:
                pmetric['u_lb']['fac0'] = iec['lowers']['fix']
                pmetric['u_ub']['fac0'] = iec['uppers']['fix']
            else:
                funcnum = self._replace_eqcs(eqcsnum, funcnum)
                V_const = np.vstack((V_const, funcnum['variable']))
                for k, param in enumerate(funcnum['parameters']):
                    v_par[k] = np.vstack((v_par[k], param))
                v_lb = np.vstack((v_lb, iec['lowers']['fix']))
                v_ub = np.vstack((v_ub, iec['uppers']['fix']))
                lagmul_veclens.append(iec['veclen'])

        if not self.is_inputbox_iec:
            for k, dummy in enumerate(v_par):
                pmetric['v_lb'][self.smb.par[k].name] = -v_par[k]
                pmetric['v_ub'][self.smb.par[k].name] = -v_par[k]
            pmetric['v_lb']['fac0'] = v_lb
            pmetric['v_ub']['fac0'] = v_ub
            pmetric['V']['fac0'] = V_const
            if self._is_lagmul_veclen_consistent(lagmul_veclens):
                lagmul_veclen = lagmul_veclens[0]  # all veclens equal to first
                lagmul_seqlen = V_const.shape[0]  # number of rows
                lagmul_len = dict(veclen=lagmul_veclen, seqlen=lagmul_seqlen)
            else:
                lagmul_len = None
                msg = "WARNING: the mixed inequality constraint sequence "
                msg += "do not consist of vectors of equal lengths. "
                msg += "Warm start of Lagrange multipliers "
                msg += "will not be possible."
                print(msg)
        else:
            lagmul_len = None

        return (pmetric, lagmul_len)

    def _set_soccs(self, eqcsnum, num):
        structnum = self.soccs.sym2num(num)
        socc_pmetric = []
        for k, socc in enumerate(structnum):
            pmetric = dict()
            pmetric['ws'] = get_pmetric_none()
            pmetric['wn'] = get_pmetric_none()
            pmetric['Wm'] = get_pmetric_none()
            pmetric['wvT'] = get_pmetric_none()
            funcnum = socc['centres']
            funcnum = self._replace_eqcs(eqcsnum, funcnum)
            pmetric['Wm']['fac0'] = funcnum['variable']
            for j, param in enumerate(funcnum['parameters']):
                pmetric['wn'][self.smb.par[j].name] = param
            pmetric['wn']['fac0'] = funcnum['fix']

            funcnum = socc['uppers']
            funcnum = self._replace_eqcs(eqcsnum, funcnum)
            pmetric['wvT']['fac0'] = funcnum['variable']
            for j, param in enumerate(funcnum['parameters']):
                pmetric['ws'][self.smb.par[j].name] = param
            pmetric['ws']['fac0'] = funcnum['fix']

            socc_pmetric.append(pmetric)

        socc = dict(pmetric=socc_pmetric)

        return socc

    def _get_is_inputbox_iec(self):
        # return True if ALL inequality constraint are input box
        for aff in self.iecs.struct['centres']:
            is_inputbox_iec = self._check_inputbox_iec(aff.struct)
            if not is_inputbox_iec:  # one non input box constraint found
                break
        return is_inputbox_iec

    def _get_is_socc_iec(self):
        is_socc_iec = False
        if len(self.soccs.struct['uppers']) != 0:
            is_socc_iec = True
        return is_socc_iec

    def _get_is_sparse_data(self):
        is_sparse_data = False
        if self.is_inputbox_iec:  # TODO Condition for use of sparse data.
            is_sparse_data = False
        return is_sparse_data


    def _check_inputbox_iec(self, struct):
        conditions = []
        for key in struct.keys():
            if key == 'variable':
                if (not struct[key].is_zero) and (struct[key].mtxs == [None]):
                    conditions.append(True)
                else:
                    conditions.append(False)
                    break
            else:
                mtxfacs = struct[key]
                if isinstance(mtxfacs, list):
                    for mtxfac in mtxfacs:
                        if mtxfac.is_zero:
                            conditions.append(True)
                        else:
                            conditions.append(False)
                            break
                else:
                    if mtxfacs.is_zero:
                        conditions.append(True)
                    else:
                        conditions.append(False)
                        break
        return np.alltrue(conditions)

    def _is_lagmul_veclen_consistent(self, lagmul_veclens):
        return lagmul_veclens[1:] == lagmul_veclens[:-1]
