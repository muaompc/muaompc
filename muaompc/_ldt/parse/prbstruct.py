"""Transform the parsed MPC problem into a structured SOCP representation.

"""
from collections import namedtuple
from copy import deepcopy

import numpy as np
import muaompc._ldt.parse.prbsym
desc_fields = ['veclen', 'seqlen', 'is_zero', 'rngend', 'rngidx', 'idxs', 'signs', 'mtxs']
MatrixFactor = namedtuple('MatrixFactor', desc_fields)

class ProblemStruct(object):
    def __init__(self, sym):
        self.smb = sym.smb
        self.eqcs = EqConstrFullAffineVectors(sym.smb, sym.sub.eqcs)
        self.min_ = QuadCostFuncFullVectors(sym.smb, sym.min_)
        self.iecs = IneqConstrFullAffineVectors(sym.smb, sym.sub.iecs)
        self.soccs = IneqConstrSecondOrderCone(sym.smb, sym.sub.soccs)

class AffineFunction(object):
    """The elementary function v[i] = c + D*u + sum(Pk*pk, k=1:n_p)
                                              + sum(Ek*ek, k=1:n_e)
    """

    def __init__(self, smb, elems):
        """
        The function v[i] is described by the given elems, and is transformed
        in the dict self.struct, with the keys:
        variable: MatrixFactor of optimization variable (D*u)
        auxs: a list of MatrixFactor of n_e auxiliary terms
        defined in equality constraint (Ek*ek)
        fix: MatrixFactor of a fix/constant term (c)
        parameters: a list of MatrixFactor of n_p parameters (Pk*pk)

        Additionally, self.row_mtx contains the symbol for a matrix that
        defines the number of rows in the vector v[i]

        :param smb: the symbols in the optimization problem
        :type smb: prbsym.SymbolElements
        :param elems: the definition of the function
        :type elems: a list of prbsym.Element
        """
        self.smb = smb
        self.struct = self._parse_expr2struct(elems)
        self.rows_mtx = self._get_row_mtx()

    def sym2num(self, num):
        structnum = deepcopy(self.struct)
        rows = self._get_rows(num)
        for key, mtxfacs in self.struct.items():
            if type(mtxfacs) is list:
                mfl = []
                for mtxfac in mtxfacs:
                    mtxfacnum = self._mtxfac_sym2num(mtxfac, num)
                    mfl.append(self._assign(mtxfacnum, rows))
                structnum[key] = mfl
            else:
                mtxfacnum = self._mtxfac_sym2num(mtxfacs, num)
                structnum[key] = self._assign(mtxfacnum, rows)
        return structnum

    def _get_rows(self, num):
        rows_mtx = eval(self.rows_mtx, num)
        if type(rows_mtx) is np.ndarray:
            rows = rows_mtx.shape[0]
        else:
            rows = rows_mtx
        return rows

    def _mtxfac_sym2num(self, mtxfac, num):
        mtxfacnum = []
        for field in mtxfac:
            if type(field) is list:
                numlist = []
                for sym in field:
                    numlist.append(eval(str(sym), num))
                mtxfacnum.append(numlist)
            else:
                mtxfacnum.append(eval(str(field), num))

        return MatrixFactor._make(mtxfacnum)

    def _parse_expr2struct(self, elems):
        addops = '+-'
        mulops = '*'
        d = dict.fromkeys(desc_fields)
        d['is_zero'] = True
        d['idxs'] = []
        d['signs'] = []
        d['mtxs'] = []
        factors = dict(variable=deepcopy(d),
                parameters=[deepcopy(d) for par in self.smb.par],
                auxs=[deepcopy(d) for aux in self.smb.aux],
                fix=deepcopy(d))
        for i, elem in enumerate(elems):
            type_ = elem.type_
            name = elem.name
            try:
                elem_postop = elems[i+1].preop
            except IndexError:
                elem_postop = '+'  # default op for no further elements

            if elem.preop in addops:
                if elem_postop in addops:
                    self._parse_struct_addadd(factors, type_, elems, i)
                elif elem_postop in mulops:
                    self._parse_struct_addmul(factors, type_, elems, i)
            elif elem.preop in mulops:
                if elem_postop in addops:
                    self._parse_struct_muladd(factors, type_, elems, i)
                elif elem_postop in mulops:
                    self._parse_struct_muladd(factors, type_, elems, i)

        factors = self._set_len(factors)
        factors['variable'] = MatrixFactor(**factors['variable'])
        factors['fix'] = MatrixFactor(**factors['fix'])
        for k, E in enumerate(factors['auxs']):
            factors['auxs'][k] =  MatrixFactor(**E)
        for k, P in enumerate(factors['parameters']):
            factors['parameters'][k] =  MatrixFactor(**P)

        return factors

    def _set_len(self, factors):
        # fix gets dummy values. The right values are to be infered later on
        factors['fix']['veclen'] = 1
        factors['fix']['seqlen'] = None
        factors['variable']['veclen'] = self.smb.opt.veclen
        factors['variable']['seqlen'] = self.smb.opt.seqlen
        for k, E in enumerate(factors['auxs']):
            factors['auxs'][k]['veclen'] = self.smb.aux[k].veclen
            factors['auxs'][k]['seqlen'] = self.smb.aux[k].seqlen
        for k, P in enumerate(factors['parameters']):
            factors['parameters'][k]['veclen'] = self.smb.par[k].veclen
            factors['parameters'][k]['seqlen'] = self.smb.par[k].seqlen

        return factors

    def _parse_struct_addadd(self, factors, type_, elems, i):
        # +e+ any type is okay
        def _set_values(dic, elems, i):
            dic['is_zero'] = False
            dic['signs'].append(elems[i].preop+'1')
            dic['mtxs'].append(None)
            dic['idxs'].append(elems[i].idx)

        if type_ == 'fix':
            _set_values(factors[type_], elems, i)
            factors[type_]['mtxs'] = [elems[i].name]
        elif type_ == 'variable':
            _set_values(factors[type_], elems, i)
        elif (type_ == 'auxs'):
            for k, aux in enumerate(self.smb.aux):
                if aux.name == elems[i].name:
                    _set_values(factors[type_][k], elems, i)
                else:
                    pass
        elif (type_ == 'parameters'):
            for k, par in enumerate(self.smb.par):
                if par.name == elems[i].name:
                    _set_values(factors[type_][k], elems, i)
                else:
                    pass

    def _parse_struct_addmul(self, factors, type_, elems, i):
        # +e* only fix matrices
        if type_ == 'fix':
            pass  # do nothing
        else:
            msg = 'Error in parsing.'
            raise TypeError(msg)  # FIXME: raise your own error

    def _parse_struct_muladd(self, factors, type_, elems, i):
        # *e+ no fixed matrices
        def _set_values(dic, elems, i):
            dic['is_zero'] = False
            dic['idxs'].append(elems[i].idx)
            dic['mtxs'].append(elems[i-1].name)
            dic['signs'].append(elems[i-1].preop+'1')  # +/- 1

        if type_ == 'fix':
            msg = 'Error in parsing.'
            raise TypeError(msg)  # FIXME: raise your own error
        elif (type_ == 'variable'):
            _set_values(factors[type_], elems, i)
        elif (type_ == 'auxs'):
            for k, aux in enumerate(self.smb.aux):
                if aux.name == elems[i].name:
                    _set_values(factors[type_][k], elems, i)
                else:
                    pass
        elif (type_ == 'parameters'):
            for k, par in enumerate(self.smb.par):
                if par.name == elems[i].name:
                    _set_values(factors[type_][k], elems, i)
                else:
                    pass

    def _parse_struct_mulmul(self, factors, type_, elems, i):
        # *e* not supported
        msg = 'Error in parsing.'
        raise TypeError(msg)  # FIXME: raise your own error

    def _get_row_mtx(self):
        for key, mtxfacs in self.struct.items():
            if key == 'variable':
                if not mtxfacs.is_zero:
                    if mtxfacs.mtxs[0] is None:
                        return self.smb.opt.veclen
                    else:
                        return mtxfacs.mtxs[0]
            elif key == 'auxs':
                for k, mtxfac in enumerate(mtxfacs):
                    if not mtxfac.is_zero:
                        if mtxfac.mtxs[0] is None:
                            return self.smb.aux[k].veclen
                        else:
                            return mtxfac.mtxs[0]
            elif key == 'parameters':
                for k, mtxfac in enumerate(mtxfacs):
                    if not mtxfac.is_zero:
                        if mtxfac.mtxs[0] is None:
                            return self.smb.par[k].veclen
                        else:
                            return mtxfac.mtxs[0]

    def _assign(self, mtxfac, rows):
        if mtxfac.is_zero:
            if mtxfac.seqlen is None:
                mtx = np.zeros((rows, mtxfac.veclen))
            else:
                mtx = np.zeros((rows, mtxfac.seqlen))
        else:
            mtx = self._assign_nonzero(mtxfac, rows)

        return mtx

    def _assign_nonzero(self, mtxfac, rows):
# consider four cases for the binary combination (idx, mtx)
        if mtxfac.seqlen is None:
            mtx = self._assign_mtx(mtxfac, rows, 0)
        else:
            mtx = np.zeros((rows, mtxfac.seqlen))
            for k, idx in enumerate(mtxfac.idxs):
                ini = idx*mtxfac.veclen
                end = (idx+1)*mtxfac.veclen
                mtx[:, ini:end] = self._assign_mtx(mtxfac, rows, k)
        return mtx

    def _assign_mtx(self, mtxfac, rows, k):
        Mk = mtxfac.mtxs[k]
        if Mk is None:
            if mtxfac.veclen != rows:
                msg = "Assigment of vectors of different shape."
                raise TypeError(msg)
            mtx = mtxfac.signs[k] * np.eye(rows)
        else:
            # FIXME: check shape compatibility of M
            mtx = mtxfac.signs[k] * Mk

        return mtx

class FunctionRanged(AffineFunction, object):
    """Extend base class with a range.

    Extend base class by
    providing basic functionality for vertically stacking a
    base class function for a given range.
    This class is meant to be subclassed.
    """
    def __init__(self, smb, elem, rng):
        pass

    def _vstack_struct_ini(self, inifuncnum, rngend):
        structnum = dict()
        for key, mtxfacs in inifuncnum.items():
            if type(mtxfacs) is list:
                fullmtx = []
                for mtxfac in mtxfacs:
                    rows, cols = mtxfac.shape
                    fullmtx.append(np.zeros((rows*(rngend), cols)))
                structnum[key] = fullmtx
            else:
                rows, cols = mtxfacs.shape
                fullmtx = np.zeros((rows*(rngend), cols))
                structnum[key] = fullmtx
        return structnum

    def _vstack_struct(self, structnum, funcnum, i, rows):
        for key, mtxfacs in funcnum.items():
            if type(mtxfacs) is list:
                for k, mtxfac in enumerate(mtxfacs):
                    structnum[key][k][i*rows:(i+1)*rows,:] = mtxfac
            else:
                structnum[key][i*rows:(i+1)*rows,:] = mtxfacs


class QuadFunction(AffineFunction, object):
    """A quadratic form v[i].T * Q * v[i]
    """
    def __init__(self, smb, elem):
        """
        The function v[i] is described by the given elem.vec, and is passed
        to the constructor of the parent class.
        Extend the parent class with self.mtx, which is the name of the
        symbol of the matrix (Q)

        :param smb: the symbols in the optimization problem
        :type smb: prbsym.SymbolElements
        :param elem: the definition of the quadratic form
        :type elem: prbsym.Qua
        """
        AffineFunction.__init__(self, smb, elem.vec.elems)
        self.mtx = elem.mtx.elems[0].name
        self.sign = elem.vec.preop + '1'

    def sym2num(self, num):
        structnum = AffineFunction.sym2num(self, num)
        sign = eval(self.sign, num)
        mtx = sign * eval(self.mtx, num)
        structnum['quadmtx'] = mtx
        return structnum


class EqConstrFullAffineVectors(FunctionRanged, object):
    """v = c + D*u + sum(Pk*pk, k=1:n_p) + sum(Ek*ek, k=1:n_e)
    """

    def __init__(self, smb, eqcs):
        """
        Similar to the base class, but self.structs is list of
        dicts (one dict for each aux in self.smb.auxs)
        with the following keys:
        inifunc: function independent of its own auxiliary variable.
        This describe the function at one end of the stack.
        rngfunc: a general function that is to be stacked iteratively
        after inifunc for the given range.
        rng: a prbsym.RangedIndex instance

        :param smb: the symbols in the optimization problem
        :type smb: prbsym.SymbolElements
        :param eqcs: the definition of the equality constraints
        :type eqcs: list of prbsym.EqConstr
        """
        self.smb = smb
        self.struct = self._parse_expr2struct(eqcs)

    def sym2num(self, num):
        eqcsnum = []
        for struct in self.struct:
            rngend = 1 + eval(struct['rng'].end, num)

            inifuncnum = struct['inifunc'].sym2num(num)
            rows = struct['inifunc']._get_rows(num)
            structnum = self._vstack_struct_ini(inifuncnum, rngend+1)

            self._vstack_struct(structnum, inifuncnum, 0, rows)
            for i in range(rngend):
                num[struct['rng'].var] = i
                rngfuncnum = struct['rngfunc'].sym2num(num)
                self._vstack_struct(structnum, rngfuncnum, i+1, rows)
            eqcsnum.append(structnum)
        return eqcsnum

    def _parse_expr2struct(self, eqcs):
        struct = []
        for aux in self.smb.aux:
            inifunc = self._get_aux_ini(aux, eqcs)
            rngfunc, rng = self._get_aux_rng(aux, eqcs)
            struct.append(dict(inifunc=inifunc, rngfunc=rngfunc, rng=rng))
        return struct

    def _get_aux_ini(self, aux, eqcs):
        eqc_aux_func = None
        for eqc in eqcs:
            eqc_aux = eqc.left.elems[0]
            if ((eqc_aux.name == aux.name) and
                    (eqc_aux.idx == '0') and
                    (eqc.rng is None)):
                eqc_aux_func = AffineFunction(self.smb, eqc.right.elems)

        if eqc_aux_func is None:
            msg = 'Not found initialization value for ' + aux.name
            raise TypeError(msg)

        return eqc_aux_func

    def _get_aux_rng(self, aux, eqcs):
        eqc_aux_func = None
        for eqc in eqcs:
            eqc_aux = eqc.left.elems[0]
            if ((eqc_aux.name == aux.name) and
                    (eqc.rng is not None)):
                eqc_aux_func = AffineFunction(self.smb, eqc.right.elems)
                eqc_aux_rng = eqc.rng
                idxvar = eqc.rng.var
                eqc_aux_stagestep = eval(eqc_aux.idx, {idxvar:0})

        if eqc_aux_func is None:
            msg = 'Not range initialization values for ' + aux.name
            raise TypeError(msg)
        if eqc_aux_stagestep != 1:
            msg = 'Stage step must be equal 1 for' + aux.name
            raise TypeError(msg)

        return eqc_aux_func, eqc_aux_rng


class IneqConstrFullAffineVectors(FunctionRanged, object):
    """c_lb <= D*u + sum(Pk*pk, k=1:n_p) + E*e <= c_ub
    """

    def __init__(self, smb, iecs):
        """
        Similar to the base class, but self.struct is a dict
        with keys:
        lowers: list of AffineFunction with only nonzero constant vector (c_lb)
        uppers: list of AffineFunction of only nonzero constant vector (c_ub)
        centres: list of AffineFunction with zero constant vector
        rngs: list of prbsym.RangedIndex instance

        :param smb: the symbols in the optimization problem
        :type smb: prbsym.SymbolElements
        :param iecs: the definition of the inequality constraints
        :type iecs: list of prbsym.IneqConstr
        """
        self.smb = smb
        self.struct = self._parse_expr2struct(iecs)

    def sym2num(self, num):
        iecs = []
        for k, rng in enumerate(self.struct['rngs']):
            iecnum = dict()
            for key in self.struct.keys():
                if key == 'rngs':
                    continue
                func = self.struct[key][k]
                if rng is not None:
                    rngend = 1 + eval(rng.end, num)
                    num[rng.var] = 0  # dummy value
                    funcnum = func.sym2num(num)
                    rows = func._get_rows(num)
                    structnum = self._vstack_struct_ini(funcnum, rngend)
                    for i in range(rngend):
                        num[rng.var] = i
                        funcnum = func.sym2num(num)
                        self._vstack_struct(structnum, funcnum, i, rows)
                else:
                    structnum = func.sym2num(num)
                    rows = func._get_rows(num)
                iecnum[key] = structnum
                iecnum['veclen'] = rows
            iecs.append(iecnum)
        return iecs

    def _parse_expr2struct(self, iecs):
        lowers = []
        centres = []
        uppers = []
        rngs = []
        for iec in iecs:
            lower = AffineFunction(self.smb, iec.lower.elems)
            centre = AffineFunction(self.smb, iec.centre.elems)
            upper = AffineFunction(self.smb, iec.upper.elems)
            lower.rows_mtx = centre.rows_mtx
            upper.rows_mtx = centre.rows_mtx
            # Only single vectors supported as fixed values
            lower.struct['fix'] = lower.struct['fix']._replace(veclen=lower.rows_mtx, seqlen=None)
            upper.struct['fix'] = upper.struct['fix']._replace(veclen=upper.rows_mtx, seqlen=None)

            rngs.append(iec.rng)
            lowers.append(lower)
            centres.append(centre)
            uppers.append(upper)

        return dict(lowers=lowers, centres=centres, uppers=uppers, rngs=rngs)


class QuadCostFuncFullVectors(FunctionRanged, object):
    """Purely quadratic cost functions, no linear terms
    """

    def __init__(self, smb, min_):
        """
        Similar to the base class, but self.struct is a dict
        with keys:
        quads: a list of QuadFunction instances
        rngs: a list of prbsym.RangedIndex instances

        :param smb: the symbols in the optimization problem
        :type smb: prbsym.SymbolElements
        :param min_: the definition of the cost function
        :type min_: prbsym.Min
        """
        self.smb = smb
        self.min_ = min_
        self.struct = self._parse_expr2struct()

    def sym2num(self, num):
        quads = []
        for k, func in enumerate(self.struct['quads']):
            rng = self.struct['rngs'][k]
            if rng is not None:
                rngend = 1 + eval(rng.end, num)
                num[rng.var] = 0  # dummy value
                funcnum = func.sym2num(num)
                rows = func._get_rows(num)
                structnum = self._vstack_struct_ini(funcnum, rngend)
                for i in range(rngend):
                    num[rng.var] = i
                    funcnum = func.sym2num(num)
                    self._vstack_struct(structnum, funcnum, i, rows)
            else:
                structnum = func.sym2num(num)
            quads.append(structnum)
        return quads

    def _parse_expr2struct(self):
        quads = []
        rngs = []
        for term in self.min_.elems:
            try:
                quads.append(QuadFunction(self.smb, term))
                rngs.append(None)
            except AttributeError:  # not a quad
                sum_ = term
                rng = sum_.rng
                for term_ in sum_.elems:
                    try:
                        quads.append(QuadFunction(self.smb, term_))
                        rngs.append(rng)
                    except AttributeError:
                        msg = "Only quadratic terms are supported"
                        raise TypeError(msg)

        return dict(quads=quads, rngs=rngs)

class IneqConstrSecondOrderCone(AffineFunction, object):
    """||v_[i]|| <= v[i]
    where v_[i] and v[i] have the form defined in AffineFunction
    and || . || is the Euclidean norm
    """

    def __init__(self, smb, soccs):
        """
        Similar to the base class, but self.struct is a dict
        with keys:
        centres: list of AffineFunction
        uppers: list of AffineFunction
        rngs: list of prbsym.RangedIndex instance

        :param smb: the symbols in the optimization problem
        :type smb: prbsym.SymbolElements
        :param soccs: the definition of the inequality constraints
        :type soccs: list of prbsym.XXXIneqConstr
        """
        self.smb = smb
        self.struct = self._parse_expr2struct(soccs)

    def sym2num(self, num):
        soccs = []
        for k, rng in enumerate(self.struct['rngs']):
            funccentres = self.struct['centres'][k]
            funcuppers = self.struct['uppers'][k]
            if rng is not None:
                rngini = eval(rng.ini, num)
                rngend = 1 + eval(rng.end, num)
                for i in range(rngini, rngend):
                    soccnum = dict()
                    num[rng.var] = i
                    soccnum['centres'] = funccentres.sym2num(num)
                    soccnum['uppers'] = funcuppers.sym2num(num)
                    soccnum['rows'] = funccentres._get_rows(num)
                    soccs.append(soccnum)
            else:
                print("Error: sym2num not implemented for not ranged socc")
        return soccs

    def _parse_expr2struct(self, soccs):
        centres = []
        uppers = []
        rngs = []

        for socc in soccs:
            centre = AffineFunction(self.smb, socc.centre.elems)
            upper = AffineFunction(self.smb, socc.upper.elems)

            rngs.append(socc.rng)
            centres.append(centre)
            uppers.append(upper)

        return dict(centres=centres, uppers=uppers, rngs=rngs)
