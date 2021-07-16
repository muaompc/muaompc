"""Make the symbolic representation of a parsed problem.

All the rules of how to create the symbolic representation of each of
the components parsed according to the DSL grammar.
It includess defining how the values of some of the problem variables
are to be handled,
in particular expanding ranges and evaluating indices.
It is independent of the shape of the matrices.
"""

import pickle

import numpy as np

from muaompc._ldt.parse.prbdsl import keywords

debug = True
symbol_list = None


class ProblemSym(object):

    """Symbolic representation of a parsed problem.
    """

    def __init__(self, prs):

        """Create a symbolic representation of a parsed problem.

        :param prs: a parsed problem using the function prbdsl.parse_file

        """
        prs = _list2dict(prs)
        global symbol_list
        symbol_list = SymbolList(prs)  # fixed symbols unknown
        self.prs = prs
        self.sub = self._get_sub(prs)
        self.min_ = self._get_min(prs)
        self.smb = SymbolElements(symbol_list)  # fixed symbols are known

    def _get_sub(self, prs):
        sub = Constr(prs[keywords['sub']])
        #ub._get_elements(prs[keywords['sub']])
        return sub

    def _get_min(self, prs):
        return Min(prs[keywords['min']])

def _list2dict(prs):
    # transform list of dicts into a dict
    dprs = dict()
    for elem in prs:
        dprs.update(elem)
    return dprs


def _load_symbolic_pickle(fname):
    with open(fname+'.pkl', 'rb') as f:
        prb = pickle.load(f)
    return prb


class SymbolList(object):
    def __init__(self, prs):

        """Get the symbols for all variables from a parsed file.

        Assign a list with the name (str) of the corresponding identifiers
        in the following fields:
        self.opt: the optimization variable
        self.aux: the auxiliary term
        self.par: the parameters
        self.fix: fix/constant elements

        :param prs: a parsed problem in dictionary form
        """

        self.par = prs[keywords['par']]  # parameters
        self.opt = prs[keywords['opt']]  # optimization variable
        self.aux = prs[keywords['aux']]  # auxiliary variable
        self.fix = set()  # fix values to be inferred. set avoids double entries


class SymbolElements(object):
    def __init__(self, lst):
        """From a set of symbols, create the corresponding classes.

        Make the following assignments:
        self.opt: a single Sequence of the optimization variable
        self.aux: a single/tuple of Sequence of the auxiliary term
        self.par: a tuple of Sequence of the parameters
        self.fix: a tuple of Elements of the fix/constant elements

        :param lst: a set of symbols
        :type lst: SymbolList
        """
        self.lst = lst
        self.opt = self.get_sequences(lst.opt)[0]
        self.aux = tuple(self.get_sequences(lst.aux))
        self.par = tuple(self.get_sequences(lst.par))
        self.fix = tuple(self.get_elements(list(lst.fix)))

    def get_elements(self, symbols):
        # parse 'x [a:b]', with the range [a:b] optional
        elements = []
        for pos, smb in enumerate(symbols):
            try:
                symbols[pos+1]
            except IndexError:
                rng = None
            else:
                if type(symbols[pos+1]) is dict:
                    rng = Ranged(symbols.pop(pos+1))
                else:
                    rng = None
            elements.append(RangedElement(smb, rng))

        return elements

    def get_sequences(self, symbols):
        # parse 'x [a:b] (n)', with the range [a:b] optional
        elements = []
        for pos, smb in enumerate(symbols):
            try:
                symbols[pos+2]
            except IndexError:
                rng = None
                veclen = symbols.pop(pos+1)
            else:
                if type(symbols[pos+1]) is dict:
                    rng = Ranged(symbols.pop(pos+1))
                    veclen = symbols.pop(pos+1)
                else:
                    rng = None
                    veclen = symbols.pop(pos+1)
            elements.append(Sequence(smb, veclen, rng))

        return elements


class Element(object):
    """Basic descriptor of all terms involved in the optimization problem."""
    def __init__(self, name, idx=None, preop=None):
        """
        :param name: the name of the element
        :type name: str
        :param idx: the identifier of the index of a vector in a sequence
        :type idx: str or None
        :param preop: left operand of the element, e.g. +, -, *
        :type preop: str or None
        """
        self.name = name
        self.idx = idx
        self.type_ = self._get_type(name)
        self.preop = preop

    def _get_type(self, name):
        global symbol_list
        if name in symbol_list.opt:
            type_ = keywords['opt']
        elif name in symbol_list.aux:
            type_ = keywords['aux']
        elif name in symbol_list.par:
            type_ = keywords['par']
        else:
            try:
                float(name)
            except ValueError:  # not a numerical value
                type_ = keywords['fix']
                symbol_list.fix.add(name)
            else:
                type_ = keywords['num']
        return type_


class RangedElement(Element, object):
    """Descriptor of terms that go throughout a range."""
    def __init__(self, name, rng=None):
        """
        Extends its parent class with the following assignments:
        self.rng: the range throughout which the elements go
        :param name: the name of the element
        :type name: str
        :param rng: a description of the range of the element
        :type rng: Ranged
        """
        self.rng = rng
        Element.__init__(self, name, idx=None)


class Sequence(RangedElement, object):
    """A sequence of vectors."""
    def __init__(self, name, veclen, rng):
        """
        Extends its parent class with the following assignments:
        self.veclen: the length of each vector in the sequence
        self.seqlen: the length of the full sequence
        Each is a simple scalar expression (str) that can be
        passed as argument to eval.

        :param name: the name of the sequence
        :type name: str
        :param veclen: a simple expression for the length of each vector
        :type veclen: str
        :param rng: a description of the range of the sequence
        :type rng: Ranged
        """
        self.veclen = str(veclen)  # length of each vector in sequence
        self.seqlen = self._get_seqlen(rng)  # number of elements in sequence
        RangedElement.__init__(self, name, rng)  # define number of vectors
        _add_elems_to_fix_symbol(self.veclen)

    def eval_seqlen(self, num):
        if self.seqlen is None:
            seqlen = eval(self.veclen, num)
        else:
            seqlen = eval(self.seqlen, num)
        return seqlen

    def _get_seqlen(self, rng):
        if rng is None:
            seqlen = None
        else:
            seqvec = str(rng.end + '-' + rng.ini + '+1')
            seqlen = self.veclen+'*('+ seqvec+')'
        return seqlen

def _add_elems_to_fix_symbol(expr):
    # add fixed values in expr (a str) to the global symbol_list
    t = expr[1:-1]  # remove parenthesis (first and last element)
    Expr(list(t))  # this adds fix to the global list of symbols

def _get_index(prs):
    return prs['idx']


class Ranged(object):
    """The range of an identifier."""
    def __init__(self, prs):
        """
        Make the following assignments:
        self.ini: a simple expression (str) of an integer describing the
        beginning of the range
        self.end: a simple expression (str) of an integer describing the
        end of the range
        These are both inclusive, i.e. the range will have
        (self.end - self.ini) + 1 elements
        :param prs: a parsed description of the range of the identifier
        """
        self.ini = self._get_ini(prs['range'])
        self.end = self._get_end(prs['range'])
        _add_elems_to_fix_symbol(self.ini)
        _add_elems_to_fix_symbol(self.end)

    def _get_end(self, prs):
        colon = prs.index(':')
        prs_end = prs[colon+1:]
        return ''.join(prs_end)

    def _get_ini(self, prs):
        colon = prs.index(':')
        prs_ini = prs[:colon]
        return ''.join(prs_ini)


class RangedIndex(Ranged, object):
    """The range of an index."""

    def __init__(self, prs):
        """Extend the parent class with the following assignments:
        self.var: a simple expression (str) of the index variable
        overload _get_ini method to consider the index variable in prs

        :param prs: a parsed description of the range of the index
        """
        prs_rng = self._get_prs(prs)
        Ranged.__init__(self, prs_rng)
        var = self._get_index_var(prs_rng['range'])
        self.var = ''.join(var)

    def _get_prs(self, prs):
        comma = prs.index(',')
        return prs[comma+1]  # after the comma

    def _get_ini(self, prs):
        colon = prs.index(':')
        eqsign = prs.index('=')
        prs_ini = prs[eqsign+1:colon]
        #ni = Expr(prs_ini)
        #eturn ini
        return ''.join(prs_ini)

    def _get_index_var(self, prs):
        equal = prs.index('=')
        return prs[:equal]


class Expr(object):
    """An expresion consisting of a list of Element instances.

    Intended to be subclassed.
    The method _get_elements should be overloaded.
    """
    def __init__(self, prs, rng=None):
        """
        self.elems: a list of instances of Element
        self.rng: the range of the expression, if any

        :param prs: a parsed description of the expression
        :param rng: a description of the range of the expression
        :type rng: RangedIndex or None
        """
        self.elems = []
        self.rng = rng
        self._get_elements(prs)

    def _get_elements(self, expr):
        ops = '+-*'
        self._prepend_plus_sign(expr)
        ns = dict()

        for pos, obj in enumerate(expr):
            if obj in ops:
                preop = obj
                name = expr.pop(pos+1)
                try:
                    expr[pos+1]
                except IndexError:
                    idx = None
                    index = ''
                else:
                    if type(expr[pos+1]) is dict:
                        idx = _get_index(expr.pop(pos+1))
                        index = '['+idx+']'
                        if self.rng is None:
                            pass
                    else:
                        idx = None
                        index = ''

                self.elems.append(
                        Element(name=name, idx=idx, preop=preop))
            else:
                print('error: unexpected element, type:', type(obj))
                print('the object itself:', obj)


    def _prepend_plus_sign(self, prs):
        ops = '+-'

        if type(prs) is not str:
            prs.insert(0, '+')
        elif prs[0] not in ops:
            prs.insert(0, '+')  # sign of first element if not specified
        self.preop = prs[0]


class Min(Expr, object):
    """The cost function to be minimized."""
    def __init__(self, prs):
        """
        Make the following assignments:
        self.elems: a list of Expr describing the cost function

        :param prs: a parsed description of the cost function
        """
        Expr.__init__(self, prs)

    def _get_elements(self, prs):
        ops = '+-'
        self._prepend_plus_sign(prs)
        for pos, obj in enumerate(prs):
            if obj in ops:
                preop = obj
                kw = prs.pop(pos+1)
                qua_expr = kw.get(keywords['qua'])
                sum_expr = kw.get(keywords['sum'])
                if qua_expr is not None:
                    t = Qua(qua_expr)
                    self.elems.append(t)
                elif sum_expr is not None:
                    t = Sum(sum_expr)
                    self.elems.append(t)
                else:
                    print('Min: keyword not understood')
            else:
                print('Min: Expression not understood')


class Sum(Expr, object):
    """Summation with range"""
    def __init__(self, prs):
        """
        It would typically look like:
        sum(f(x[i]), i=1:N)

        Make the following assignments:
        self.elems: a list of Expr describing f(x[i]) together with the range

        :param prs: a parsed description of the summation
        """
        Expr.__init__(self, prs, RangedIndex(prs))

    def _get_elements(self, prs):
        comma = prs.index(',')
        kw_sum = prs[:comma]  # before the comma
        ops = '+-'
        self._prepend_plus_sign(kw_sum)

        for pos, obj in enumerate(kw_sum):
            if obj in ops:
                preop = obj
                kw = kw_sum.pop(pos+1)
                qua_expr = kw.get(keywords['qua'])
                if qua_expr is not None:
                    q = Qua(qua_expr, rng=self.rng)
                    self.elems.append(q)
                else:
                    print('Sum: keyword not understood.',
                        'Only sum of quadratic terms is supported.')
            else:
                print('Sum: Expression not understood')


class Qua(Expr, object):
    """Quadratic form."""
    def __init__(self, prs, rng=None):
        """
        It would typically look like:
        f(x).T * g(y) * f(x)

        Make the following assignments:
        self.vec: the vector of the quadratic form, an instance of Expr
        self.mtx: the matrix of the quadratic form, an instance of Expr

        :param prs: a parsed description of the quadratic form
        :param rng: a description of the range of the expression
        :type rng: RangedIndex or None
        """
        Expr.__init__(self, prs, rng)

    def _get_elements(self, prs):
        comma = prs.index(',')
        prs_vec = prs[:comma]
        prs_mtx = prs[comma+1:]
        self.vec = Expr(prs_vec, rng=self.rng)
        self.mtx = Expr(prs_mtx, rng=self.rng)


class Constr(object):
    """Equality and inequality constraints."""
    def __init__(self, prs):
        """Sort the constraints into equality and inequality constraints.

        Make the following assignments:
        self.eqcs: a list with instances of equality constraints class
        self.iecs: a list with instances of inequality constraints class

        :param prs: a parsed description of all constraints
        """
        self.eqcs = []
        self.iecs = []
        self.soccs = []
        for c in prs:
            if c.get(keywords['eqc']) is not None:
                c = c[keywords['eqc']]
                try:
                    c.index(',')
                except ValueError:
                    self.eqcs.append(EqConstr(c))
                else:
                    self.eqcs.append(EqConstrRanged(c))

            elif c.get(keywords['iec']) is not None:
                c = c[keywords['iec']]
                try:
                    c.index(',')
                except ValueError:
                    self.iecs.append(IneqConstr(c))
                else:
                    self.iecs.append(IneqConstrRanged(c))
            elif c.get(keywords['socc']) is not None:
                c = c[keywords['socc']]
                try:
                    c.index(',')
                except ValueError:
                    self.soccs.append(IneqConstrSOC(c))
                else:
                    self.soccs.append(IneqConstrSOCRanged(c))


class EqConstr(Expr, object):
    """Simple equality constraint."""
    # x[0] = x_k
    def __init__(self, prs):
        """
        The typical form would be:
        x[0]=x_k
        With the following assignments:
        self.left = the left hand side of the equality
        self.right = the right hand side of the equality

        See the implementation of self._get_elements for details

        :param prs: a parsed description of the expression
        """
        Expr.__init__(self, prs)

    def _get_elements(self, prs):
        equal = prs.index('=')
        self.left = Expr(prs[:equal])
        self.right = Expr(prs[equal+1:])


class EqConstrRanged(Expr, object):
    """Ranged equality constraint."""
    # x[i+1]=A*x[i]+B*u[i], i=0:N-1
    def __init__(self, prs):
        """
        The typical form would be:
        x[i+1]=A*x[i]+B*u[i], i=0:N-1
        With the following assignments:
        self.left = the left hand side of the equality
        self.right = the right hand side of the equality
        self.rng = the range of the expression

        See the implementation of self._get_elements for details

        :param prs: a parsed description of the expression
        """
        Expr.__init__(self, prs, RangedIndex(prs))

    def _get_elements(self, prs):
        equal = prs.index('=')
        comma = prs.index(',')
        self.left = Expr(prs[:equal], rng=self.rng)
        self.right = Expr(prs[equal+1:comma], rng=self.rng)


class IneqConstr(Expr, object):
    """Simple doubled-sided inequality constraint."""
    def __init__(self, prs):
        """
        The typical form would be:
        a <= f(x[k]) <= b
        With the following assignments:
        self.lower = the lower bound of the inequality (a)
        self.upper = the upper bound of the inequality (b)
        self.centre = the centre or main body of the inequality (f(.))

        See the implementation of self._get_elements for details

        :param prs: a parsed description of the expression
        """
        Expr.__init__(self, prs)

    def _get_elements(self, prs):
        ineq = prs.index('<=')
        self.lower = Expr(prs[:ineq])
        if prs.count('<=') == 1:
            self.centre = Expr(prs[ineq+1:])
            self.upper = Expr('inf')  # no constraint, bound is infinity
        elif prs.count('<=') == 2:
            centre_upper = prs[ineq+1:]
            ineq = centre_upper.index('<=')
            self.centre = Expr(centre_upper[:ineq])
            self.upper = Expr(centre_upper[ineq+1:])


class IneqConstrRanged(Expr, object):
    """Ranged doubled-sided inequality constraint."""
# lower must be present, upper is optional. Both must be fixed or parameters
# centre must be an affine function of the optimization variable
# TODO: accept quadratic constraints

    def __init__(self, prs):
        """
        The typical form would be:
        a <= f(x[i]) <= b, i=0:N-1
        With the following assignments:
        self.lower = the lower bound of the inequality (a)
        self.upper = the upper bound of the inequality (b)
        self.centre = the centre or main body of the inequality (f(.))
        self.rng = the description of the range of the index

        See the implementation of self._get_elements for details

        :param prs: a parsed description of the expression and range
        """
        Expr.__init__(self, prs, RangedIndex(prs))

    def _get_elements(self, prs):
        ineq = prs.index('<=')
        comma = prs.index(',')
        self.lower = Expr(prs[:ineq], rng=self.rng)
        if prs.count('<=') == 1:
            self.centre = Expr(prs[ineq+1:comma], rng=self.rng)
            self.upper = Expr('inf')  # no constraint, bound is infinity
        elif prs.count('<=') == 2:
            centre_upper = prs[ineq+1:comma]
            ineq = centre_upper.index('<=')
            self.centre = Expr(centre_upper[:ineq], rng=self.rng)
            self.upper = Expr(centre_upper[ineq+1:], rng=self.rng)

class IneqConstrSOC(Expr, object):
    """Secon-order cone constraint."""
    def __init__(self, prs):
        """
        The typical form would be:

        With the following assignments:

        See the implementation of self._get_elements for details

        :param prs: a parsed description of the expression
        """
        Expr.__init__(self, prs)

    def _get_elements(self, prs):
        ineq = prs.index('<=')
        self.centre = Expr(prs[:ineq])
        self.upper = Expr(prs[ineq+1:])


class IneqConstrSOCRanged(Expr, object):

    def __init__(self, prs):
        """
        The typical form would be:
        With the following assignments:

        See the implementation of self._get_elements for details

        :param prs: a parsed description of the expression and range
        """
        Expr.__init__(self, prs, RangedIndex(prs))

    def _get_elements(self, prs):
        ineq = prs.index('<=')
        comma = prs.index(',')
        # FIXME: this is a quick and dirty implementation
# This assumes that centre consist only of one Norm2 element
        norm2_prs = prs[:ineq][0]
        self.centre = Expr(norm2_prs.get(keywords['norm2']), rng=self.rng)
        self.upper = Expr(prs[ineq+1:comma], rng=self.rng)
