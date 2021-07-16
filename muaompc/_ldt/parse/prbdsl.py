"""Provide the grammar for the optimization problem. """
from copy import deepcopy

from pyparsing import (Keyword, Word, Forward, alphas, nums, alphanums, delimitedList,
        Regex, Literal, OneOrMore, ZeroOrMore, Empty, Group, Dict, ParseException, LineEnd,
                       Suppress)
from pyparsing import ParseFatalException, infixNotation, opAssoc

class SetupError(Exception): pass

keywords = dict(opt="variable", par="parameters", aux="auxs",
        stp='parameter', min="minimize", sub="subject to", sum='sum',
        qua='quad', eqc='eqconstr', iec='ineqconstr', fix='fix', num='numeric',
        norm2='norm2', socc='socc')

def parse_file(fname, verbose=False):
    bnf = BNF()
    dsl, ddsl = bnf.get_dsl()
    prs = None
    check_file(fname, ddsl, verbose)
    prs = dsl.parseFile(fname)
    return prs

def check_file(fname, ddsl, verbose):
    with open(fname) as f:
        lf = f.readlines()
        sf = ''.join(lf)
    optvar = ddsl['optvar']
    params = ddsl['params']
    auxvar = ddsl['auxvar']

    # Stage 1 - scan to check all expressions are there in the right amount
    n_optvar = count_dslrule(optvar, sf, verbose)
    if (n_optvar == 0):
        msg = ("There must be at least one definition of the optimization variable.")
        raise SetupError(msg)
    elif (n_optvar > 1):
        msg = ("There must be only one definition of the optimization variable.")
        raise SetupError(msg)

    # Stage 2 - check that each line contains a valid parsing expression
    for k, sf in enumerate(lf):
        n = 0
        for rule in ddsl.values():
            if is_dslrule_valid(rule, sf, verbose):
                n += 1
        if n == 0:
            msg = ('\nCould not parse line %d \nOffending line: %s' % (k, sf))
            raise SetupError(msg)
    return
    
def is_dslrule_valid(dslrule, sf, verbose):
    try:
        t = dslrule.parseString(sf)
    except:
        valid = False
    else:
        valid = True
        if verbose:
            print('DSL rule ', t)
    return valid 

def count_dslrule(dslrule, sf, verbose):
    count = 0
    for t, s, e in dslrule.scanString(sf):
        count += 1
        if verbose:
            print('Found DSL rule ', t)
    return count


class BNF(object):
    """The Backus-Naur form defining the domain-specific language."""
    def __init__(self):
        self.plus  = Literal( "+" )
        self.minus = Literal( "-" )
        self.mult  = Literal( "*" )
        self.div   = Literal( "/" )
        self.lpar  = Literal( "(" )
        self.rpar  = Literal( ")" )
        self.addop  = self.plus | self.minus
        self.multop = self.mult | self.div
        self.fnumber = Regex(
                r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")  # float
        self.var = Word(alphas, alphas+nums+"_")

    def get_dsl(self):
        # the rules of muaompc's dsl.
        term = self._term_arithmetic()
        term.setParseAction(lambda s,l,t: '('+t[0]+')')
        var = Word(alphas, alphanums+"_[]")
        seq = Word(alphas, alphanums+"_") + self._expr_var_range(term)*(0,1)
        vecseq = seq + self.lpar.suppress() + term + self.rpar.suppress()
        eos = Literal(";").suppress()  # end of statement
        optvar = Keyword(keywords['opt']) + vecseq + eos
        optvar.setParseAction(lambda s,l,t: {keywords['opt']:t[1:]})
        auxvar = Keyword(keywords['aux']) + Group(OneOrMore(vecseq)) + eos
        auxvar.setParseAction(lambda s,l,t: {keywords['aux']:t[1][:]})
        params = Keyword(keywords['par']) + Group(OneOrMore(vecseq)) + eos
        params.setParseAction(lambda s,l,t: {keywords['par']:t[1][:]})
        cost = Keyword(keywords['min']) + self._expr_cost() + eos
        cost.setParseAction(lambda s,l,t: {keywords['min']:t[1:]})
        rng = self._expr_time_range(term)
        expr = self._expr_simple()
        eqconstr = expr + Literal("=") + expr + ("," + rng)*(0,1) + eos
        eqconstr.setParseAction(lambda s,l,t: {keywords['eqc']:t[:]})
        ineqconstr = expr + Literal("<=") + expr + (Literal("<=") + expr)*(0,1) + ("," + rng)*(0,1) + eos
        ineqconstr.setParseAction(lambda s,l,t: {keywords['iec']:t[:]})
        socconstr = self._expr_norm2(expr) + (Literal("<=") + expr)*(0,1) + ("," + rng)*(0,1) + eos
        socconstr.setParseAction(lambda s,l,t: {keywords['socc']:t[:]})

        constr = Keyword(keywords['sub']) + OneOrMore(eqconstr | ineqconstr|
                socconstr)
        constr.setParseAction(lambda s,l,t: {keywords['sub']:t[1:]})
        dsl = optvar & auxvar & params & cost & constr
        ddsl = dict(optvar=optvar, params=params, auxvar=auxvar, 
                cost=cost, constr=constr, 
                eqconstr=eqconstr, ineqconstr=ineqconstr,
                socconstr=socconstr)
        return dsl, ddsl

    def _term_arithmetic(self):
        # term that is computed by simple arithmetic operations
        atom = ((0,None)*self.addop + (self.fnumber |
                self.var))
        term = infixNotation(atom, [
               #('-',1,opAssoc.RIGHT),
            (self.multop, 2, opAssoc.RIGHT, parse_action),
            (self.addop, 2, opAssoc.RIGHT, parse_action)
            ], self.lpar, self.rpar)
        return term

    def _term_index(self, term):
        index = Literal("[") + term  + Literal("]")
        index.setParseAction(lambda s,l,t: {'idx':''.join(t[1:-1])})
        return index

    def _expr_time_range(self, expr):
        idx = self.var
        rng = idx + "=" + expr + ":" + expr
        rng.setParseAction(lambda s,l,t: {'range':t[:]})
        rng.setFailAction(fail_action)  # FIXME: does not work as expected
        return rng

    def _expr_var_range(self, term):
        vrng = '[' + term + ":" + term + ']'
        vrng.setParseAction(lambda s,l,t: {'range':t[1:-1]})
        return vrng

    def _expr_sum(self):
        expr = self._expr_complex()
        term = self._term_arithmetic()
        term.setParseAction(lambda s,l,t: '('+t[0]+')')
        rng = self._expr_time_range(term)
        _sum = Keyword("sum") + self.lpar + expr + Literal(",") + rng + self.rpar
        _sum.setParseAction(lambda s,l,t: {'sum':t[2:-1]})
        _sum.setFailAction(fail_action)
        return _sum

    def _expr_quad(self, term):
        _quad = (Keyword("quad") + self.lpar + term + Literal(",") +
                    term + self.rpar)
        _quad.setParseAction(lambda s,l,t: {'quad':t[2:-1]})
        return _quad

    def _expr_norm2(self, expr):
        _norm2 = (Keyword("norm2") + self.lpar + expr + self.rpar)
        _norm2.setParseAction(lambda s,l,t: {'norm2':t[2:-1]})
        return _norm2

    def _expr_simple(self):
        term = self._term_arithmetic()
        index = self._term_index(term)
        var_elem = self.var + index*(0,1)
        expr = Forward()
        atom = ((0,None)*self.minus + (self.fnumber | var_elem))
        factor = Forward()
        factor << atom
        term_ = factor + ZeroOrMore((self.multop + factor))
        expr << term_ + ZeroOrMore((self.addop + term_))
        return expr

    def _expr_complex(self):
        expr_s = self._expr_simple()
        _quad = self._expr_quad(expr_s)
        atom = ((0,None)*self.minus + (_quad | expr_s))
        factor = Forward()
        factor << atom
# FIXME: this allows to have e.g. quad(.) * quad(.)
        term = factor + ZeroOrMore((self.multop + factor))
        expr = Forward()
        expr << term + ZeroOrMore((self.addop + term))
        return expr

    def _expr_cost(self):
        expr_c = self._expr_complex()
        sum_ = self._expr_sum()
        # FIXME: only sum is allowed to be the first element
        atom = ((0,None)*self.minus + (sum_ | expr_c))
        factor = Forward()
        factor << atom
        term = factor
        expr = Forward()
        expr << term + ZeroOrMore((self.addop + term))
        return expr

def fail_action(s, loc, expr, err):
    print("Error parsing: ", s.splitlines()[loc])
    print("Offending expression: ", expr)
    raise ParseFatalException

def fail_action_soft(s, loc, expr, err):
    print("Error parsing: ", s.splitlines()[loc])
    print("Offending expression: ", expr)
    return

def parse_action(s, l, t):
    return ''.join(t[0])

def get_keywords():
    # FIXME: there should only be one parsing function
    #text = OneOrMore(Word(alphas, alphanums+"_[]();:,.+-*/= \n"))
    ddsl = dict()
    for key in keywords.values():
        dsl = Keyword(key)
        ddsl[key] = dsl
    return ddsl

def get_syntax_highlight(text):
    #ith open(fname) as f:
    #   lf = f.readlines()
    lf = text.splitlines()

    ddsl = get_keywords()
    txt = [[{'normal':txt}] for txt in lf]

    hgltxt = []
    for k, sf in enumerate(lf):
        for rule in ddsl.values():
            for text, start, end in rule.scanString(sf):
                # This will only save last rule found
                txt[k] = [{'normal':sf[0:start]},
                        {'keyword':text[0]},
                        {'normal':sf[end:]}]
        txt[k].append({'endline': ""})
        for dic in txt[k]:
            hgltxt.append(dic)

    return hgltxt 
