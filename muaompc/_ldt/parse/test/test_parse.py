#!/usr/bin/python3
from collections import namedtuple
from sys import path
import unittest
import shutil
from subprocess import call

import numpy as np
from numpy.testing import assert_allclose, assert_raises

from muaompc._ldt.parse import prbstruct
from muaompc._ldt.parse import prbdsl
from muaompc._ldt.parse.prbsym import Sequence as Sequence_
from muaompc._ldt.codegen.former.cvp import cvp


class TestErrorChecking(unittest.TestCase):

    def setUp(self):
        bnf = prbdsl.BNF()
        dsl, ddsl = bnf.get_dsl()
        self.ddsl = ddsl
        self.fname = 'testerrorchecking.prb'
        self.prb = ["variable u[0:N-1](m);",
                    "auxs x[0:N](n);",
                    "parameters x_k(n) xr[0:N](n) ur[0:N-1](m);",
                    "minimize sum(quad(x[i]-xr[i],Q)+quad(u[i]-ur[i], R), i=0:N-1)+quad(x[N]-xr[N],P);",
                    "subject to x[i+1] = A*x[i]+B*u[i], i=0:N-1;",
                    "x[0]=x_k;",
                    "u_lb <= u[i] <= u_ub, i=0:N-1;"]

    def tearDown(self):
        pass

    def test_setup_okay(self):
        with open(self.fname, 'w') as f:
            f.writelines([line+'\n' for line in self.prb])
        prbdsl.check_file(self.fname, self.ddsl, False) 

    def test_setup_errors(self):
        with open(self.fname, 'w') as f:
            f.writelines([line+'\n' for line in self.prb])
            f.writelines(self.prb[0]+'\n')  # variable defined twice
        assert_raises(prbdsl.SetupError,
                prbdsl.check_file, self.fname, self.ddsl, False) 
        with open(self.fname, 'w') as f:
            self.prb[0] = "asnouhntoenthuont"  # variable not defined
            f.writelines([line+'\n' for line in self.prb])
        assert_raises(prbdsl.SetupError,
                prbdsl.check_file, self.fname, self.ddsl, False) 


class TestSyntaxHighlight(unittest.TestCase):

    def setUp(self):
        self.fname = 'testerrorchecking.prb'
        self.prb = ["variable u[0:N-1](m);",
                    "auxs x[0:N](n);",
                    "parameters x_k(n) xr[0:N](n) ur[0:N-1](m);",
                    "minimize sum(quad(x[i]-xr[i],Q)+quad(u[i]-ur[i], R), i=0:N-1)+quad(x[N]-xr[N],P);",
                    "subject to x[i+1] = A*x[i]+B*u[i], i=0:N-1;",
                    "x[0]=x_k;",
                    "u_lb <= u[i] <= u_ub, i=0:N-1;"]

    def tearDown(self):
        pass

    def test_setup_okay(self):
        with open(self.fname, 'w') as f:
            f.writelines([line+'\n' for line in self.prb])
        with open(self.fname) as f:
            txt = ''.join(f.readlines())
        hgltxt = prbdsl.get_syntax_highlight(txt)


if __name__ == '__main__':
    unittest.main()
