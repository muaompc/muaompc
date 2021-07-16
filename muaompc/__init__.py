"""
This package automatically generates self-contained C-code specific for 
a model predictive control (MPC) problem. 
The generated code consist of C-code for the problem, 
and C-code for data corresponding to that problem.
The problem, as well as the data, are specified by the user
in two separate text files using a high-level description language. 
The generated C-code is a ready-to-use fast MPC
implementation.

The generated C code is fully compatible with the ISO C89/C90 standard, 
and is platform independent. The code can be directly used in embedded 
applications, using popular platforms like
`Arduino <http://www.arduino.cc/>`_, and
`Raspberry Pi <https://www.raspberrypi.org/>`_,
or any other application on which C/C++ code is accepted, like many current 
generation programable logic controllers (PLC).
Additionally, MATLAB/Simulink interfaces to the generated code are provided.

At the moment, we consider the following types of systems:

* linear discrete-time sytems (the ``ldt`` module)
"""

from muaompc import ldt 
from muaompc.version import version
