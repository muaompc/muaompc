#!/usr/bin/python3
"""Create the MPC controller MATLAB interface.

For the complete example, run the main_motor.m MATLAB script.
"""

from muaompc import ldt

mpc = ldt.setup_mpc_problem('myprb.prb')
ldt.generate_mpc_data(mpc, 'mydat.dat')
