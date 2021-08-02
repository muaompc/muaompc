.. _tutor.advanced:

**************************
A more complex MPC problem
**************************

In this section we consider a more elaborated example. However, the procedure to
follow is the same: describe the problem, generate C-code from it, and finally 
use the generated code.


We now consider a problem that presents many of the features
available. The code for this example can be found 
inside the *tutorial* directory ``muaompc_root/examples/ldt/tutorial``, 
where ``muaompc_root`` is the path to the root directory of ``muaompc``.


The MPC setup description
=========================

.. default-role:: math

The system considered is the Cessna Citation 500
aircraft presented in ([M02]_, p.64).  A continuous-time
linear model is given by `\dot{x} = A_c x + B_c u, y = C x`, where

.. math::
   A_c = \left[ \begin{matrix}
   -1.2822 & 0 & 0.98 & 0 \\
   0 & 0 & 1 & 0 \\
   -5.4293 & 0 & -1.8366 & 0 \\
   -128.2 & 128.2 & 0 & 0 \\
   \end{matrix} \right], \;\;
   B_c = \left[ \begin{matrix}
   -0.3 \\
   0 \\
   -17 \\
   0 \\
   \end{matrix} \right],
   C = \left[ \begin{matrix}
   0 & 1 & 0 & 0 \\
   0 & 0 & 0 & 1 \\
   -128.2 & 128.2 & 0 & 0 \\
   \end{matrix} \right],

and the state vector is given by `x = [x_1 \; x_2 \; x_3 \; x_4]^T`, where:

* `x_1` is the angle of attack (rad),
* `x_2` is the pitch angle (rad),
* `x_3` is the pitch angle rate (rad/s), and
* `x_4` is the altitude (m).

The only input `u_1` is the elevator angle (rad).
The outputs are `y_1 = x_2`,  `y_2 = x_4`, and `y_3 = -128.2 x_1 + 128.2 x_2`
is the altitude rate (m/s)

The system is subject to the following constraints:

* input constraints `-0.262 \leq u_1 \leq 0.262`,
* slew rate constraint in the input `-0.524 \leq \dot{u}_1 \leq 0.524`
* state constraints `-0.349 \leq x_2 \leq 0.349`,
* output constraints `-30.0 \leq y3 \leq 30.0`.

To consider the slew rate constraint in the input, we introduce an additional
state `x_5`. The sampling interval is `dt = 0.5` s, and the
horizon length is `N = 10` steps.

The controller parameters
-------------------------

The book [M02]_ proposes to use identity matrices of appropriate size for
the weighting matrices `Q` and `R`. We instead select them diagonal
with values that give a similar controller performance and much lower
condition number of the Hessian of the MPC quadratic program,
a desirable property for any numerical algorithm.


The MPC problem file
====================

The MPC setup can be rather intuitively described in the problem file.
In your favorite text editor write the following::

    variable u[0:N-1](m);
    auxs x[0:N](n);
    parameters x_bar(n);
    minimize sum(quad(x[i],Q)+quad(u[i], R), i=0:N-1)+quad(x[N],P);
    subject to x[i+1] = A*x[i]+B*u[i], i=0:N-1;
    u_lb <= u[i] <= u_ub, i=0:N-1;
    e_lb <= Kx*x[i] + Ku*u[i] <= e_ub, i=0:N-1;
    f_lb <= Kf*x[N] <= f_ub;
    x[0]=x_bar;


With the problem file already finished, we can now write the data file.

The MPC data file
=================

Without going into further details, let us write data file. In your favourite text editor write::

    # weighting matrices
    # Q is a 5x5 diagonal matrix
    Q = [1014.7, 0, 0, 0, 0; 0, 3.2407, 0, 0, 0; 0, 0, 5674.8, 0, 0; 0, 0, 0, 0.3695, 0; 0, 0, 0, 0, 471.75]
    R = [471.65]
    # P is a copy of Q
    P = [1014.7, 0, 0, 0, 0; 0, 3.2407, 0, 0, 0; 0, 0, 5674.8, 0, 0; 0, 0, 0, 0.3695, 0; 0, 0, 0, 0, 471.75]
    # system matrices (discrete time)
    A = [  0.23996015,   0., 0.17871287, 0., 0.; -0.37221757, 1., 0.27026411, 0., 0.; -0.99008755, 0., 0.13885973, 0., 0.; -48.93540655, 64.1, 2.39923411, 1., 0.; 0., 0., 0., 0., 0.]
    B = [-1.2346445; -1.43828223; -4.48282454; -1.79989043; 1.]
    # input constraints
    u_lb = [-0.262]
    u_ub = [0.262]
    # state constraints
    e_lb = [-0.349; -30; -0.25]
    e_ub = [0.349; 30; 0.25]
    Kx = [0, 1, 0, 0, 0; -128.2, 128.2, 0, 0, 0; 0., 0., 0., 0., -1.]
    Ku = [0; 0; 1]
    f_lb = [-0.349; -30; -0.25]
    f_ub = [0.349; 30; 0.25]
    Kf = [0, 1, 0, 0, 0; -128.2, 128.2, 0, 0, 0; 0., 0., 0., 0., -1.]
    # dimensions
    N = 10
    n = 5
    m = 1


.. note::

    At the moment, each matrix or column vector must be described in a single line.


Generating the C-code
=====================

Similarly to the :ref:`tutor.basic`, we proceed to create an ``mpc`` object.
In the directory containing ``myprb.prb``,
launch your Python interpreter 
and in it type::

   from muaompc import ldt

   mpc = ldt.setup_mpc_problem('myprb.prb')

This will generate code specific for the problem described
by ``myprb.prb``.
The next step is to generate code for data 
that can be used with the problem code 
for ``myprb.prb`` we just generated. 
To generate code that represents the data in ``mydat.dat``, 
continue typing in your Python interpreter::

   ldt.generate_mpc_data(mpc, 'mydat.dat')

And that's it! If everything went allright, you should now see inside current 
directory a new folder called ``myprb_mpc``. As an alternative to typing the 
above code, 
you can execute the file ``main.py`` found in the *tutorial_advanced* directory, 
which contains exactly that code. The *tutorial* directory already contains
the files ``myprb.prb`` and ``mydat.dat``.
In the next section, you will learn how to use the generated C code.


Using the generated C-code
==========================

The next step is to make use of the generated code. For further
details on the generated code see :ref:`tutor.basic`.

Algorithm configuration
-----------------------

The next step is to configure the algorithm. In this case, we have a system
with input and state constraints. The only parameters to configure are the 
number of iterations of the algorithm. The state constrained algorithm is an 
augmented Lagrangian method, which means it requires a double iteration loop 
(an *internal* and an *external* loop). From simulation
we determine that 24 *internal* iterations,
and 2 *external* iterations provide an acceptable approximation of the MPC problem using the warmstart strategy::

   ctl.conf.in_iter = 24; /* number of internal iterations */
   ctl.conf.ex_iter = 2; /* number of external iterations */
   ctl.conf.warmstart = 1; /* automatically warmstart algorithm */


Using the generated code in Python 
----------------------------------

Just as in the simpler tutorial example, we can use the 
Python interface to test our algorithm. 
Let's try doing the same using the Python interface.
Go to the to the *tutorial_advanced* directory,
change to the generated code folder ``myprb_mpc``, 
and install the Python extension::

   python mpcsetup.py install --force

Finally launch your Python interpreter, and in it type::

  from mpc import mpcctl
  ctl = mpcctl.Ctl('data/mydat/mpcmydat.json')
  # controller solver configuration
  ctl.conf.in_iter = 24; 
  ctl.conf.ex_iter = 2; 
  ctl.conf.warmstart = 1; 
  # set current state
  c.x_bar[:] = [0., 0., 0., -400., 0.]
  # get solution
  ctl.solve_problem();
  ctl.u_opt
