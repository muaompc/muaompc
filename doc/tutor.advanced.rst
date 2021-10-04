.. _tutor.advanced:

**************************
A more complex MPC problem
**************************

In this section we consider a more elaborated example. However, the procedure to
follow is just as discussed in :ref:`tutor.basic`: describe the problem, generate C-code from it, and finally 
use the generated code.

We now consider a problem that presents many of the features
available in ``muaompc``. The code for this example can be found 
inside the *tutorial advanced* directory ``muaompc_root/examples/ldt/tutorial_advanced``, 
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

Observe that the system in the ``muaompc`` problem description is in discrete-time,
i.e. ``x[i+1] = A*x[i]+B*u[i], i=0:N-1;``.
The continuous-time matrices `A_c` and `B_c` of the initial formulation 
need therefore to be discretized. Besides the ``.dat`` file format presented in the
basic tutorial, a python or matlab file can be used as data input.
In this tutorial, we use a Python ``.py`` module as our data file.
As it is a regular Python file, we can, among other things, import ``scipy`` to helps discretized
the continuous-time system matrices using a zero-order hold.
All matrices must be defined as 2-dimensional numpy arrays (e.g. ``R = np.array([[472.]])``).
For this example, the data file looks like::

   import numpy as np
   from scipy.signal import cont2discrete as c2d
   # weighting matrices
   Q = np.diag([1014.7, 3.2407, 5674.8, 0.3695, 471.75])
   R = np.array([[472.]])
   P = Q
   # system matrices (continuos time)
   Ac = np.array([[-1.2822, 0, 0.98, 0], [0, 0, 1, 0], [-5.4293, 0, -1.8366, 0], [-128.2, 128.2, 0, 0]])
   Bc = np.array([[0.3], [0], [-17], [0]])
   Cc = np.array([[0, 1, 0, 0], [0, 0, 0, 1], [-128.2, 128.2, 0, 0]])
   Dc = np.zeros((3,1))
   # discretization
   dt = 0.5
   (A, B, C, D, dt) = c2d((Ac, Bc, Cc, Dc), dt)
   # The system is extended with a 5th state, to account for slew rate constraints
   # Extend A from a 4x4 matrix, to a 5x5 matrix
   # Add first a row of zeros:
   A = np.concatenate((A, np.zeros((1,4))))
   # then a column of zeros
   A = np.concatenate((A, np.zeros((5,1))), axis=1)
   # Extend B from a 4x1 column vector, to a 5x1 column vector
   # Add an new row at the bottom, with the element = 1:
   B = np.concatenate((B, np.ones((1,1))))
   # input constraints
   u_lb = np.array([[-0.262]])
   u_ub = np.array([[0.262]])
   # state constraints
   e_lb = np.array([[-0.349, -30, -0.25]]).T
   e_ub = -1*e_lb
   Kx = np.array([[0, 1, 0, 0, 0], [-128.2, 128.2, 0, 0, 0], [0., 0., 0., 0., -1.]])
   Ku = np.array([[0, 0, 1]]).T
   # terminal state constraints
   f_lb = e_lb
   f_ub = e_ub
   Kf = Kx 
   # dimensions
   N = 10  # horizon length
   n = 5  # number of states
   m = 1  # number of inputs



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
To generate code that represents the data in ``mydat.py``, 
continue typing in your Python interpreter::

   ldt.generate_mpc_data(mpc, 'mydat.py', safe_mode=False)


Note that to use a Python module like ``mydat.py`` as input, we must set ``safe_mode=False``.
By default, ``safe_mode=True``, and only the simple ``.dat`` format is accepted.

If everything went allright, you should now see inside current 
directory a new folder called ``mpc_myprb``. As an alternative to typing the 
above code, 
you can execute the file ``main.py`` found in the *tutorial_advanced* directory, 
which contains exactly that code. The *tutorial advanced* directory already contains
the files ``myprb.prb`` and ``mydat.py``.
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
and 2 *external* iterations provide an acceptable approximation of the MPC problem using the warm-start strategy::

   ctl.conf.in_iter = 24; /* number of internal iterations */
   ctl.conf.ex_iter = 2; /* number of external iterations */
   ctl.conf.warm_start = 1; /* automatically warm-start algorithm */


Using the generated code in Python 
----------------------------------

Just as in the *tutorial* example, we can use the 
Python interface to test our algorithm. 
Let's try doing the same using the Python interface.
Go to the to the *tutorial_advanced* directory,
change to the generated code folder ``mpc_myprb``, 
and install the Python extension::

   python mpcsetup.py install --force

Finally launch your Python interpreter, and in it type::

  from mpc import mpcctl
  ctl = mpcctl.Ctl('data/mydat/mpcmydat.json')
  # controller solver configuration
  ctl.conf.in_iter = 24 
  ctl.conf.ex_iter = 2 
  ctl.conf.warm_start = 1 
  # set current state
  ctl.parameters.x_bar[:] = [0., 0., 0., -400., 0.]
  # get solution
  ctl.solve_problem()

The optimal input should be::

  print(ctl.u_opt)
  [-0.262      -0.13715448  0.00099704  0.03492525  0.05084397  0.04797471
  0.04002588  0.03289448  0.02473329  0.01478007]

This concludes the advanced tutorial.