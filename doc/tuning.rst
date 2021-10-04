.. _tuning:

******
Tuning
******


Basics of tuning
================

There are only two tuning parameters for the default 
optimization algorithm used by 
:math:`\muaompc`: 
the number of *internal* and *external*  iterations. We find that in 
many cases the tuning procedure is easy and intuitive.  For problems
without state constraints, only the number of internal iterations is 
of importance. These parameters are specified online. 

At the moment, the selection of these parameters is made entirely by the user.
In many embedded systems, the number of iterations may be limited by the 
processor computational power. More generally, the user may need to compare 
the MPC controller performance given by the solution
of an exact solver (like CVXOPT) against that given by the solution of 
:math:`\muaompc` for
a given number of iterations. For example, the comparison could be made
using the stage cost at each point of a given trajectory
(see [ZKF13]_). In the end, the precise number of iterations strongly depends on
the application.

.. _tuning.mu:

The penalty parameter
=====================

An optional third tuning value is the *penalty parameter* :math:`\mu`, 
which is selected off-line (i.e. specified in the data file).
:math:`\muaompc` will by default automatically compute a *good* value for 
:math:`\mu` if none is specified (recommended). 
Roughly speaking, a large penalty parameter implies that 
a low number of external iterations are required to reach good performance, 
especially when state constraints are active. However, more internal iterations
are necessary, because the condition number of the internal problem increases. 
The opposite is also true, a small :math:`\mu` makes the internal problem 
easier to solve, especially if no state constraints are active.  When 
the state constraint are active, however, the required number of external 
iterations is higher.  

By now it should be clear that the selection of an appropriate value of 
:math:`\mu` (not too low, not too high) is crucial. 

Although in general not recommended, :math:`\muaompc` allows
experienced users to explicitely set a value for :math:`\mu` 
in the data file. 
The selection of the penalty parameter :math:`\mu` is easily done via 
the function ``find_penalty_parameters`` in the ``ldt`` module. 
For example, using the the problem and data files 
from the tutorial::
   
   from muaompc import ldt
   res = ldt.find_penalty_parameters('myprb.prb', 'mydat.dat')

``res`` is a dictionary that contains the keys ``params`` and ``condnums``. 
In this case, each of them is a list
consisting of a single element. For ``params``, it is the value of :math:`\mu`
used by default, and for ``condnums`` is the condition number for the algorithm
corresponding to that parameter. Thus, to get default value of :math:`\mu` type::

   mu = res['params'][0]

The default value of the penalty parameter is one that is not too high 
but not too low. By using the parameter ``factors``, 
several values of multiples of the default :math:`\mu` can be tried at once::

   res = ldt.find_penalty_parameters('myprb.prb', 'mydat.dat', factors=[1, 4])

Now, ``res`` will contain the two penalty parameters with their corresponding 
condition number. 
For example, by typing::

   mu = res['params'][0]
   mu4 = res['params'][1]
   cn4 = res['condnums'][1]

we get in ``mu`` the default value for :math:`\mu` (i.e. :math:`\mu * 1`),
and ``mu4`` has the value :math:`\mu * 4`, corresponding to the second factor
in the list given as input via ``factors``.
This may help to check that the parameter ``mu4`` does not make the 
value of ``cn4`` too high (i.e. the internal 
problem is ill-conditioned). 

Once you have found a new value of ``mu`` that better suit your needs,
you need to include that information in your data file.
For example, if the new value of the penalty parameter is ``123``,
modify the ``mydat.dat`` data file by adding the line::

   mu = 123

Save the data file, and generate the data again as explained in 
the Section :ref:`tutor.basic` .
