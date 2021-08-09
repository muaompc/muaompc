.. _tutor.basic:

***************************
Code generation at a glance
***************************

The MPC problem description is written in a file called
the *problem* file. The data file for a specific problem
is given in a different text file called the *data* file.

After writing these files, the next step is to actually
auto-generate the C code. This is done in two easy steps:

#. create an mpc object from the problem file. 
#. generate code for the data using the newly created mpc object 
   and the data file. 
   
The first step will automatically generate the C-code for the problem. 
The second step will generate code for the data for two cases: 
static memory allocation and dynamic memory allocation. In the first case, 
the data consists on several C files that statically allocate memory and 
need to be compiled.
In the second case, a single data file in json format is written.
This data contained in the json file that can be dynamically loaded.
The first case (C-code) is useful for deployment in embedded systems,
whereas the second case (json file) allows more flexibility during
simulation (no need to compile the data).

*******************
A basic MPC problem
*******************

The code generation described in this section 
basically consist of the following steps:

#. write the problem file with the MPC problem description,
#. write a data file corresponding to the problem,
#. create a ``muaompc`` object out of that problem, and
#. create data from that object based on the data file.

The simplest problem that can be setup with the ``ldt`` module is an input
constrained problem. The code for this example can be found inside the 
*tutorial* directory ``muaompc_root/examples/ldt/tutorial``, 
where ``muaompc_root`` is the path to the root directory where 
``muaompc`` sources were extracted.

The MPC setup description 
=========================

.. default-role:: math

Consider the following setup. The plant to be controlled is described by the 
prediction model `x^+ = A x + B u`, where 
`x  \in \mathbb{R}^n`, and `u \in \mathbb{R}^m` are the 
current state and input vector, respectively. The state at the next
sampling time is denoted by `x^+`.
The discrete-time system and input matrices are denoted as 
`A` and `B`, respectively. 

The inputs are constrained to be in a box set `\mathcal{C}_u = \{u \; | \; u\_lb \leq u \leq u\_ub\}`.

The MPC setup is thus as follows:

.. math::
   \underset{\vb{u}}{\text{minimize}} & \;\; 
   \sum\limits_{i=0}^{N-1} (\| x_i\|^2_Q + \|u_i\|^2_R)  + \|x_{N}\|^2_P \\
   \text{subject to} 
   & \;\; x_{i+1}=A x_i + B u_i, \;\; i = 0, \cdots, N-1 \\
   & \;\; u\_lb \leq u_i \leq u\_ub, \;\; i = 0, \cdots, N-1   \\
   & \;\; x_0 = \bar{x} \\

where the integer `N \geq 2` is the prediction horizon.  The symmetric matrices
`Q`, `R`, and `P` are the state, input, and final state weighting matrices, 
respectively. The vector `\bar{x}` represents the system state at the current sampling time, which is given online as a parameter to the optimization problem.

The optimization variable `\vb{u} \in \mathbb{R}^{Nm}` 
is defined as the sequence `\vb{u} = \{u_0 \; u_1 \; \ldots \; u_{N-1}\}`.
Similarly, we define the (auxiliary) state sequence `\vb{x} = \{x_0\;  x_1\;  \ldots \; x_{N}\}`, with `\vb{x} \in \mathbb{R}^{(N+1)n}`.


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
    x[0]=x_bar;

Save the file as ``myprb.prb``.
The resemblance of this file to the mathematical description should be apparent.
Let us make a few remarks about the notation.

The problem description is based on *sequences*, which are defined
using the format `\texttt{v[a:b](m)}`, where `\texttt{v}` is the name of the
sequence, `\texttt{m}` is the length of each vector in the sequence, and `\texttt{a}`
and `\texttt{b}` denote the index of the first and last element of the sequence.
To refer to the element `\texttt{i}` of this sequence, we use the notation `\texttt{v[i]}`.
In the case that a sequence consist of a single element, only the vector length need
to be specified, i.e. `\texttt{v(m)}`.
To refer to this vector no indexing is required in later parts of the problem specification,
i.e. it suffices to write `\texttt{v}`.

The `\texttt{variable}` keyword identifies the optimization variable.
The `\texttt{parameters}` keyword identify sequences that are to be specified
online. 
The keyword `\texttt{auxs}` is used to specify the dimensions of the sequence `\texttt{x}`.

The keyword `\texttt{minimize}` identifies the text following it as the cost function to be minimized. Several
special keywords are accepted, for example the `\texttt{sum(h[i], i=a:b)}` denotes the summation of the
real valued functions `h_{i}` for `$i=a,\ldots, b$`. The keyword `\texttt{quad(v,M)}`
denotes the quadratic form `\|v\|_M^2`. 


As seen in problem file, the
equality and inequality constraints 
follow after the keyword `\texttt{subject to}`.
The constraints optionally accept an index variable and its range.
For example, the prediction model is written as `\texttt{x[i+1] = A*x[i] + B*u[i], i=0:N-1}`.


With the problem file already finished, we can now write the data file.

The MPC data file
=================

The data file must contain the numerical values for the matrices, vectors, 
and scalars found in the problem file. 
Values for the sequences defined by the keywords 
`\texttt{variable}`, `\texttt{parameters}` and `\texttt{aux}` 
are not required.
Thus, in our example, values for `\texttt{u}`, `\texttt{x\_bar}`, 
and `\texttt{x}` do *not* need to be specified.
The following elements need to be specified:
the matrices
`\texttt{Q}`, `\texttt{R}`, `\texttt{P}`,
`\texttt{A}`, `\texttt{B}`, 
the vectors
`\texttt{u\_lb}`, `\texttt{u\_ub}`, 
and the scalars
`\texttt{N}`, `\texttt{m}`, `\texttt{n}`.

The matrices are specified using MATLAB syntax. For example, an identity matrix `I \in \mathbb{R}^{2 \times 2}` could be written in the following ways::

    I = [1 0; 0 1]
    I = [1, 0; 0, 1]


Without going into further details, let us write data file. In your favourite text editor write::

   # weighting matrices
   Q = [1, 0; 0, 1]
   R = [1]
   P = [1, 0; 0, 1]
   # system matrices
   A = [1.,  0.01; 0.,  0.9]
   B = [1.e-04; 0.02]
   # input constraints
   u_lb = [-100]
   u_ub = [100]
   # dimensions
   N = 5
   n = 2
   m = 1


.. note::

    At the moment, each matrix or column vector must be described in a single line.


Save this file as ``mydat.dat``. 
The matrices `A` and `B` represent the discrete time model of a DC-motor.
The state vector is given by `x = [x_1 \; x_2]^T \in \mathbb{R}^n`, where `x_1` and `x_2` are the rotor
position and angular speed, respectively. The input is 
constrained to be between -100% and 100%. 

For this example, we chose the weighting matrices to be
identity matrices of appropriate size, i.e. `P = Q = I \in \mathbb{R}^{n
\times n}`, and `R = 1`. 
Clearly, the value of dimension of the state and input vector are 
`n = 2` and `m = 1`. 
The horizon length is specified as steps through the parameter `N=5`.


Generating the C-code
=====================

Now that we have written the ``myprb.prb`` problem file, 
we proceed to create an ``mpc`` object.
In the directory containing ``myprb.prb``,
launch your Python interpreter 
and in it type::

   from muaompc import ldt

   mpc = ldt.setup_mpc_problem('myprb.prb')

This will generate code specific for the problem described
by ``myprb.prb``.
By itself, the code we just generated is very not useful.
It only contains and abstract description of an MPC problem
without any data.
The next step is to generate code for data 
that can be used with the problem code 
for ``myprb.prb`` we just generated. 
To generate code that represents the data in ``mydat.dat``, 
continue typing in your Python interpreter::

   ldt.generate_mpc_data(mpc, 'mydat.dat')

And that's it! If everything went allright, you should now see inside current 
directory a new folder called ``mpc_myprb``. As an alternative to typing the 
above code, 
you can execute the file ``main.py`` found in the *tutorial* directory, 
which contains exactly that code. The *tutorial* directory already contains
the files ``myprb.prb`` and ``mydat.dat``.
In the next section, you will learn how to use the generated C code.

.. tip::
   If the code generation was not succesful, try passing the ``verbose=True``
   input parameter to the function ``setup_mpc_problem``. It will print extra
   information about the code generation procedure. For example:
   
   ``mpc = ldt.setup_mpc_problem('myprb', verbose=True)``

.. tip::
   By default, the generated code uses double precision float (64-bit) for all
   computations. You can specify a different numeric representation via
   the input parameter ``numeric`` of the function ``setup_mpc_problem``.
   For example, to use single precision (32-bit) floating point numbers type:
   
   ``mpc = ldt.setup_mpc_problem('myprb.prb', numeric='float32')``

Structure of the generated code
-------------------------------

In general terms, the generated code is structured as follows::

   + <prefix>_<prb_name>
     + src
       - C code + interfaces
     - <prefix>setup.py
     + data
       + <dat_name>
       + <dat_name_1>
       + ...

The folder where all the generated code is placed has a name in the form ``<prefix>_<prb_name>``, where ``<prefix>`` is a 
keyword argument passed to ``setup_mpc_problem``, and ``<prb_name>`` is the 
name of the problem file used to generate the code. The default prefix value is ``'mpc'``.
For example, to change the default prefix to something like ``xyz``, call:

   ``mpc = ldt.setup_mpc_problem('myprb.prb', prefix='xyz')``

For instance, in this tutorial the problem file is called ``myprb.prb``, 
and no prefix is specified (i.e. ``<prefix>=mpc``),
then ``<prb_name>=myprb``, and the directory for the generated code 
is ``mpc_myprb``.

Inside the ``src`` folder, the code for solving a problem are generated: the C-code, and
the Cython and MATLAB interfaces. All C-file names start with ``<prefix>``, which creates a sort of *name space*. 
This allows you to have several generated code coexist in a single application, as long as
each ``<prefix>`` is unique.
Similarly, the Cython and Matlab interfaces use the prefix as part of the interface name.

The ``<prefix>setup.py`` file is used to compile the Cython interface (see next section).

Finally, inside the ``<prefix>_<prb_name>`` folder, you will find the ``data`` folder. In it, you will find
the ``<dat_name>`` folder, which contains the generated data files for the ``<dat_name>.dat`` file. 
For each ``<dat_name>.dat`` MPC data file for which the call ``ldt.generate_mpc_data(mpc, '<dat_name>.dat')``
is made, a folder ``<dat_name>`` will be generated inside the ``data`` subfolder.  
This allows to generate different data sets 
(e.g. a ``<dat_name_1>.dat`` with different weighting matrices) for the same
problem.  This can be useful for controller tuning.

For instance, in this tutorial, inside the ``mpc_myprb`` folder, you will find the ``data`` folder, which in turn
contains the ``mydat`` folder.  ``mydat`` stores the generated data files for the data file ``mydat.dat``  
that correspond to the MPC problem ``myprb.prb``. 


Using the generated code
========================

In the folder ``mpc_myprb`` you will find all the automatically 
generated code for the current example.  
To use the generated code in a control loop, the following steps are to be followed:

#. setup a MPC controller
#. configure the optimization algorithm
#. set the parameters for the MPC controller
#. solve the MPC problem
#. apply the control input
#. repeat from step 3

We now proceed to exemplify the use of the generated code from
steps 1 to 5.
We start our tutorial using the Python interface, as it is simpler to
explain. Later we show how it is done in pure C, and using the MATLAB interface.


Using the generated code in Python 
----------------------------------
   
The Python interface makes it possible to 
directly make use of the generated code and data (i.e. the MPC controller)
from within Python.

Once the code has been generated,
the next step is to compile the Python interface.
Technically, we use Cython to define a C-extension to Python. 

In a console/terminal change to *tutorial* directory ``muaompc_root/examples/ldt/tutorial``. 
Change then to the generated code folder ``mpc_myprb``. 
To install the Python extension, execute the ``mpcsetup.py`` installation script::

   python mpcsetup.py install --force

If everything went ok, you should see no errors, and the last three lines 
should be (tested in Ubuntu 20.04):: 

   Installed <>.egg
   Processing dependencies for mpc==1.0
   Finished processing dependencies for mpc==1.0

where  ``<>`` is a general place holder.

Now you can use the interface which is encapsulated in a package called 
``mpc``  which represents the MPC controller.  In general, the Python package's name
is the same as the ``<prefix>`` used during code generation.

While in the folder ``mpc_myprb``, fire up your Python interpreter, and type::

   from mpc import mpcctl

The next step is to declare an 
instance of the class ``mpcctl.Ctl``, which we usually call ``ctl`` (*controller*). 
The input parameter for the constructor of the class is the name
of a json file contaning the generated data.
In this example, the data is saved in the folder
``mpc_myprb/data/mydat``. In our example,
the generated json data file is called ``mpcmydat.json``.
Continue typing in the console::

   ctl = mpcctl.Ctl('data/mydat/mpcmydat.json')

The next step is to configure the optimization algorithm. 
In this case, we have an input
constrained problem. The only parameter to configure is the number of iterations of
the algorithm
(see section :ref:`tuning` for details).
For this simple case, let's set it to 10 iterations::

   ctl.conf.in_iter = 10; 
   
Let us assume that the current state is `\bar{x} = [0.1 \; -0.5]^T`. 
The controller object has a field for the parameters defined in the problem file. The parameter ``x_bar`` can be set as follows::

   ctl.parameters.x_bar[:] = [0.1, -0.5]

We can finally
solve our MPC problem for this state by calling::

   ctl.solve_problem();
   
The solution is stored in an array ``ctl.u_opt``, whose first ``m`` elements are
commonly applied to the controlled plant.
Print the optimal input vector ``ctl.u_opt``, and if everything went okay, 
you should see the following::

   print(ctl.u_opt)
   array([0.03056814, 0.02406793, 0.0178332 , 0.01179073, 0.00586953])


Using the generated code in MATLAB 
----------------------------------

The MATLAB interface makes it possible to 
directly make of the generated code and data (i.e. the MPC controller)
from within MATLAB.

Once the code has been generated,
the next step is to compile the MATLAB interface. 

Start MATLAB, and switch to the folder
``mpc_myprb/src/matlab``.
In the MATLAB console type ``mpcmake``, which will execute the ``mpcmake.m`` script. 
The last step is to add the ``matlab`` directory to the PATH environment in 
MATLAB.
For example, 
assuming the MATLAB current directory is 
the tutorial directory ``muaompc_root/examples/ldt/tutorial``, in the MATLAB console type::

   cd mpc_myprb/src/matlab 
   mpcmake
   cd ..
   addpath matlab

Now you can use the interface which is encapsulated in a class called 
``mpcctl``  which represents the MPC controller. Simply declare an 
instance of that class, which we usually call ``ctl`` (*controller*). 
The input parameter for the constructor of the class is the name
of a json file contaning the generated data.
``muaompc`` by default saves the data in the folder
``mpc_myprb/data/mydat``. In our example,
the generated json data file is called ``mpcmydat.json``.
Continue typing in the console::

    ctl = mpcctl('mpc_myprb/data/mydat/mpcmydat.json'); 

The next step is to configure the optimization algorithm. 
In this case, we have an input
constrained problem. The only parameter to configure is the number of iterations of
the algorithm
(see section :ref:`tuning` for details).
For this simple case, let's set it to 10 iterations::

   ctl.conf.in_iter = 10; 
   
Let us assume that the current state is `\bar{x} = [0.1 \; -0.5]^T`. 
The controller object has a field for the parameters defined in the problem file. The parameter ``x_bar`` can be set as follows::

   ctl.parameters.x_bar = [0.1; -0.5];

We can finally
solve our MPC problem for this state by calling::

   ctl.solve_problem();
   
The solution is stored in an array ``ctl.u_opt``, whose first `m` elements are
commonly applied to the controlled plant.
The complete MATLAB example can be found in the tutorial folder under ``main.m``. 


Using the generated code in C 
-----------------------------

The folder ``mpc_myprb/data/mydat`` already contains a template 
for a main file, called ``mpcmydatmain.c``.
Switch to the the folder ``mydat`` and open ``mpcmydatmain.c``
in your favourite editor.
This template file shows how to solve an MPC problem using dynamic 
or static memory allocation. 
This file might look at bit daunting at first, but it just a template
you can modify to fit your needs.

In the current directory you will find two main files with just the basics,
that are based on the template file.
The file ``mpcmydatmain_dynmem.c`` exemplifies how the dynamic memory allocation
is done in C code. The file ``mpcmydatmain_staticmem.c`` exemplifies the
how the static memory allocation version of the data can be used. 
Both files follow the same structure as the MATLAB 
tutorial above. 

Let us take the file ``mpcmydatmain_staticmem.c`` as example. 

The first thing to include
is the header file of the library called ``mpcctl.h``.

We need to have access to some of the algorithm's variables, for example the MPC 
system input, the parameters, and the algorithm settings. This is done through the fields of the 
``struct mpc_ctl`` structure, which we denote the *controller* structure.
We first create an instance of this controller structure, and we set the controller by passing a pointer to the structure to the function ``mpcmydat_ctl_setup_ctl``, which is found in ``mpcmydatctldata.h``. 
For example, after including the corresponding headers, in the body of the main function we type::

    struct mpc_ctl ctlst;  /* Structure for static memory allocation */
    struct mpc_ctl *ctl;  /* pointer to the an allocated structure */

    ctl = &ctlst;
    mpcmydat_ctl_setup_ctl(ctl);

Once we controller is setup, we can continue in a similar fashion to the MATLAB case, that is first we setup the parameters, then we configure the algorithm, solver the problem::

    ctl->parameters->x_bar[0] = 0.1;
    ctl->parameters->x_bar[1] = -0.5; 
    ctl->solver->conf->in_iter = 10;
    mpc_ctl_solve_problem(ctl);

Finally, the computed control input is found in the array ``ctl->u_opt``.

.. note::

   At the moment the user needs to know the length of the different arrays in the controller structure. This information can be infered by the user from the problem and data files. The length of the different arrays will be available in the controller structure in future releases.

To run an compile this code do the following. Copy the file ``mpcmydatmain_staticmem.c`` into the folder ``mpc_myprb/data/mydat`` and remove the main template file  ``mpcmydatmain.c`` found in ``mpc_myprb/data/mydat`` (otherwise you will end up with two ``main`` functions in two different files, and the compilation will fail).  In that folder you will also find example 
Makefiles, called ``mpcmydatMakefile.*``, which compiles the
generated code. 
The Makefile ``mpcmydatMakefile.mk`` compiles the code
using the GNU Compiler Collection (*gcc*).
Adapt the Makefile to your compiler if necessary.

For example, to generate, compile and run the code in Linux you would type 
in a console::
   
   cd muaompc_root/examples/ldt/tutorial  # the tutorial folder
   python main.py  # generates code and data
   cp mpcmydatmain_staticmem.c mpc_myprb/data/mydat  # the tutorial main file
   cd mpc_myprb/data/mydat
   rm mpcmydatmain.c  # remove template main file
   make -f mpcmydatMakefile.mk  # compile
   ./main  # run the controller

If everything went okay, you will see the output::

   ctl->u_opt[0] = 3.056814e-02

This concludes our tutorial!

