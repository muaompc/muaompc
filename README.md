# muaompc

## Overview

μAO-MPC is a code generation software package for linear model predictive control. μAO-MPC generates highly portable C code tailored for embedded applications. The underlying optimization algorithm and its implementation explicitly consider many of the limitations and requirements of real-time embedded applications, and in particular of microcontroller applications: low memory footprint, deterministic execution time, only additions and multiplications are performed (no divisions, square roots, etc.), and support for fixed-point and floating- point arithmetic. μAO-MPC is developed at the Laboratory for Systems Theory and Automatic Control, is written in Python, and provides MATLAB/Simulink interfaces to the generated C code.

The MPC optimization algorithm is a quadratic program (QP) solver based on an augmented Lagrangian method (also called method of multipliers) combined with Nesterov’s fast gradient method. Other QP solvers can also be easily used.

The generated code has been tested in several platforms, including x86/AMD64 PCs, ARM Cortex-M microcontrollers, Lego Mindstorms NXT, and Arduino microcontrollers. 

## Installation

We have thoroughly tested μAO-MPC under Python 3.  To install type in a console:
```
pip install muaompc
```

Alternatively, to install from sources:
```
python setup.py install
```

For further details, see [the documentation](https://muaompc.readthedocs.io/en/latest/index.html).
You can find examples and tutorials in the [examples\ldt](https://github.com/muaompc/muaompc/tree/main/examples/ldt) folder.

## μAO-MPC 1.x

 This version is a complete reimplementation of the μAO-MPC core. This new version is backwards incompatible with, yet very similar to, the μAO-MPC 0.4.x (no longer being maintaned). Version 1.x includes a new type of optimization algorithm, offers more flexibility, and deals with many more types of MPC setups.

Here is the link to [Older versions](http://ifatwww.et.uni-magdeburg.de/syst/muAO-MPC/) (no longer being maintained).

## Citation

If you find this software useful, interesting, etc. please consider using citation [1]. Here is a corresponding BibTeX entry:

```
@INPROCEEDINGS{ifat:Zometa_12b,
author = {Zometa, P. and K\"ogel, M. and Findeisen, R.},
title = {{muAO-MPC}: A Free Code Generation Tool for Embedded Real-Time Linear Model Predictive Control},
booktitle = {Proc. American Control Conference ({ACC}), 2013},
year = {2013},
pages = {5340-5345},
address = {Washington D.C., USA},
pubtype = {proceedings}
} 
```

## License

This software is released under the 3-Clause BSD license.

