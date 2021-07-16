# muaompc

## Overview

μAO-MPC is a code generation software package for linear model predictive control. μAO-MPC generates highly portable C code tailored for embedded applications. The underlying optimization algorithm and its implementation explicitly consider many of the limitations and requirements of real-time embedded applications, and in particular of microcontroller applications: low memory footprint, deterministic execution time, only additions and multiplications are performed (no divisions, square roots, etc.), and support for fixed-point and floating- point arithmetic. μAO-MPC is developed at the Laboratory for Systems Theory and Automatic Control, is written in Python, and provides MATLAB/Simulink interfaces to the generated C code.

The MPC optimization algorithm is a quadratic program (QP) solver based on an augmented Lagrangian method (also called method of multipliers) combined with Nesterov’s fast gradient method. Other QP solvers can also be easily used.

The generated code has been tested in several platforms, including x86/AMD64 PCs, ARM Cortex-M microcontrollers, Lego Mindstorms NXT, and Arduino microcontrollers. 

## μAO-MPC 1.x

 This version is a complete reimplementation of the μAO-MPC core. This new version is backwards incompatible with, yet very similar to, the currently available μAO-MPC 0.4.x and older. Version 1.x includes a new type of optimization algorithm, offers more flexibility, and deals with many more types of MPC setups.

If your MPC problem fits into a setpoint stabilization or trajectory tracking problem and want a stable well tested code generation tool, use μAO-MPC 0.4.x. If you have a problem that does not fit those two MPC problems, or you are in for an adventure, try version 1.x.

Before you continue, please read the documentation.

[Full documentation (HTML)](http://ifatwww.et.uni-magdeburg.de/syst/research/muAO-MPC/doc/html-1.0/index.html) - Here you will find the required information to get you started, the supported MPC setup, installation instructions, tutorials, and function references.

[Full documentation (PDF)](http://ifatwww.et.uni-magdeburg.de/syst/research/muAO-MPC/doc/muaompc-1.0.pdf) - same as above, but in PDF.
    

## Older versions

Here is the link to [Older versions](http://ifatwww.et.uni-magdeburg.de/syst/muAO-MPC/)

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
