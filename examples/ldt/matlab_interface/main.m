%% Basic example
% This script first calls the (extenal) Python code generator.
% Afterwards, it compiles and uses the generated MPC controller.
% It uses the default solver, and shows how to use a different
% solver (quadprog) using the optimization problem data.
%% create the interface
!python main.py
%% compile it and make it available
cd myprb_mpc/src/matlab/
mpcmake
cd ../../..
addpath myprb_mpc/src/matlab/
%% Usage example
ctl = mpcctl('myprb_mpc/data/mydat/mpcmydat.json'); % create an instance of the class
ctl.conf.in_iter = 10;  % configure the algorithm
ctl.parameters.x_k = [0.1; -0.5]; % current state
% forming a QP is only necessary if a different QP solver is used
ctl.form_problem();  % form the QP for the current state
prb = ctl.prb;  % prb contains the QP in common form
u = quadprog(prb.H, prb.g, [], [], [], [], prb.u_lb, prb.u_ub);
% the QP is automatically formed when using its own algorithm (ALM+FGM) 
ctl.solve_problem();  % solve the MPC problem for current state
ctl.u_opt  % display the computed input sequence
norm(u - ctl.u_opt)  % ctl.u_opt is the MPC approximation of u
