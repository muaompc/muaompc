%% Tutorial example
%% Call the (extenal) Python code generator to create the interface
% If the following line does not work, 
% you might need to generate the code from a console outside MATLAB
!python main.py  
%% compile the generated code and make it available
cd myprb_mpc/src/matlab/
mpcmake
cd ../../..
addpath myprb_mpc/src/matlab/
%% Usage example
% create an instance of the class from generated data
ctl = mpcctl('myprb_mpc/data/mydat/mpcmydat.json'); 
ctl.conf.in_iter = 10; 
ctl.parameters.x_bar = [0.1; -0.5]; % current state
ctl.solve_problem();  % solve the MPC problem for current state x_bar
ctl.u_opt  % display the computed input sequence
