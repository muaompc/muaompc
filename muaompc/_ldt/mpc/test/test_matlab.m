cd ~/workspace/muaompc/install_mpcmat/mpcmat/matlab
!make -f mpcmatMakefile.mk
cd ..
addpath matlab/
c = mpcmatctl;
c.parameters.x_k = [2,3];
c.form_problem;
if (abs(c.prb.g-[0.1342, 0.1241]') < 1e-4)
    disp('Test OK')
else
    disp('Test Fail')
end