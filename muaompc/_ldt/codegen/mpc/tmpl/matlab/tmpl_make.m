mex -outdir ./@{prefix}ctl -I../include -I./include {prefix}_matlab_create_ctl.c {c_src}

mex -outdir ./@{prefix}ctl -I../include -I./include {prefix}_matlab_form_problem.c {c_src}

mex -outdir ./@{prefix}ctl -I../include -I./include {prefix}_matlab_solve_problem.c {c_src}
