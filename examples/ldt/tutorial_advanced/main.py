from muaompc import ldt 

mpc = ldt.setup_mpc_problem('myprb.prb')
ldt.generate_mpc_data(mpc, 'mydat.py', safe_mode=False)
