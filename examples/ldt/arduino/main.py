r"""
Generate the MPC code and copy it
into the muaompc directory. 
If the file muaompc/muaompc.ino does not exist,
copy a sample file. An existing muaompc/muaompc.ino will not
be replaced.
"""
import shutil

from muaompc import ldt
mpc = ldt.setup_mpc_problem('myprb.prb', numeric='float32')
ldt.generate_mpc_data(mpc, 'mydat.dat', muc=True)

try:
    shutil.os.mkdir('muaompc_arduino')
except FileExistsError:
    pass
else:
    shutil.copy('tmpl_muaompc_arduino.ino', shutil.os.path.join(
        'muaompc_arduino', 'muaompc_arduino.ino'))

for file in shutil.os.listdir(mpc.path['muc']):
    fname = shutil.os.path.join(mpc.path['muc'], file)
    shutil.copy(fname, 'muaompc_arduino')

print('The arduino sketch can be found in the folder muaompc_arduino')
