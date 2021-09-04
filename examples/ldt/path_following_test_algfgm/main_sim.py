from mpc import mpcctl
from mydat import Ax, Bx, Az, Bz, dt
import numpy as np
from pandas import DataFrame as df
from matplotlib import pyplot as plt

print(df(Ax))
print(df(Bx))

ctl = mpcctl.Ctl('mpc_myprb/data/mydat/mpcmydat.json')
(n, m)  = Bx.shape
xk = np.zeros((n,1))
zk = np.zeros((n,1))
xk[0,0] = 0
zk[0,0] = 0
zk[1,0] = .05
Ns = 3000
x = np.zeros([n,Ns])
z = np.zeros([n,Ns])
u = np.zeros([m,Ns])
ctl.conf.warm_start = 1
ctl.conf.in_iter = 12
ctl.conf.ex_iter = 24
for k in range(Ns):
    x[:,k] = xk[:,0]
    z[:,k] = zk[:,0]
    ctl.parameters.x_k[:] = xk[:,0]
    ctl.parameters.z_k[:] = zk[:,0]
    ctl.solve_problem()
    uk = ctl.u_opt[0:m]
    uk = np.reshape(uk, (2,1))
    u[:,k] = uk[:,0]
    xk = Ax @ xk + Bx @ uk
    zk = Az @ zk + Bz @ uk
t = np.arange(0, Ns*dt, dt)
plt.subplot(3,1,1)
plt.plot(t, x[0,:], t, z[0,:])
plt.ylabel('position')
plt.subplot(3,1,2)
plt.ylabel('velocity')
plt.plot(t, x[1,:], t, z[1,:])
plt.subplot(3,1,3)
plt.ylabel('Input')
plt.plot(t, u[0,:], t, u[1,:])
plt.legend(['x', 'z'])
plt.show()