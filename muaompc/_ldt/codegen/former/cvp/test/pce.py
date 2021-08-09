from os import path
import numpy as np
from scipy import io
from numpy import diag
dt = 0.5
N = 10
matpath = path.join(path.dirname(__file__), 'qp_no_pce.mat')
mm = io.loadmat(matpath)
Q = mm['Q']
Q = np.array(diag([1014.7, 3.2407, 5674.8, 0.3695, 471.75]))
R = mm['R']
R = np.array(diag([4716.5]))
Ad = mm['A']
Ad = np.array([[  0.23996015,   0., 0.17871287,   0., 0.],
                   [ -0.37221757,   1., 0.27026411,   0., 0.],
                   [ -0.99008755,   0., 0.13885973,   0., 0.],
                   [-48.93540655, 64.1, 2.39923411,   1., 0.],
                   [0., 0., 0., 0., 0.]])
Bd = mm['B']
Bd = np.array([[-1.2346445 ],
                   [-1.43828223],
                   [-4.48282454],
                   [-1.79989043],
                   [1.]])
(n, m) = np.array(Bd).shape

P=Q
# input constraints
eui = 0.262  # rad (15 degrees). Elevator angle.
u_lb = [[-eui]]
u_ub =  [[eui]]

# mixed constraints
ex2 = 0.349  # rad/s (20 degrees). Pitch angle constraint.
ex5 = 0.524 * dt  # rad/s * dt input slew rate constraint in discrete time
ey3 = 30.
# bounds
e_lb = [[-ex2], [-ey3], [-ex5]]
e_ub = np.array([[ex2], [ex5]])
e_lb = -e_ub
(ncx, dummy) = e_ub.shape
# constraint matrices
# ss = io.loadmat('socp_matrices.mat')
# ssVv = ss['V21']
# ssVv2 = ssVv[42:47,42:47]  # This are the covariance weights, which are in an odd position
# beta = 0.05
# kappa_inv = 1/np.sqrt((1.-beta)/beta)
# nsocc = 1
Kx = np.zeros((ncx, n))
Kx[0, 1] = 1.
Kx[1, 4] = -1.
# Mm = np.zeros((nsocc, n))
# Mm[0,6] = 1.*kappa_inv  # The median of x2
# Vv = np.zeros((n, n))
# Vv[7:12,7:12] = ssVv2  # The covariance of x2 (5 terms)
# neg_c_max = 0.3490*kappa_inv
# terminal state constraints
f_lb = e_lb
f_ub = e_ub
F = Kx
data = dict(A=np.array(Ad), B=np.array(Bd), P=np.array(P), Q=np.array(Q), R=np.array(R), n=n, m=m, N=N,
u_lb=np.array(u_lb), u_ub=np.array(u_ub),
e_lb=np.array(e_lb), e_ub=np.array(e_ub),
Kx=np.array(Kx),
f_lb=f_lb, f_ub=f_ub, F=F,
# Vv=np.array(Vv), Mm=np.array(Mm), neg_c_max=neg_c_max,
zero=np.zeros((n,1)))
