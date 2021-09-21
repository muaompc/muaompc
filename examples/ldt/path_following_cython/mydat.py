import numpy as np
from scipy.signal import cont2discrete as c2d
# we start with a simple mass, that moves in 1-dimension. Position given by x_1
# F = a; F = u; a = \dot(x)_2
Axc = np.array([[0, 1], [0, 0]])
Buc = np.array([[0], [1.]])
nx = Axc.shape[1]
mu = Buc.shape[1]
# The parameter dynamics: double integrator  \dot(z) = Azc*z + Bzc*v,  \theta = z1
Azc = np.array([[0, 1], [0, 0]])
Bvc = np.array([[0], [1.]])
nz = Azc.shape[1]
mv = Bvc.shape[1]

# We have thus the optimization variable, r.T = (u, v)
Bxc = np.concatenate((Buc, np.zeros([nx, mv])), axis=1)
Bzc = np.concatenate((np.zeros([nz, mu]), Bvc), axis=1)

# The error matrices
Ex = np.array([[1., 0, 0, 0], [0, 1, 0, 0]])
Ez = np.array([[0., 0, 1, 0], [0, 0, 0, 1]])
ne = Ex.shape[0]

# weighting matrices
Q = np.diag([1e3, 1.])
R = np.diag([1., 1.])
P = Q

# discretization
dt = 0.1
(Ax, Bx, C, D, dt) = c2d((Axc, Bxc, np.eye(nx), np.zeros(Bxc.shape)), dt)
(Az, Bz, C, D, dt) = c2d((Azc, Bzc, np.eye(nz), np.zeros(Bzc.shape)), dt)
# input constraints
# constraint on v must have 0 in its interior
r_lb = np.array([[-.1, -.01]]).T
r_ub = np.array([[.1, .1]]).T
# state constraints
# only constraint the path parameter, and it's speed (i.e. in z)
z_lb = np.array([[0., 0.]]).T
z_ub = np.array([[10, 10.]]).T
Kz = np.array([[0, 0, 1, 0], [0, 0, 0, 1]])
# dimensions
N = 10  # horizon length
n = nx  # muaompc only supports single character dimensions
m = Bxc.shape[1]  # same for Bx and Bz