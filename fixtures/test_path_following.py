import numpy as np
# we start with a simple mass, that moves in 1-dimension. Position given by x_1
# F = a; F = u; a = \dot(x)_2
Axc = np.array([[0., 1], [0, 1]])
Buc = np.array([[0.], [1]])
nx = Axc.shape[1]
mu = Buc.shape[1]
# The parameter dynamics: double integrator  \dot(z) = Azc*z + Bzc*v,  \theta = z1
Azc = np.array([[0., 2], [0, 1]])
Bvc = np.array([[0.], [1]])
nz = Azc.shape[1]
mv = Bvc.shape[1]

# We have thus the optimization variable, r.T = (u, v)
Bxc = np.concatenate((Buc, np.zeros([nx, mv])), axis=1)
Bzc = np.concatenate((np.zeros([nz, mu]), Bvc), axis=1)

nx = Axc.shape[1]  # number of states in x
nz = Azc.shape[1]  # number of states in z

# weighting matrices
Q = np.diag([1., 1.])
R = np.diag([1., 1.])
P = np.diag([1., 1.])

(Ax, Bx, Az, Bz) = (Axc, Bxc, Azc, Bzc)
# input constraints
# constraint on v must have 0 in its interior
r_lb = np.array([[-0.1, -0.1]]).T
r_ub = np.array([[0.1, 0.1]]).T
# state constraints
# only constraint the path parameter, and it's speed (i.e. in z)
z_lb = np.array([[0., 0.]]).T
z_ub = np.array([[10, 10.]]).T
Kz = np.array([[0, 0, 1, 0], [0, 0, 0, 1]])
# dimensions
N = 3  # horizon length
n = nx  # muaompc only supports single character dimensions
m = Bx.shape[1]