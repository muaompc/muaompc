import numpy as np
from scipy.signal import cont2discrete as c2d
# weighting matrices
Q = np.diag([1014.7, 3.2407, 5674.8, 0.3695, 471.75])
R = np.array([[472.]])
P = Q
# system matrices (continuos time)
Ac = np.array([[-1.2822, 0, 0.98, 0], [0, 0, 1, 0], [-5.4293, 0, -1.8366, 0], [-128.2, 128.2, 0, 0]])
Bc = np.array([[0.3], [0], [-17], [0]])
Cc = np.array([[0, 1, 0, 0], [0, 0, 0, 1], [-128.2, 128.2, 0, 0]])
Dc = np.zeros((3,1))
# discretization
dt = 0.5
(A, B, C, D, dt) = c2d((Ac, Bc, Cc, Dc), dt)
# The system is extended with a 5th state, to account for slew rate constraints
# Extend A from a 4x4 matrix, to a 5x5 matrix
# Add first a row of zeros:
A = np.concatenate((A, np.zeros((1,4))))
# then a column of zeros
A = np.concatenate((A, np.zeros((5,1))), axis=1)
# Extend B from a 4x1 column vector, to a 5x1 column vector
# Add an new row at the bottom, with the element = 1:
B = np.concatenate((B, np.ones((1,1))))
# input constraints
u_lb = np.array([[-0.262]])
u_ub = np.array([[0.262]])
# state constraints
e_lb = np.array([[-0.349, -30, -0.25]]).T
e_ub = -1*e_lb
Kx = np.array([[0, 1, 0, 0, 0], [-128.2, 128.2, 0, 0, 0], [0., 0., 0., 0., -1.]])
Ku = np.array([[0, 0, 1]]).T
# terminal state constraints
f_lb = e_lb
f_ub = e_ub
Kf = Kx 
# dimensions
N = 10  # horizon length
n = 5  # number of states
m = 1  # number of inputs