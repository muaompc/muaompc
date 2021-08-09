import numpy as np
# input constraints
u_lb = np.array([[-100]])
u_ub = np.array([[100]])
# weighting matrices
Q = np.array([[1, 0], [0, 1]])
R = np.array([[1]])
P = np.array([[2, 3], [3, 4]])
A = np.array([[ 1.,  0.0095162581964040437],
        [ 0.,  0.90483741803595941]])
B = np.array([[9.6748360719191485e-05],
        [0.019032516392808087]])
N = 2
mu = 256.
e_lb = np.array([[-1]])
e_ub = np.array([[1]])
f_lb = np.array([[-0.5]])
f_ub = np.array([[0.5]])
Kx = np.array([[3.5, 5.5]])
F = np.array([[0.2, 1.1]])
Ad = A
Bd = B
dt = 0.001
n, m = B.shape
neg_c_max = -0.3490
Mm = np.array([[1, 0]])
Vv = np.array([[1, 0], [0, 1]])
data = dict(A=A, B=B, P=P, Q=Q, R=R, Kx=Kx, N=N, n=n, m=m,
u_lb=u_lb, u_ub=u_ub,
e_lb=e_lb, e_ub=e_ub,
f_lb=f_lb, f_ub=f_ub, F=F,
mu=mu,
neg_c_max=neg_c_max, Mm=Mm, Vv=Vv
)
