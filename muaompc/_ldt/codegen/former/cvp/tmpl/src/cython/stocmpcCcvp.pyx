import numpy as np

cimport stocmpcCcvp as Ccvp

cdef _py2c_pardata(Ccvp.stocmpc_cvp *ccvp, pardata):
    _py2c_term(ccvp.par[Ccvp.STOCMPC_X_K], pardata["x_k"])


cdef _c2py_qp(Ccvp.stocmpc_cvp_prb *cprb):
    prb = dict()
    prb['ws1'] = _c2py_term(cprb.ws1)
    prb['ws0'] = _c2py_term(cprb.ws0)
    prb['wn0'] = _c2py_term(cprb.wn0)
    prb['wn1'] = _c2py_term(cprb.wn1)
    prb['g'] = _c2py_term(cprb.g)
    prb['u_lb'] = _c2py_term(cprb.u_lb)
    prb['H'] = _c2py_term(cprb.H)
    prb['Wm1'] = _c2py_term(cprb.Wm1)
    prb['Wm0'] = _c2py_term(cprb.Wm0)
    prb['wvT1'] = _c2py_term(cprb.wvT1)
    prb['wvT0'] = _c2py_term(cprb.wvT0)
    prb['u_ub'] = _c2py_term(cprb.u_ub)

    return prb

cdef _py2c_term(Ccvp.stocmpc_term *term, pydata):
    cdef int i
    cdef int j
    if pydata.shape != (term.rows, term.cols):
        print('Error: incompatible shapes')
        return

    for i in range(term.rows):
        for j in range(term.cols):
            term.data[i*term.cols + j] = pydata[i,j]
    return

cdef _c2py_term(Ccvp.stocmpc_term *term):
    a = np.zeros((term.rows, term.cols))

    cdef int i
    cdef int j
    for i in range(term.rows):
        for j in range(term.cols):
            a[i,j] = term.data[i*term.cols + j]
    return a

