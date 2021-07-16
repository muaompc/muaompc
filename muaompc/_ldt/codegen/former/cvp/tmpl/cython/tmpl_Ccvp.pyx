import numpy as np

cimport {prefix}Ccvp as Ccvp

include "{prefix}Ccvp.pxi"


cdef _py2c_pardata(Ccvp.{prefix}_cvp *ccvp, pardata):
{pardata_py2c}

cdef _c2py_prb(Ccvp.{prefix}_cvp_prb *cprb):
    prb = dict()
{prb_c2py}
    return prb

cdef _py2c_term(Ccvp.{prefix}_term *term, pydata):
    cdef int i
    cdef int j
    if pydata.shape != (term.rows, term.cols):
        print('Error: incompatible shapes')
        return

    for i in range(term.rows):
        for j in range(term.cols):
            term.data[i*term.cols + j] = pydata[i,j]
    return

cdef _c2py_term(Ccvp.{prefix}_term *term):
    a = np.zeros((term.rows, term.cols))

    cdef int i
    cdef int j
    for i in range(term.rows):
        for j in range(term.cols):
            a[i,j] = term.data[i*term.cols + j]
    return a

IF CVP_PRB_SOCC == 1:
  cdef _c2py_socc(Ccvp.{prefix}_cvp_prb *cprb):
      cdef int k
      socc = []
      for k in range(cprb.socc_num[0]):
          terms = dict()
          terms['Wm'] = _c2py_term(cprb.socc[k].Wm)
          terms['wvT'] = _c2py_term(cprb.socc[k].wvT)
          terms['wn'] = _c2py_term(cprb.socc[k].wn)
          terms['ws'] = _c2py_term(cprb.socc[k].ws)
          socc.append(terms)
      return socc
