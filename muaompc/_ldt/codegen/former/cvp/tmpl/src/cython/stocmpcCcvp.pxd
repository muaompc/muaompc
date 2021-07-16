cdef extern from "../include/stocmpccvp.h":

    cdef enum:
        STOCMPC_X_K = 0,

        STOCMPC_PAR_NUM

    cdef struct stocmpc_term:
        int rows
        int cols
        double *data

    cdef struct stocmpc_cvp_prb:
        stocmpc_term *ws1
        stocmpc_term *ws0
        stocmpc_term *wn0
        stocmpc_term *wn1
        stocmpc_term *g
        stocmpc_term *u_lb
        stocmpc_term *H
        stocmpc_term *Wm1
        stocmpc_term *Wm0
        stocmpc_term *wvT1
        stocmpc_term *wvT0
        stocmpc_term *u_ub


    cdef struct stocmpc_cvp:
        stocmpc_term *par[STOCMPC_PAR_NUM]
        stocmpc_cvp_prb *prb

    void stocmpc_cvp_form_problem(stocmpc_cvp *cvp)

cdef extern from "../include/stocmpccvpdynmem.h":

    int stocmpc_cvp_setup_former(stocmpc_cvp *cvp, char *fname)
    stocmpc_cvp *stocmpc_cvp_allocate_former()

cdef _py2c_pardata(stocmpc_cvp *ccvp, pardata)
cdef _c2py_qp(stocmpc_cvp_prb *cprb)

