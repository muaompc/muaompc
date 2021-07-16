
cdef extern from "../include/{prefix}cvp.h":

    cdef enum:
{par_inds}
        {PREFIX}_PAR_NUM

    cdef struct {prefix}_term:
        int rows
        int cols
        double *data

    cdef struct {prefix}_cvp_prb_socc:
        {prefix}_term *Wm;
        {prefix}_term *wn;
        {prefix}_term *wvT;
        {prefix}_term *ws;

    cdef struct {prefix}_cvp_prb:
{prb_terms}

    cdef struct {prefix}_cvp:
        {prefix}_term *par[{PREFIX}_PAR_NUM]
        {prefix}_cvp_prb *prb

    void {prefix}_cvp_form_problem({prefix}_cvp *cvp)

cdef extern from "../include/{prefix}cvpdynmem.h":

    int {prefix}_cvp_setup_former({prefix}_cvp *cvp, char *fname)
    {prefix}_cvp *{prefix}_cvp_allocate_former()

cdef _py2c_pardata({prefix}_cvp *ccvp, pardata)
cdef _c2py_prb({prefix}_cvp_prb *cprb)

