#ifndef STOCMPC_CVPDYNMEM_H
#define STOCMPC_CVPDYNMEM_H

#include "stocmpcdynmem.h"
#include "stocmpccvp.h"

extern stocmpc_dynmem_error_t stocmpc_cvp_setup_former(struct stocmpc_cvp *cvp, char *fname);
extern stocmpc_dynmem_error_t stocmpc_cvp_setup_parameters(struct stocmpc_cvp_parameters *parameters,
                struct stocmpc_cvp *cvp);
extern struct stocmpc_cvp *stocmpc_cvp_allocate_former(void);
extern struct stocmpc_cvp_parameters *stocmpc_cvp_allocate_parameters(void);
extern struct stocmpc_cvp *stocmpc_cvp_free_former(struct stocmpc_cvp *cvp);
struct stocmpc_cvp_parameters *stocmpc_cvp_free_parameters(struct stocmpc_cvp_parameters *p);

#endif /* STOCMPC_CVPDYNMEM_H */
