#ifndef {PREFIX}_ALMDYNMEM_H
#define {PREFIX}_ALMDYNMEM_H

#include "{prefix}dynmem.h"
#include "{prefix}cvp.h"
#include "{prefix}alm.h"

extern struct {prefix}_alm *{prefix}_alm_allocate_solver(void);

extern {prefix}_dynmem_error_t {prefix}_alm_setup_solver(
                struct {prefix}_alm *alm,
                struct {prefix}_cvp_prb *prb, char *fname);
extern struct {prefix}_alm *{prefix}_alm_free_solver(struct {prefix}_alm *alm);

#endif /* {PREFIX}_ALMDYNMEM_H */
