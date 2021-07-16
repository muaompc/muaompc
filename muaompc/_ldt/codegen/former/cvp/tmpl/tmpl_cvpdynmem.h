#ifndef {PREFIX}_CVPDYNMEM_H
#define {PREFIX}_CVPDYNMEM_H

#include "{prefix}dynmem.h"
#include "{prefix}cvp.h"

extern {prefix}_dynmem_error_t {prefix}_cvp_setup_former(struct {prefix}_cvp *cvp, char *fname);
extern {prefix}_dynmem_error_t {prefix}_cvp_setup_parameters(struct {prefix}_cvp_parameters *parameters,
                struct {prefix}_cvp *cvp);
extern struct {prefix}_cvp *{prefix}_cvp_allocate_former(void);
extern struct {prefix}_cvp_parameters *{prefix}_cvp_allocate_parameters(void);
extern struct {prefix}_cvp *{prefix}_cvp_free_former(struct {prefix}_cvp *cvp);
struct {prefix}_cvp_parameters *{prefix}_cvp_free_parameters(struct {prefix}_cvp_parameters *p);

#endif /* {PREFIX}_CVPDYNMEM_H */
