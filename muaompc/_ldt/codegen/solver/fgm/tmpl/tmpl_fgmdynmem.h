#ifndef {PREFIX}_FGMDYNMEM_H
#define {PREFIX}_FGMDYNMEM_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}dynmem.h"
#include "{prefix}cvp.h"
#include "{prefix}fgm.h"

extern struct {prefix}_fgm *{prefix}_fgm_allocate_solver(void);
extern struct {prefix}_fgm *{prefix}_fgm_free_solver(struct {prefix}_fgm *fgm);
extern {prefix}_dynmem_error_t {prefix}_fgm_setup_solver(
                struct {prefix}_fgm *fgm,
                struct {prefix}_cvp_prb *prb, char *fname);

#ifdef __cplusplus
}}
#endif

#endif /* {PREFIX}_FGMDYNMEM_H */
