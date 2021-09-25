#ifndef {PREFIX}_PBMDYNMEM_H
#define {PREFIX}_PBMDYNMEM_H

#ifdef __cplusplus
extern "C" {{
#endif

/* TODO implement socp.h */
//#include "../{prefix}socp.h"
#include "{prefix}dynmem.h"
#include "{prefix}cvp.h"
#include "{prefix}pbm.h"

/* External function declarations */

extern struct {prefix}_pbm *{prefix}_pbm_allocate_solver(void);

/* TODO implement socp.h */
/*
extern {prefix}_dynmem_error_t {prefix}_pbm_setup_solver(
                struct {prefix}_pbm *pbm, struct {prefix}_socp_prb *prb,
                char *fname); 
*/
extern {prefix}_dynmem_error_t {prefix}_pbm_setup_solver(
                struct {prefix}_pbm *pbm, struct {prefix}_cvp_prb *prb,
                char *fname);

#ifdef __cplusplus
}}
#endif

#endif /* {PREFIX}_PBMDYNMEM_H */
