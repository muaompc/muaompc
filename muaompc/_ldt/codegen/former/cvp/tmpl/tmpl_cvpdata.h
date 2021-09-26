#ifndef {PREFIX}{DATANAME}CVPDATA_H
#define {PREFIX}{DATANAME}CVPDATA_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}cvp.h"

extern void {prefix}{dataname}_cvp_setup_former(struct {prefix}_cvp *cvp);

extern void {prefix}{dataname}_cvp_setup_parameters(struct {prefix}_cvp_parameters *parameters, 
		struct {prefix}_cvp *cvp);
#ifdef __cplusplus
}}
#endif

#endif
