#ifndef {PREFIX}{DATANAME}CVPDATA_H
#define {PREFIX}{DATANAME}CVPDATA_H

#include "{prefix}cvp.h"

extern void {prefix}{dataname}_cvp_setup_former(struct {prefix}_cvp *cvp);

extern void {prefix}{dataname}_cvp_setup_parameters(struct {prefix}_cvp_parameters *parameters, 
		struct {prefix}_cvp *cvp);
#endif
