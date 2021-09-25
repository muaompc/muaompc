#ifndef {PREFIX}{DATANAME}CTLDATA_H
#define {PREFIX}{DATANAME}CTLDATA_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}ctl.h"

extern void {prefix}{dataname}_ctl_setup_ctl(
		struct {prefix}_ctl *ctl);

#ifdef __cplusplus
}}
#endif

#endif
