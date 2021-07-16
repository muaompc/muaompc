#ifndef {PREFIX}CTLDYNMEM_H
#define {PREFIX}CTLDYNMEM_H

#include "{prefix}{solver}dynmem.h"
#include "{prefix}{former}dynmem.h"
#include "{prefix}ctl.h"

extern struct {prefix}_ctl *{prefix}_ctl_allocate_ctl(void);
extern {prefix}_dynmem_error_t {prefix}_ctl_setup_ctl(struct {prefix}_ctl *ctl, char *fname);
extern struct {prefix}_ctl *{prefix}_ctl_free_ctl(struct {prefix}_ctl *ctl);

#endif /* {PREFIX}CTLDYNMEM_H */
