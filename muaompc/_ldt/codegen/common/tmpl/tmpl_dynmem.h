#ifndef {PREFIX}_DYNMEM_H
#define {PREFIX}_DYNMEM_H

#include "cjson.h"

typedef enum {prefix}_dynmem_error {{
    {PREFIX}_DYNMEM_OK = 0,
    {PREFIX}_DYNMEM_FAIL = 1
}} {prefix}_dynmem_error_t;

cJSON *{prefix}_dynmem_get_data(char *fname);
void {prefix}_free(void *p);

#endif
