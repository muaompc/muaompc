#ifndef STOCMPC_DYNMEM_H
#define STOCMPC_DYNMEM_H

#include "cjson.h"

typedef enum stocmpc_dynmem_error {
    STOCMPC_DYNMEM_OK = 0,
    STOCMPC_DYNMEM_FAIL = 1
} stocmpc_dynmem_error_t;

cJSON *stocmpc_dynmem_get_data(char *fname);
void stocmpc_free(void *p);

#endif
