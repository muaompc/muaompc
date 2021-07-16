#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */

#include "cjson.h"

#include "{prefix}dynmem.h"

cJSON *{prefix}_dynmem_get_data(char *fname) {{
    cJSON *data;
    char *fdata;
    long len, bytesread;
    FILE *f=fopen(fname,"rb");
    if (NULL == f) {{
        printf("ERROR: could not open file %s \n", fname);
        return NULL;
    }}

    fseek(f,0,SEEK_END);
    len=ftell(f);
    fseek(f,0,SEEK_SET);
    fdata=(char*)malloc(len+1);
    if (NULL == fdata) {{
        printf("ERROR: could not allocate memory for file %s \n", fname);
        return NULL;
    }}

    bytesread = (long)fread(fdata, 1, len, f);
    if (bytesread != len) {{
        printf("ERROR: could not read file %s \n", fname);
        return NULL;
    }}
    fclose(f);

    data = cJSON_Parse(fdata);
    free(fdata);
    if (NULL == data) {{
        printf("ERROR: could not parse file %s \n", fname);
        return NULL;
    }}

    return data;
}}

void {prefix}_free(void *p) {{
  if (p != NULL) {{
    free(p);
    p = NULL;
  }}
  return;
}}
