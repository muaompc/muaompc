#ifndef {PREFIX}_CONST_H
#define {PREFIX}_CONST_H

#ifdef __cplusplus
extern "C" {{
#endif

#include "{prefix}_base.h"

{extra_defines}

enum {{ 
{PREFIX}_OPTVAR_VECLEN = {OPTVAR_VECLEN},  /**< Number of system inputs. */
{PREFIX}_MXCONSTRS_VECLEN = {MXCONSTRS_VECLEN}, /**< Number of mixed stage constraints. */
{PREFIX}_OPTVAR_SEQLEN = {OPTVAR_SEQLEN},  /**< Horizon times number of inputs. */
{PREFIX}_MXCONSTRS_SEQLEN = {MXCONSTRS_SEQLEN}  /**< Horizon times number of mixed constrained
plus the number of end state constraints. */
}}; 

#ifdef __cplusplus
}}
#endif

#endif /* {PREFIX}_CONST_H */
