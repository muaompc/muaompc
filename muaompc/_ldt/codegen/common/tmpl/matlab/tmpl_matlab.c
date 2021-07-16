#include "{prefix}matlab.h"

uint32_t {prefix}_matlab_cell_index(
    uint32_t row,
    uint32_t col,
    uint32_t cols
) {{
    return row * cols + col;
}}

void {prefix}_matlab_move_matrix(
    real_t* destination,
    const real_t* source,
    uint32_t rows,
    uint32_t cols
) {{
    unsigned int row = 0;
    unsigned int col = 0;

    /* loops are neccessary in order to transpose matrices */
    for (row = 0; row < rows; row++) {{
        for (col = 0; col < cols; col++) {{
            destination[{prefix}_matlab_cell_index(col, row, rows)] =
                (double)source[{prefix}_matlab_cell_index(row, col, cols)];
        }}
        col = 0;
    }}
}}

void {prefix}_matlab_add_matrix(
    mxArray* structure,
    const char* name,
    const real_t* data,
    uint32_t rows,
    uint32_t cols
) {{
    mxArray* matrix = mxCreateDoubleMatrix(rows, cols, 0);
    {prefix}_matlab_move_matrix(mxGetPr(matrix), data, rows, cols);
    mxAddField(structure, name);
    mxSetField(structure, 0, name, matrix);
}}

void {prefix}_matlab_add_scalar(
    mxArray* structure,
    const char* name,
    real_t data
) {{
    mxArray* scalar = mxCreateDoubleScalar(data);
    mxAddField(structure, name);
    mxSetField(structure, 0, name, scalar);
}}

void {prefix}_matlab_read_matrix(
    const mxArray* structure,
    const char* name,
    real_t* data,
    uint32_t rows,
    uint32_t cols
) {{
    mxArray* matrix = mxGetField(structure, 0, name);
    {prefix}_matlab_move_matrix(data, mxGetPr(matrix), rows, cols);
}}

void {prefix}_matlab_read_scalar(
    const mxArray* structure,
    const char* name,
    real_t* data
) {{
    mxArray* scalar = mxGetField(structure, 0, name);
    *data = mxGetScalar(scalar);
}}

void {prefix}_matlab_read_uint(
    const mxArray* structure,
    const char* name,
    uint32_t* data
) {{
    double temporary;
    {prefix}_matlab_read_scalar(structure, name, &temporary);
    *data = (uint32_t)temporary;
}}

