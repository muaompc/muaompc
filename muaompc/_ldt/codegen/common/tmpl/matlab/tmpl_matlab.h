#ifndef MPC_{PREFIX}MATLAB_H
#define MPC_{PREFIX}MATLAB_H

#include "mex.h"

#include "arithmetic.h"

uint32_t {prefix}_matlab_cell_index(
    uint32_t row, 
    uint32_t col, 
    uint32_t cols
);

void {prefix}_matlab_move_matrix(
    real_t* destination, 
    const real_t* source, 
    uint32_t rows,
    uint32_t cols
);

void {prefix}_matlab_add_matrix(
    mxArray* structure,
    const char* name,
    const real_t* data, 
    uint32_t rows,
    uint32_t cols
);

void {prefix}_matlab_add_scalar(
    mxArray* structure, 
    const char* name, 
    real_t data
);

void {prefix}_matlab_read_matrix(
    const mxArray* structure,
    const char* name,
    real_t* data,
    uint32_t rows,
    uint32_t cols
);

void {prefix}_matlab_read_scalar(
    const mxArray* structure,
    const char* name,
    real_t* data
);

void {prefix}_matlab_read_uint(
    const mxArray* structure,
    const char* name,
    uint32_t* data
);

#endif
