#ifndef ZFILL_H
#define ZFILL_H

#include "magma_tally3.h"

void zfill_matrix(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3DoubleComplex *A, magma_tally3_int_t lda );

void zfill_rhs(
    magma_tally3_int_t m, magma_tally3_int_t nrhs, magma_tally3DoubleComplex *X, magma_tally3_int_t ldx );

void zfill_matrix_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3DoubleComplex *dA, magma_tally3_int_t ldda );

void zfill_rhs_gpu(
    magma_tally3_int_t m, magma_tally3_int_t nrhs, magma_tally3DoubleComplex *dX, magma_tally3_int_t lddx );

#endif        //  #ifndef ZFILL_H
