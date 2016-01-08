#ifndef ZFILL_H
#define ZFILL_H

#include "magma_tally4.h"

void zfill_matrix(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4DoubleComplex *A, magma_tally4_int_t lda );

void zfill_rhs(
    magma_tally4_int_t m, magma_tally4_int_t nrhs, magma_tally4DoubleComplex *X, magma_tally4_int_t ldx );

void zfill_matrix_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4DoubleComplex *dA, magma_tally4_int_t ldda );

void zfill_rhs_gpu(
    magma_tally4_int_t m, magma_tally4_int_t nrhs, magma_tally4DoubleComplex *dX, magma_tally4_int_t lddx );

#endif        //  #ifndef ZFILL_H
