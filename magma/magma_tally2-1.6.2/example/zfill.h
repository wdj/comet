#ifndef ZFILL_H
#define ZFILL_H

#include "magma_tally2.h"

void zfill_matrix(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2DoubleComplex *A, magma_tally2_int_t lda );

void zfill_rhs(
    magma_tally2_int_t m, magma_tally2_int_t nrhs, magma_tally2DoubleComplex *X, magma_tally2_int_t ldx );

void zfill_matrix_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2DoubleComplex *dA, magma_tally2_int_t ldda );

void zfill_rhs_gpu(
    magma_tally2_int_t m, magma_tally2_int_t nrhs, magma_tally2DoubleComplex *dX, magma_tally2_int_t lddx );

#endif        //  #ifndef ZFILL_H
