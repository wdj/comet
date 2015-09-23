#ifndef ZFILL_H
#define ZFILL_H

#include "magma_minproduct.h"

void zfill_matrix(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductDoubleComplex *A, magma_minproduct_int_t lda );

void zfill_rhs(
    magma_minproduct_int_t m, magma_minproduct_int_t nrhs, magma_minproductDoubleComplex *X, magma_minproduct_int_t ldx );

void zfill_matrix_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductDoubleComplex *dA, magma_minproduct_int_t ldda );

void zfill_rhs_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t nrhs, magma_minproductDoubleComplex *dX, magma_minproduct_int_t lddx );

#endif        //  #ifndef ZFILL_H
