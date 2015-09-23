/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_minproduct_ZC_H
#define MAGMA_minproduct_ZC_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_minproduct_int_t
magma_minproduct_zcgesv_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *ipiv,
    magma_minproductInt_ptr dipiv,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductDoubleComplex_ptr dworkd, magma_minproductFloatComplex_ptr dworks,
    magma_minproduct_int_t *iter,
    magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_zcgetrs_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr  dA, magma_minproduct_int_t ldda,
    magma_minproductInt_ptr        dipiv,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductFloatComplex_ptr dSX,
    magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_zcposv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductDoubleComplex_ptr dworkd, magma_minproductFloatComplex_ptr dworks,
    magma_minproduct_int_t *iter,
    magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_zcgeqrsv_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex_ptr dX, magma_minproduct_int_t lddx,
    magma_minproduct_int_t *iter,
    magma_minproduct_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproduct_ZC_H */
