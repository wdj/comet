/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_zc.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproduct_DS_H
#define MAGMA_minproduct_DS_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_minproduct_int_t
magma_minproduct_dsgesv_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *ipiv,
    magma_minproductInt_ptr dipiv,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDouble_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductDouble_ptr dworkd, magma_minproductFloat_ptr dworks,
    magma_minproduct_int_t *iter,
    magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_dsgetrs_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr  dA, magma_minproduct_int_t ldda,
    magma_minproductInt_ptr        dipiv,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDouble_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductFloat_ptr dSX,
    magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_dsposv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDouble_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductDouble_ptr dworkd, magma_minproductFloat_ptr dworks,
    magma_minproduct_int_t *iter,
    magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_dsgeqrsv_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDouble_ptr dX, magma_minproduct_int_t lddx,
    magma_minproduct_int_t *iter,
    magma_minproduct_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproduct_DS_H */
