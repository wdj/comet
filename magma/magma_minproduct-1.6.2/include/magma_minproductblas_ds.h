/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproductblas_zc.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproductBLAS_DS_H
#define MAGMA_minproductBLAS_DS_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void
magma_minproductblas_dsaxpycp(
    magma_minproduct_int_t m,
    magma_minproductFloat_ptr  r,
    magma_minproductDouble_ptr x,
    magma_minproductDouble_const_ptr b,
    magma_minproductDouble_ptr w );

void
magma_minproductblas_daxpycp(
    magma_minproduct_int_t m,
    magma_minproductDouble_ptr r,
    magma_minproductDouble_ptr x,
    magma_minproductDouble_const_ptr b );

// TODO add ldsa
void
magma_minproductblas_dslaswp(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr  A, magma_minproduct_int_t lda,
    magma_minproductFloat_ptr  SA,
    magma_minproduct_int_t m, const magma_minproduct_int_t *ipiv, magma_minproduct_int_t incx );

void
magma_minproductblas_dlag2s(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr  A, magma_minproduct_int_t lda,
    magma_minproductFloat_ptr        SA, magma_minproduct_int_t ldsa,
    magma_minproduct_int_t *info );

void
magma_minproductblas_slag2d(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr  SA, magma_minproduct_int_t ldsa,
    magma_minproductDouble_ptr        A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_dlat2s(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr  A, magma_minproduct_int_t lda,
    magma_minproductFloat_ptr        SA, magma_minproduct_int_t ldsa,
    magma_minproduct_int_t *info );

void
magma_minproductblas_slat2d(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr  SA, magma_minproduct_int_t ldsa,
    magma_minproductDouble_ptr        A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproductBLAS_DS_H */
