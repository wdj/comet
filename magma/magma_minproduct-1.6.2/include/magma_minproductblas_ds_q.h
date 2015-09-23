/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproductblas_zc_q.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproductBLAS_DS_Q_H
#define MAGMA_minproductBLAS_DS_Q_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magma_minproductblas_dsaxpycp_q(
    magma_minproduct_int_t m,
    magma_minproductFloat_ptr  r,
    magma_minproductDouble_ptr x,
    magma_minproductDouble_const_ptr b,
    magma_minproductDouble_ptr w,
    magma_minproduct_queue_t queue );

void magma_minproductblas_daxpycp_q(
    magma_minproduct_int_t m,
    magma_minproductDouble_ptr r,
    magma_minproductDouble_ptr x,
    magma_minproductDouble_const_ptr b,
    magma_minproduct_queue_t queue  );

void magma_minproductblas_dslaswp_q(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr A, magma_minproduct_int_t lda,
    magma_minproductFloat_ptr SA, magma_minproduct_int_t m,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t incx,
    magma_minproduct_queue_t queue );

void magma_minproductblas_dlag2s_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr A,  magma_minproduct_int_t lda,
    magma_minproductFloat_ptr       SA, magma_minproduct_int_t ldsa,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void magma_minproductblas_slag2d_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr SA, magma_minproduct_int_t ldsa,
    magma_minproductDouble_ptr       A,  magma_minproduct_int_t lda,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void magma_minproductblas_dlat2s_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr A,  magma_minproduct_int_t lda,
    magma_minproductFloat_ptr       SA, magma_minproduct_int_t ldsa,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void magma_minproductblas_slat2d_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr SA, magma_minproduct_int_t ldsa,
    magma_minproductDouble_ptr       A,  magma_minproduct_int_t lda,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproductBLAS_DS_H */
