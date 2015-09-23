/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_minproductBLAS_ZC_Q_H
#define MAGMA_minproductBLAS_ZC_Q_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magma_minproductblas_zcaxpycp_q(
    magma_minproduct_int_t m,
    magma_minproductFloatComplex_ptr  r,
    magma_minproductDoubleComplex_ptr x,
    magma_minproductDoubleComplex_const_ptr b,
    magma_minproductDoubleComplex_ptr w,
    magma_minproduct_queue_t queue );

void magma_minproductblas_zaxpycp_q(
    magma_minproduct_int_t m,
    magma_minproductDoubleComplex_ptr r,
    magma_minproductDoubleComplex_ptr x,
    magma_minproductDoubleComplex_const_ptr b,
    magma_minproduct_queue_t queue  );

void magma_minproductblas_zclaswp_q(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr SA, magma_minproduct_int_t m,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t incx,
    magma_minproduct_queue_t queue );

void magma_minproductblas_zlag2c_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr A,  magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr       SA, magma_minproduct_int_t ldsa,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void magma_minproductblas_clag2z_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr SA, magma_minproduct_int_t ldsa,
    magma_minproductDoubleComplex_ptr       A,  magma_minproduct_int_t lda,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void magma_minproductblas_zlat2c_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr A,  magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr       SA, magma_minproduct_int_t ldsa,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void magma_minproductblas_clat2z_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr SA, magma_minproduct_int_t ldsa,
    magma_minproductDoubleComplex_ptr       A,  magma_minproduct_int_t lda,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproductBLAS_ZC_H */
