/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_tally2BLAS_ZC_Q_H
#define MAGMA_tally2BLAS_ZC_Q_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magma_tally2blas_zcaxpycp_q(
    magma_tally2_int_t m,
    magma_tally2FloatComplex_ptr  r,
    magma_tally2DoubleComplex_ptr x,
    magma_tally2DoubleComplex_const_ptr b,
    magma_tally2DoubleComplex_ptr w,
    magma_tally2_queue_t queue );

void magma_tally2blas_zaxpycp_q(
    magma_tally2_int_t m,
    magma_tally2DoubleComplex_ptr r,
    magma_tally2DoubleComplex_ptr x,
    magma_tally2DoubleComplex_const_ptr b,
    magma_tally2_queue_t queue  );

void magma_tally2blas_zclaswp_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr A, magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr SA, magma_tally2_int_t m,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t incx,
    magma_tally2_queue_t queue );

void magma_tally2blas_zlag2c_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr A,  magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr       SA, magma_tally2_int_t ldsa,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void magma_tally2blas_clag2z_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr SA, magma_tally2_int_t ldsa,
    magma_tally2DoubleComplex_ptr       A,  magma_tally2_int_t lda,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void magma_tally2blas_zlat2c_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr A,  magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr       SA, magma_tally2_int_t ldsa,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void magma_tally2blas_clat2z_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr SA, magma_tally2_int_t ldsa,
    magma_tally2DoubleComplex_ptr       A,  magma_tally2_int_t lda,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally2BLAS_ZC_H */
