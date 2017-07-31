/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_tally2BLAS_ZC_H
#define MAGMA_tally2BLAS_ZC_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void
magma_tally2blas_zcaxpycp(
    magma_tally2_int_t m,
    magma_tally2FloatComplex_ptr  r,
    magma_tally2DoubleComplex_ptr x,
    magma_tally2DoubleComplex_const_ptr b,
    magma_tally2DoubleComplex_ptr w );

void
magma_tally2blas_zaxpycp(
    magma_tally2_int_t m,
    magma_tally2DoubleComplex_ptr r,
    magma_tally2DoubleComplex_ptr x,
    magma_tally2DoubleComplex_const_ptr b );

// TODO add ldsa
void
magma_tally2blas_zclaswp(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr  A, magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr  SA,
    magma_tally2_int_t m, const magma_tally2_int_t *ipiv, magma_tally2_int_t incx );

void
magma_tally2blas_zlag2c(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr  A, magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr        SA, magma_tally2_int_t ldsa,
    magma_tally2_int_t *info );

void
magma_tally2blas_clag2z(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr  SA, magma_tally2_int_t ldsa,
    magma_tally2DoubleComplex_ptr        A, magma_tally2_int_t lda,
    magma_tally2_int_t *info );

void
magma_tally2blas_zlat2c(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr  A, magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr        SA, magma_tally2_int_t ldsa,
    magma_tally2_int_t *info );

void
magma_tally2blas_clat2z(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr  SA, magma_tally2_int_t ldsa,
    magma_tally2DoubleComplex_ptr        A, magma_tally2_int_t lda,
    magma_tally2_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally2BLAS_ZC_H */
