/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_tally4BLAS_ZC_H
#define MAGMA_tally4BLAS_ZC_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void
magma_tally4blas_zcaxpycp(
    magma_tally4_int_t m,
    magma_tally4FloatComplex_ptr  r,
    magma_tally4DoubleComplex_ptr x,
    magma_tally4DoubleComplex_const_ptr b,
    magma_tally4DoubleComplex_ptr w );

void
magma_tally4blas_zaxpycp(
    magma_tally4_int_t m,
    magma_tally4DoubleComplex_ptr r,
    magma_tally4DoubleComplex_ptr x,
    magma_tally4DoubleComplex_const_ptr b );

// TODO add ldsa
void
magma_tally4blas_zclaswp(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr  A, magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr  SA,
    magma_tally4_int_t m, const magma_tally4_int_t *ipiv, magma_tally4_int_t incx );

void
magma_tally4blas_zlag2c(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr  A, magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr        SA, magma_tally4_int_t ldsa,
    magma_tally4_int_t *info );

void
magma_tally4blas_clag2z(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr  SA, magma_tally4_int_t ldsa,
    magma_tally4DoubleComplex_ptr        A, magma_tally4_int_t lda,
    magma_tally4_int_t *info );

void
magma_tally4blas_zlat2c(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr  A, magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr        SA, magma_tally4_int_t ldsa,
    magma_tally4_int_t *info );

void
magma_tally4blas_clat2z(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr  SA, magma_tally4_int_t ldsa,
    magma_tally4DoubleComplex_ptr        A, magma_tally4_int_t lda,
    magma_tally4_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally4BLAS_ZC_H */
