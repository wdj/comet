/*
 *   -- MAGMA_tally2 (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions mixed zc -> ds
 */

#ifndef _MAGMA_tally2BLAS_ZC_H_
#define _MAGMA_tally2BLAS_ZC_H_

#include "cublas.h"
#include "cuda.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magma_tally2blas_zcaxpycp(cuFloatComplex *, cuDoubleComplex *, magma_tally2_int_t, cuDoubleComplex *, cuDoubleComplex *);
void magma_tally2blas_zaxpycp(cuDoubleComplex *, cuDoubleComplex *, magma_tally2_int_t, cuDoubleComplex *);
void magma_tally2blas_zclaswp(magma_tally2_int_t, cuDoubleComplex *, magma_tally2_int_t, cuFloatComplex *, magma_tally2_int_t, magma_tally2_int_t *, magma_tally2_int_t);
void magma_tally2blas_zlag2c(magma_tally2_int_t M, magma_tally2_int_t N, const cuDoubleComplex *A, magma_tally2_int_t lda,  cuFloatComplex *SA, magma_tally2_int_t ldsa, magma_tally2_int_t *info);

void magma_tally2blas_clag2z(magma_tally2_int_t M, magma_tally2_int_t N, 
                      cuFloatComplex  *SA, magma_tally2_int_t ldsa, 
                      cuDoubleComplex *A,  magma_tally2_int_t lda, 
                      magma_tally2_int_t *info);
void magma_tally2blas_zlat2c(char uplo, magma_tally2_int_t n, 
                      cuDoubleComplex *A,  magma_tally2_int_t lda, 
                      cuFloatComplex  *SA, magma_tally2_int_t ldsa, 
                      magma_tally2_int_t *info);

#ifdef __cplusplus
}
#endif

#endif
