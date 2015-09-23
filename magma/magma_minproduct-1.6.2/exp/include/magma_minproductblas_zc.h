/*
 *   -- MAGMA_minproduct (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions mixed zc -> ds
 */

#ifndef _MAGMA_minproductBLAS_ZC_H_
#define _MAGMA_minproductBLAS_ZC_H_

#include "cublas.h"
#include "cuda.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magma_minproductblas_zcaxpycp(cuFloatComplex *, cuDoubleComplex *, magma_minproduct_int_t, cuDoubleComplex *, cuDoubleComplex *);
void magma_minproductblas_zaxpycp(cuDoubleComplex *, cuDoubleComplex *, magma_minproduct_int_t, cuDoubleComplex *);
void magma_minproductblas_zclaswp(magma_minproduct_int_t, cuDoubleComplex *, magma_minproduct_int_t, cuFloatComplex *, magma_minproduct_int_t, magma_minproduct_int_t *, magma_minproduct_int_t);
void magma_minproductblas_zlag2c(magma_minproduct_int_t M, magma_minproduct_int_t N, const cuDoubleComplex *A, magma_minproduct_int_t lda,  cuFloatComplex *SA, magma_minproduct_int_t ldsa, magma_minproduct_int_t *info);

void magma_minproductblas_clag2z(magma_minproduct_int_t M, magma_minproduct_int_t N, 
                      cuFloatComplex  *SA, magma_minproduct_int_t ldsa, 
                      cuDoubleComplex *A,  magma_minproduct_int_t lda, 
                      magma_minproduct_int_t *info);
void magma_minproductblas_zlat2c(char uplo, magma_minproduct_int_t n, 
                      cuDoubleComplex *A,  magma_minproduct_int_t lda, 
                      cuFloatComplex  *SA, magma_minproduct_int_t ldsa, 
                      magma_minproduct_int_t *info);

#ifdef __cplusplus
}
#endif

#endif
