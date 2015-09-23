/*
 *   -- MAGMA_minproduct (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions mixed zc -> ds
 */

#ifndef _MAGMA_minproduct_ZC_H_
#define _MAGMA_minproduct_ZC_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_minproduct_int_t magma_minproduct_zcgesv_gpu(  char trans, magma_minproduct_int_t N, magma_minproduct_int_t NRHS, 
                   cuDoubleComplex *dA, magma_minproduct_int_t ldda, 
                   magma_minproduct_int_t *IPIV, magma_minproduct_int_t *dIPIV,
                   cuDoubleComplex *dB, magma_minproduct_int_t lddb, 
                   cuDoubleComplex *dX, magma_minproduct_int_t lddx, 
                   cuDoubleComplex *dworkd, cuFloatComplex *dworks,
                   magma_minproduct_int_t *iter, magma_minproduct_int_t *info);

magma_minproduct_int_t magma_minproduct_zcgetrs_gpu( char trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
                   cuFloatComplex  *dA, magma_minproduct_int_t ldda,
                               magma_minproduct_int_t *ipiv, 
                   cuDoubleComplex *dB, magma_minproduct_int_t lddb,
                   cuDoubleComplex *dX, magma_minproduct_int_t lddx,
                               cuFloatComplex  *dSX, 
                   magma_minproduct_int_t *info );

magma_minproduct_int_t magma_minproduct_zcposv_gpu( char uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
                              cuDoubleComplex *dA, magma_minproduct_int_t ldda, 
                              cuDoubleComplex *dB, magma_minproduct_int_t lddb, 
                              cuDoubleComplex *dX, magma_minproduct_int_t lddx, 
                              cuDoubleComplex *dworkd, cuFloatComplex *dworks,
                              magma_minproduct_int_t *iter, magma_minproduct_int_t *info);

magma_minproduct_int_t magma_minproduct_zcgeqrsv_gpu(magma_minproduct_int_t M, magma_minproduct_int_t N, magma_minproduct_int_t NRHS, 
                   cuDoubleComplex *dA,  magma_minproduct_int_t ldda, 
                   cuDoubleComplex *dB,  magma_minproduct_int_t lddb, 
                   cuDoubleComplex *dX,  magma_minproduct_int_t lddx,
                   magma_minproduct_int_t *iter,    magma_minproduct_int_t *info);
  

#ifdef __cplusplus
}
#endif

#endif /* _MAGMA_minproduct_Z_H_ */
