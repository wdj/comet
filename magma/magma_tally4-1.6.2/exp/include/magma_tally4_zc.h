/*
 *   -- MAGMA_tally4 (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions mixed zc -> ds
 */

#ifndef _MAGMA_tally4_ZC_H_
#define _MAGMA_tally4_ZC_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_tally4_int_t magma_tally4_zcgesv_gpu(  char trans, magma_tally4_int_t N, magma_tally4_int_t NRHS, 
                   cuDoubleComplex *dA, magma_tally4_int_t ldda, 
                   magma_tally4_int_t *IPIV, magma_tally4_int_t *dIPIV,
                   cuDoubleComplex *dB, magma_tally4_int_t lddb, 
                   cuDoubleComplex *dX, magma_tally4_int_t lddx, 
                   cuDoubleComplex *dworkd, cuFloatComplex *dworks,
                   magma_tally4_int_t *iter, magma_tally4_int_t *info);

magma_tally4_int_t magma_tally4_zcgetrs_gpu( char trans, magma_tally4_int_t n, magma_tally4_int_t nrhs, 
                   cuFloatComplex  *dA, magma_tally4_int_t ldda,
                               magma_tally4_int_t *ipiv, 
                   cuDoubleComplex *dB, magma_tally4_int_t lddb,
                   cuDoubleComplex *dX, magma_tally4_int_t lddx,
                               cuFloatComplex  *dSX, 
                   magma_tally4_int_t *info );

magma_tally4_int_t magma_tally4_zcposv_gpu( char uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs, 
                              cuDoubleComplex *dA, magma_tally4_int_t ldda, 
                              cuDoubleComplex *dB, magma_tally4_int_t lddb, 
                              cuDoubleComplex *dX, magma_tally4_int_t lddx, 
                              cuDoubleComplex *dworkd, cuFloatComplex *dworks,
                              magma_tally4_int_t *iter, magma_tally4_int_t *info);

magma_tally4_int_t magma_tally4_zcgeqrsv_gpu(magma_tally4_int_t M, magma_tally4_int_t N, magma_tally4_int_t NRHS, 
                   cuDoubleComplex *dA,  magma_tally4_int_t ldda, 
                   cuDoubleComplex *dB,  magma_tally4_int_t lddb, 
                   cuDoubleComplex *dX,  magma_tally4_int_t lddx,
                   magma_tally4_int_t *iter,    magma_tally4_int_t *info);
  

#ifdef __cplusplus
}
#endif

#endif /* _MAGMA_tally4_Z_H_ */
