/**
 *
 * @file transpose.h
 *
 *  MAGMA_tally3 (version 1.6.1) --
 *  Univ. of Tennessee, Knoxville
 *  Univ. of California, Berkeley
 *  Univ. of Colorado, Denver
 *  @date January 2015
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date January 2015
 *
 * Macro to transpose matrices before and after computation
 * in LU kernels
 *
 **/

#ifndef _MAGMA_tally3_TRANSPOSE_H_
#define _MAGMA_tally3_TRANSPOSE_H_

#define magma_tally3blas_sgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally3blas_sinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(float), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magma_tally3blas_stranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magma_tally3blas_sgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally3blas_sinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magma_tally3blas_stranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magma_tally3blas_dgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally3blas_dinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(double), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magma_tally3blas_dtranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magma_tally3blas_dgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally3blas_dinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magma_tally3blas_dtranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magma_tally3blas_cgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally3blas_cinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(cuFloatComplex), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magma_tally3blas_ctranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magma_tally3blas_cgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally3blas_cinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magma_tally3blas_ctranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magma_tally3blas_zgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally3blas_zinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(cuDoubleComplex), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magma_tally3blas_ztranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magma_tally3blas_zgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally3blas_zinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magma_tally3blas_ztranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#endif /* _MAGMA_tally3_TRANSPOSE_H_ */
