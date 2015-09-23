/**
 *
 * @file transpose.h
 *
 *  MAGMA_minproduct (version 1.6.1) --
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

#ifndef _MAGMA_minproduct_TRANSPOSE_H_
#define _MAGMA_minproduct_TRANSPOSE_H_

#define magma_minproductblas_sgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_sinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(float), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magma_minproductblas_stranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magma_minproductblas_sgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_sinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magma_minproductblas_stranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magma_minproductblas_dgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_dinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(double), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magma_minproductblas_dtranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magma_minproductblas_dgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_dinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magma_minproductblas_dtranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magma_minproductblas_cgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_cinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(cuFloatComplex), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magma_minproductblas_ctranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magma_minproductblas_cgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_cinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magma_minproductblas_ctranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#define magma_minproductblas_zgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_zinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    cublasStatus status = cublasAlloc( (m)*(n), sizeof(cuDoubleComplex), (void**)&dAT); \
    if (status != CUBLAS_STATUS_SUCCESS)                                \
      return -7;                                                        \
    magma_minproductblas_ztranspose2( dAT, ldda, dA, ldda, m, n );                 \
  }

#define magma_minproductblas_zgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_zinplace_transpose( dAT, ldda, ldda );            \
  } else {                                                      \
    magma_minproductblas_ztranspose2( dA, ldda, dAT, ldda, n, m );         \
    cublasFree(dAT);                                            \
  }

#endif /* _MAGMA_minproduct_TRANSPOSE_H_ */
