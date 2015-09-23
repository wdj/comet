/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mathieu Faverge
 
       Macro to transpose matrices before and after computation
       in LU kernels
*/

#ifndef MAGMA_minproduct_TRANSPOSE_H
#define MAGMA_minproduct_TRANSPOSE_H

// TODO shouldn't need % 32 == 0 checks anymore; inplace works for any m == n.

#define magma_minproductblas_sgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_stranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_smalloc(&dAT, (m)*(n)) )         \
      return MAGMA_minproduct_ERR_DEVICE_ALLOC;                            \
    magma_minproductblas_stranspose( m, n, dA, ldda, dAT, ldda );          \
  }

#define magma_minproductblas_sgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_stranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    magma_minproductblas_stranspose( n, m, dAT, ldda, dA, ldda );          \
    magma_minproduct_free(dAT);                                            \
  }

#define magma_minproductblas_dgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_dtranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_dmalloc(&dAT, (m)*(n)))          \
      return MAGMA_minproduct_ERR_DEVICE_ALLOC;                            \
    magma_minproductblas_dtranspose( m, n, dA, ldda, dAT, ldda );          \
  }

#define magma_minproductblas_dgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_dtranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    magma_minproductblas_dtranspose( n, m, dAT, ldda, dA, ldda );          \
    magma_minproduct_free(dAT);                                            \
  }

#define magma_minproductblas_cgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_ctranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc(&dAT, (m)*(n)))          \
      return MAGMA_minproduct_ERR_DEVICE_ALLOC;                            \
    magma_minproductblas_ctranspose( m, n, dA, ldda, dAT, ldda );          \
  }

#define magma_minproductblas_cgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_ctranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    magma_minproductblas_ctranspose( n, m, dAT, ldda, dA, ldda );          \
    magma_minproduct_free(dAT);                                            \
  }

#define magma_minproductblas_zgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_ztranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc(&dAT, (m)*(n)))          \
      return MAGMA_minproduct_ERR_DEVICE_ALLOC;                            \
    magma_minproductblas_ztranspose( m, n, dA, ldda, dAT, ldda );          \
  }

#define magma_minproductblas_zgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_minproductblas_ztranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    magma_minproductblas_ztranspose( n, m, dAT, ldda, dA, ldda );          \
    magma_minproduct_free(dAT);                                            \
  }

#endif /* MAGMA_minproduct_TRANSPOSE_H */
