/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mathieu Faverge
 
       Macro to transpose matrices before and after computation
       in LU kernels
*/

#ifndef MAGMA_tally4_TRANSPOSE_H
#define MAGMA_tally4_TRANSPOSE_H

// TODO shouldn't need % 32 == 0 checks anymore; inplace works for any m == n.

#define magma_tally4blas_sgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally4blas_stranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    if (MAGMA_tally4_SUCCESS != magma_tally4_smalloc(&dAT, (m)*(n)) )         \
      return MAGMA_tally4_ERR_DEVICE_ALLOC;                            \
    magma_tally4blas_stranspose( m, n, dA, ldda, dAT, ldda );          \
  }

#define magma_tally4blas_sgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally4blas_stranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    magma_tally4blas_stranspose( n, m, dAT, ldda, dA, ldda );          \
    magma_tally4_free(dAT);                                            \
  }

#define magma_tally4blas_dgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally4blas_dtranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc(&dAT, (m)*(n)))          \
      return MAGMA_tally4_ERR_DEVICE_ALLOC;                            \
    magma_tally4blas_dtranspose( m, n, dA, ldda, dAT, ldda );          \
  }

#define magma_tally4blas_dgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally4blas_dtranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    magma_tally4blas_dtranspose( n, m, dAT, ldda, dA, ldda );          \
    magma_tally4_free(dAT);                                            \
  }

#define magma_tally4blas_cgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally4blas_ctranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    if (MAGMA_tally4_SUCCESS != magma_tally4_cmalloc(&dAT, (m)*(n)))          \
      return MAGMA_tally4_ERR_DEVICE_ALLOC;                            \
    magma_tally4blas_ctranspose( m, n, dA, ldda, dAT, ldda );          \
  }

#define magma_tally4blas_cgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally4blas_ctranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    magma_tally4blas_ctranspose( n, m, dAT, ldda, dA, ldda );          \
    magma_tally4_free(dAT);                                            \
  }

#define magma_tally4blas_zgetmo_in( dA, dAT, ldda, m, n )              \
  dAT = dA;                                                     \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally4blas_ztranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc(&dAT, (m)*(n)))          \
      return MAGMA_tally4_ERR_DEVICE_ALLOC;                            \
    magma_tally4blas_ztranspose( m, n, dA, ldda, dAT, ldda );          \
  }

#define magma_tally4blas_zgetmo_out( dA, dAT, ldda, m, n )             \
  if ( ( (m) == (n) ) && ( (m)%32 == 0) && ( (ldda)%32 == 0) ){ \
    magma_tally4blas_ztranspose_inplace( ldda, dAT, ldda );            \
  } else {                                                      \
    magma_tally4blas_ztranspose( n, m, dAT, ldda, dA, ldda );          \
    magma_tally4_free(dAT);                                            \
  }

#endif /* MAGMA_tally4_TRANSPOSE_H */
