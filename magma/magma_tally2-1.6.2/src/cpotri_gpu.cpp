/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zpotri_gpu.cpp normal z -> c, Fri Jan 30 19:00:12 2015

*/
#include "common_magma_tally2.h"

#define PRECISION_c

/**
    Purpose
    -------
    CPOTRI computes the inverse of a real symmetric positive definite
    matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
    computed by CPOTRF.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangle of A is stored;
      -     = Magma_tally2Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX array on the GPU, dimension (LDA,N)
            On entry, the triangular factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T, as computed by
            CPOTRF.
            On exit, the upper or lower triangle of the (symmetric)
            inverse of A, overwriting the input factor U or L.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, the (i,i) element of the factor U or L is
                  zero, and the inverse could not be computed.

    @ingroup magma_tally2_cposv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_cpotri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info)
{
    /* Local variables */
    *info = 0;
    if ((uplo != Magma_tally2Upper) && (uplo != Magma_tally2Lower))
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,n))
        *info = -4;

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return *info;
    
    /* Invert the triangular Cholesky factor U or L */
    magma_tally2_ctrtri_gpu( uplo, Magma_tally2NonUnit, n, dA, ldda, info );
    if ( *info == 0 ) {
        /* Form inv(U) * inv(U)**T or inv(L)**T * inv(L) */
        magma_tally2_clauum_gpu( uplo, n, dA, ldda, info );
    }
    
    return *info;
} /* magma_tally2_cpotri */