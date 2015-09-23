/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zpotri.cpp normal z -> d, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_minproduct.h"

#define PRECISION_d

/**
    Purpose
    -------
    DPOTRI computes the inverse of a real symmetric positive definite
    matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
    computed by DPOTRF.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  Upper triangle of A is stored;
      -     = Magma_minproductLower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the triangular factor U or L from the Cholesky
            factorization A = U**T*U or A = L*L**T, as computed by
            DPOTRF.
            On exit, the upper or lower triangle of the (symmetric)
            inverse of A, overwriting the input factor U or L.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, the (i,i) element of the factor U or L is
                  zero, and the inverse could not be computed.

    @ingroup magma_minproduct_dposv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_dpotri(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info)
{
    /* Local variables */
    *info = 0;
    if ((uplo != Magma_minproductUpper) && (uplo != Magma_minproductLower))
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,n))
        *info = -4;

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return *info;
    
    /* Invert the triangular Cholesky factor U or L */
    magma_minproduct_dtrtri( uplo, Magma_minproductNonUnit, n, A, lda, info );
    if ( *info == 0 ) {
        /* Form inv(U) * inv(U)**T or inv(L)**T * inv(L) */
        magma_minproduct_dlauum( uplo, n, A, lda, info );
    }
    
    return *info;
} /* magma_minproduct_dpotri */
