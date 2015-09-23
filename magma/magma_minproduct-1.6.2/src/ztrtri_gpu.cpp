/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    ZTRTRI computes the inverse of a real upper or lower triangular
    matrix dA.

    This is the Level 3 BLAS version of the algorithm.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  A is upper triangular;
      -     = Magma_minproductLower:  A is lower triangular.

    @param[in]
    diag    magma_minproduct_diag_t
      -     = Magma_minproductNonUnit:  A is non-unit triangular;
      -     = Magma_minproductUnit:     A is unit triangular.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX_16 array ON THE GPU, dimension (LDDA,N)
            On entry, the triangular matrix A.  If UPLO = Magma_minproductUpper, the
            leading N-by-N upper triangular part of the array dA contains
            the upper triangular matrix, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_minproductLower, the
            leading N-by-N lower triangular part of the array dA contains
            the lower triangular matrix, and the strictly upper
            triangular part of A is not referenced.  If DIAG = Magma_minproductUnit, the
            diagonal elements of A are also not referenced and are
            assumed to be 1.
            On exit, the (triangular) inverse of the original matrix, in
            the same storage format.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -i, the i-th argument had an illegal value
      -     > 0: if INFO = i, dA(i,i) is exactly zero.  The triangular
                    matrix is singular and its inverse cannot be computed.
                 (Singularity check is currently disabled.)

    @ingroup magma_minproduct_zgesv_aux
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_ztrtri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info)
{
#define dA(i, j) (dA+(j)*ldda + (i))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const( uplo );
    const char* diag_ = lapack_diag_const( diag );
    magma_minproduct_int_t nb, nn, j, jb;
    //magma_minproductDoubleComplex c_zero     = MAGMA_minproduct_Z_ZERO;
    magma_minproductDoubleComplex c_one      = MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex c_neg_one  = MAGMA_minproduct_Z_NEG_ONE;
    magma_minproductDoubleComplex *work;

    int upper  = (uplo == Magma_minproductUpper);
    int nounit = (diag == Magma_minproductNonUnit);

    *info = 0;

    if (! upper && uplo != Magma_minproductLower)
        *info = -1;
    else if (! nounit && diag != Magma_minproductUnit)
        *info = -2;
    else if (n < 0)
        *info = -3;
    else if (ldda < max(1,n))
        *info = -5;

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Check for singularity if non-unit */
    /* cannot do here with matrix dA on GPU -- need kernel */
    /*
    if (nounit) {
        for (j=0; j < n; ++j) {
            if ( MAGMA_minproduct_Z_EQUAL( *dA(j,j), c_zero )) {
                *info = j+1;  // Fortran index
                return *info;
            }
        }
    }
    */

    /* Determine the block size for this environment */
    nb = magma_minproduct_get_zpotrf_nb(n);

    if (MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc_pinned( &work, nb*nb )) {
        *info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return *info;
    }

    magma_minproduct_queue_t stream[2];
    magma_minproduct_queue_create( &stream[0] );
    magma_minproduct_queue_create( &stream[1] );

    if (nb <= 1 || nb >= n) {
        magma_minproduct_zgetmatrix( n, n, dA, ldda, work, n );
        lapackf77_ztrtri( uplo_, diag_, &n, work, &n, info );
        magma_minproduct_zsetmatrix( n, n, work, n, dA, ldda );
    }
    else {
        if (upper) {
            /* Compute inverse of upper triangular matrix */
            for (j=0; j < n; j += nb) {
                jb = min(nb, (n-j));

                /* Compute rows 1:j-1 of current block column */
                magma_minproduct_ztrmm( Magma_minproductLeft, Magma_minproductUpper,
                             Magma_minproductNoTrans, Magma_minproductNonUnit, j, jb,
                             c_one, dA(0,0), ldda, dA(0, j), ldda );

                magma_minproduct_ztrsm( Magma_minproductRight, Magma_minproductUpper,
                             Magma_minproductNoTrans, Magma_minproductNonUnit, j, jb,
                             c_neg_one, dA(j,j), ldda, dA(0, j), ldda );

                magma_minproduct_zgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[1] );

                magma_minproduct_queue_sync( stream[1] );

                /* Compute inverse of current diagonal block */
                lapackf77_ztrtri( Magma_minproductUpperStr, diag_, &jb, work, &jb, info );

                magma_minproduct_zsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[0] );
            }
        }
        else {
            /* Compute inverse of lower triangular matrix */
            nn = ((n-1)/nb)*nb+1;

            for (j=nn-1; j >= 0; j -= nb) {
                jb = min(nb,(n-j));

                if ((j+jb) < n) {
                    /* Compute rows j+jb:n of current block column */
                    magma_minproduct_ztrmm( Magma_minproductLeft, Magma_minproductLower,
                                 Magma_minproductNoTrans, Magma_minproductNonUnit, (n-j-jb), jb,
                                 c_one, dA(j+jb,j+jb), ldda, dA(j+jb, j), ldda );

                    magma_minproduct_ztrsm( Magma_minproductRight, Magma_minproductLower,
                                 Magma_minproductNoTrans, Magma_minproductNonUnit, (n-j-jb), jb,
                                 c_neg_one, dA(j,j), ldda, dA(j+jb, j), ldda );
                }

                magma_minproduct_zgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[1] );

                magma_minproduct_queue_sync( stream[1] );

                /* Compute inverse of current diagonal block */
                lapackf77_ztrtri( Magma_minproductLowerStr, diag_, &jb, work, &jb, info );

                magma_minproduct_zsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[0] );
            }
        }
    }

    magma_minproduct_queue_destroy( stream[0] );
    magma_minproduct_queue_destroy( stream[1] );
    magma_minproduct_free_pinned( work );

    return *info;
}
