/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from ztrtri_gpu.cpp normal z -> d, Fri Jan 30 19:00:13 2015

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    DTRTRI computes the inverse of a real upper or lower triangular
    matrix dA.

    This is the Level 3 BLAS version of the algorithm.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  A is upper triangular;
      -     = Magma_tally4Lower:  A is lower triangular.

    @param[in]
    diag    magma_tally4_diag_t
      -     = Magma_tally4NonUnit:  A is non-unit triangular;
      -     = Magma_tally4Unit:     A is unit triangular.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array ON THE GPU, dimension (LDDA,N)
            On entry, the triangular matrix A.  If UPLO = Magma_tally4Upper, the
            leading N-by-N upper triangular part of the array dA contains
            the upper triangular matrix, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally4Lower, the
            leading N-by-N lower triangular part of the array dA contains
            the lower triangular matrix, and the strictly upper
            triangular part of A is not referenced.  If DIAG = Magma_tally4Unit, the
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

    @ingroup magma_tally4_dgesv_aux
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dtrtri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info)
{
#define dA(i, j) (dA+(j)*ldda + (i))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const_tally4( uplo );
    const char* diag_ = lapack_diag_const_tally4( diag );
    magma_tally4_int_t nb, nn, j, jb;
    //double c_zero     = MAGMA_tally4_D_ZERO;
    double c_one      = MAGMA_tally4_D_ONE;
    double c_neg_one  = MAGMA_tally4_D_NEG_ONE;
    double *work;

    int upper  = (uplo == Magma_tally4Upper);
    int nounit = (diag == Magma_tally4NonUnit);

    *info = 0;

    if (! upper && uplo != Magma_tally4Lower)
        *info = -1;
    else if (! nounit && diag != Magma_tally4Unit)
        *info = -2;
    else if (n < 0)
        *info = -3;
    else if (ldda < max(1,n))
        *info = -5;

    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Check for singularity if non-unit */
    /* cannot do here with matrix dA on GPU -- need kernel */
    /*
    if (nounit) {
        for (j=0; j < n; ++j) {
            if ( MAGMA_tally4_D_EQUAL( *dA(j,j), c_zero )) {
                *info = j+1;  // Fortran index
                return *info;
            }
        }
    }
    */

    /* Determine the block size for this environment */
    nb = magma_tally4_get_dpotrf_nb(n);

    if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc_pinned( &work, nb*nb )) {
        *info = MAGMA_tally4_ERR_HOST_ALLOC;
        return *info;
    }

    magma_tally4_queue_t stream[2];
    magma_tally4_queue_create( &stream[0] );
    magma_tally4_queue_create( &stream[1] );

    if (nb <= 1 || nb >= n) {
        magma_tally4_dgetmatrix( n, n, dA, ldda, work, n );
        lapackf77_dtrtri( uplo_, diag_, &n, work, &n, info );
        magma_tally4_dsetmatrix( n, n, work, n, dA, ldda );
    }
    else {
        if (upper) {
            /* Compute inverse of upper triangular matrix */
            for (j=0; j < n; j += nb) {
                jb = min(nb, (n-j));

                /* Compute rows 1:j-1 of current block column */
                magma_tally4_dtrmm( Magma_tally4Left, Magma_tally4Upper,
                             Magma_tally4NoTrans, Magma_tally4NonUnit, j, jb,
                             c_one, dA(0,0), ldda, dA(0, j), ldda );

                magma_tally4_dtrsm( Magma_tally4Right, Magma_tally4Upper,
                             Magma_tally4NoTrans, Magma_tally4NonUnit, j, jb,
                             c_neg_one, dA(j,j), ldda, dA(0, j), ldda );

                magma_tally4_dgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[1] );

                magma_tally4_queue_sync( stream[1] );

                /* Compute inverse of current diagonal block */
                lapackf77_dtrtri( Magma_tally4UpperStr, diag_, &jb, work, &jb, info );

                magma_tally4_dsetmatrix_async( jb, jb,
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
                    magma_tally4_dtrmm( Magma_tally4Left, Magma_tally4Lower,
                                 Magma_tally4NoTrans, Magma_tally4NonUnit, (n-j-jb), jb,
                                 c_one, dA(j+jb,j+jb), ldda, dA(j+jb, j), ldda );

                    magma_tally4_dtrsm( Magma_tally4Right, Magma_tally4Lower,
                                 Magma_tally4NoTrans, Magma_tally4NonUnit, (n-j-jb), jb,
                                 c_neg_one, dA(j,j), ldda, dA(j+jb, j), ldda );
                }

                magma_tally4_dgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[1] );

                magma_tally4_queue_sync( stream[1] );

                /* Compute inverse of current diagonal block */
                lapackf77_dtrtri( Magma_tally4LowerStr, diag_, &jb, work, &jb, info );

                magma_tally4_dsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[0] );
            }
        }
    }

    magma_tally4_queue_destroy( stream[0] );
    magma_tally4_queue_destroy( stream[1] );
    magma_tally4_free_pinned( work );

    return *info;
}
