/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from ztrtri.cpp normal z -> c, Fri Jan 30 19:00:13 2015

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    CTRTRI computes the inverse of a real upper or lower triangular
    matrix A.

    This is the Level 3 BLAS version of the algorithm.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  A is upper triangular;
      -     = Magma_tally3Lower:  A is lower triangular.

    @param[in]
    diag    magma_tally3_diag_t
      -     = Magma_tally3NonUnit:  A is non-unit triangular;
      -     = Magma_tally3Unit:     A is unit triangular.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       COMPLEX array, dimension (LDA,N)
            On entry, the triangular matrix A.  If UPLO = Magma_tally3Upper, the
            leading N-by-N upper triangular part of the array A contains
            the upper triangular matrix, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally3Lower, the
            leading N-by-N lower triangular part of the array A contains
            the lower triangular matrix, and the strictly upper
            triangular part of A is not referenced.  If DIAG = Magma_tally3Unit, the
            diagonal elements of A are also not referenced and are
            assumed to be 1.
            On exit, the (triangular) inverse of the original matrix, in
            the same storage format.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -i, the i-th argument had an illegal value
      -     > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
                    matrix is singular and its inverse cannot be computed.

    @ingroup magma_tally3_cgesv_aux
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_ctrtri(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info)
{
    #define  A(i, j) ( A + (i) + (j)*lda )
    #define dA(i, j) (dA + (i) + (j)*ldda)

    /* Local variables */
    const char* uplo_ = lapack_uplo_const_tally3( uplo );
    const char* diag_ = lapack_diag_const_tally3( diag );
    magma_tally3_int_t     ldda, nb, nn, j, jb;
    magma_tally3FloatComplex c_zero     = MAGMA_tally3_C_ZERO;
    magma_tally3FloatComplex c_one      = MAGMA_tally3_C_ONE;
    magma_tally3FloatComplex c_neg_one  = MAGMA_tally3_C_NEG_ONE;
    magma_tally3FloatComplex *dA;

    int upper  = (uplo == Magma_tally3Upper);
    int nounit = (diag == Magma_tally3NonUnit);

    *info = 0;

    if (! upper && uplo != Magma_tally3Lower)
        *info = -1;
    else if (! nounit && diag != Magma_tally3Unit)
        *info = -2;
    else if (n < 0)
        *info = -3;
    else if (lda < max(1,n))
        *info = -5;

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ( n == 0 )
        return *info;

    /* Check for singularity if non-unit */
    if (nounit) {
        for (j=0; j < n; ++j) {
            if ( MAGMA_tally3_C_EQUAL( *A(j,j), c_zero )) {
                *info = j+1;  // Fortran index
                return *info;
            }
        }
    }

    /* Determine the block size for this environment */
    nb = magma_tally3_get_cpotrf_nb(n);

    ldda = ((n+31)/32)*32;
    if (MAGMA_tally3_SUCCESS != magma_tally3_cmalloc( &dA, (n)*ldda )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }

    magma_tally3_queue_t stream[2];
    magma_tally3_queue_create( &stream[0] );
    magma_tally3_queue_create( &stream[1] );

    if (nb <= 1 || nb >= n)
        lapackf77_ctrtri(uplo_, diag_, &n, A, &lda, info);
    else {
        if (upper) {
            /* Compute inverse of upper triangular matrix */
            for (j=0; j < n; j += nb) {
                jb = min(nb, (n-j));
                magma_tally3_csetmatrix( jb, (n-j),
                                  A(j, j),  lda,
                                  dA(j, j), ldda );

                /* Compute rows 1:j-1 of current block column */
                magma_tally3_ctrmm( Magma_tally3Left, Magma_tally3Upper,
                             Magma_tally3NoTrans, Magma_tally3NonUnit, j, jb,
                             c_one, dA(0,0), ldda, dA(0, j),ldda);

                magma_tally3_ctrsm( Magma_tally3Right, Magma_tally3Upper,
                             Magma_tally3NoTrans, Magma_tally3NonUnit, j, jb,
                             c_neg_one, dA(j,j), ldda, dA(0, j),ldda);

                magma_tally3_cgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        A(j, j),  lda, stream[1] );

                magma_tally3_cgetmatrix_async( j, jb,
                                        dA(0, j), ldda,
                                        A(0, j),  lda, stream[0] );

                magma_tally3_queue_sync( stream[1] );

                /* Compute inverse of current diagonal block */
                lapackf77_ctrtri(Magma_tally3UpperStr, diag_, &jb, A(j,j), &lda, info);

                magma_tally3_csetmatrix( jb, jb,
                                  A(j, j),  lda,
                                  dA(j, j), ldda );
            }
        }
        else {
            /* Compute inverse of lower triangular matrix */
            nn=((n-1)/nb)*nb+1;

            for (j=nn-1; j >= 0; j -= nb) {
                jb=min(nb,(n-j));

                if ((j+jb) < n) {
                    magma_tally3_csetmatrix( (n-j), jb,
                                      A(j, j),  lda,
                                      dA(j, j), ldda );

                    /* Compute rows j+jb:n of current block column */
                    magma_tally3_ctrmm( Magma_tally3Left, Magma_tally3Lower,
                                 Magma_tally3NoTrans, Magma_tally3NonUnit, (n-j-jb), jb,
                                 c_one, dA(j+jb,j+jb), ldda, dA(j+jb, j), ldda );

                    magma_tally3_ctrsm( Magma_tally3Right, Magma_tally3Lower,
                                 Magma_tally3NoTrans, Magma_tally3NonUnit, (n-j-jb), jb,
                                 c_neg_one, dA(j,j), ldda, dA(j+jb, j), ldda );

                    magma_tally3_cgetmatrix_async( n-j-jb, jb,
                                            dA(j+jb, j), ldda,
                                            A(j+jb, j),  lda, stream[1] );

                    magma_tally3_cgetmatrix_async( jb, jb,
                                            dA(j,j), ldda,
                                            A(j,j),  lda, stream[0] );

                    magma_tally3_queue_sync( stream[0] );
                }

                /* Compute inverse of current diagonal block */
                lapackf77_ctrtri(Magma_tally3LowerStr, diag_, &jb, A(j,j), &lda, info);

                magma_tally3_csetmatrix( jb, jb,
                                  A(j, j),  lda,
                                  dA(j, j), ldda );
            }
        }
    }

    magma_tally3_queue_destroy( stream[0] );
    magma_tally3_queue_destroy( stream[1] );
    magma_tally3_free( dA );

    return *info;
}