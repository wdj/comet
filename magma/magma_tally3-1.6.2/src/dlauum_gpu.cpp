/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlauum_gpu.cpp normal z -> d, Fri Jan 30 19:00:13 2015

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    DLAUUM computes the product U * U' or L' * L, where the triangular
    factor U or L is stored in the upper or lower triangular part of
    the array dA.

    If UPLO = Magma_tally3Upper then the upper triangle of the result is stored,
    overwriting the factor U in dA.
    If UPLO = Magma_tally3Lower then the lower triangle of the result is stored,
    overwriting the factor L in dA.
    This is the blocked form of the algorithm, calling Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally3_uplo_t
            Specifies whether the triangular factor stored in the array dA
            is upper or lower triangular:
      -     = Magma_tally3Upper:  Upper triangular
      -     = Magma_tally3Lower:  Lower triangular

    @param[in]
    n       INTEGER
            The order of the triangular factor U or L.  N >= 0.

    @param[in,out]
    dA      DOUBLE PRECISION array on the GPU, dimension (LDDA,N)
            On entry, the triangular factor U or L.
            On exit, if UPLO = Magma_tally3Upper, the upper triangle of dA is
            overwritten with the upper triangle of the product U * U';
            if UPLO = Magma_tally3Lower, the lower triangle of dA is overwritten with
            the lower triangle of the product L' * L.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -k, the k-th argument had an illegal value

    @ingroup magma_tally3_dposv_aux
    ***************************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_dlauum_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info)
{
#define dA(i, j) (dA + (j)*ldda + (i))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const_tally3( uplo );
    magma_tally3_int_t         nb, i, ib;
    double              d_one = MAGMA_tally3_D_ONE;
    double  c_one = MAGMA_tally3_D_ONE;
    double  *work;

    int upper  = (uplo == Magma_tally3Upper);

    *info = 0;

    if (! upper && uplo != Magma_tally3Lower)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,n))
        *info = -4;

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    nb = magma_tally3_get_dpotrf_nb(n);

    if (MAGMA_tally3_SUCCESS != magma_tally3_dmalloc_pinned( &work, nb*nb )) {
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        return *info;
    }

    magma_tally3_queue_t stream[2];
    magma_tally3_queue_create( &stream[0] );
    magma_tally3_queue_create( &stream[1] );

    if (nb <= 1 || nb >= n) {
        magma_tally3_dgetmatrix( n, n, dA, ldda, work, n );
        lapackf77_dlauum(uplo_, &n, work, &n, info);
        magma_tally3_dsetmatrix( n, n, work, n, dA, ldda );
    }
    else {
        if (upper) {
            /* Compute inverse of upper triangular matrix */
            for (i=0; i < n; i += nb) {
                ib = min(nb, (n-i));

                /* Compute the product U * U'. */
                magma_tally3_dtrmm( Magma_tally3Right, Magma_tally3Upper,
                             Magma_tally3ConjTrans, Magma_tally3NonUnit, i, ib,
                             c_one, dA(i,i), ldda, dA(0, i),ldda);

                magma_tally3_dgetmatrix( ib, ib,
                                  dA(i, i), ldda,
                                  work,     ib );

                lapackf77_dlauum(Magma_tally3UpperStr, &ib, work, &ib, info);

                magma_tally3_dsetmatrix( ib, ib,
                                  work,     ib,
                                  dA(i, i), ldda );

                if (i+ib < n) {
                    magma_tally3_dgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans,
                                 i, ib, (n-i-ib), c_one, dA(0,i+ib),
                                 ldda, dA(i, i+ib), ldda, c_one,
                                 dA(0,i), ldda);

                    magma_tally3_dsyrk( Magma_tally3Upper, Magma_tally3NoTrans, ib,(n-i-ib),
                                 d_one, dA(i, i+ib), ldda,
                                 d_one, dA(i, i),    ldda);
                }
            }
        }
        else {
            /* Compute the product L' * L. */
            for (i=0; i < n; i += nb) {
                ib=min(nb,(n-i));

                magma_tally3_dtrmm( Magma_tally3Left, Magma_tally3Lower,
                             Magma_tally3ConjTrans, Magma_tally3NonUnit, ib,
                             i, c_one, dA(i,i), ldda,
                             dA(i, 0),ldda);

                magma_tally3_dgetmatrix( ib, ib,
                                  dA(i, i), ldda,
                                  work,     ib );

                lapackf77_dlauum(Magma_tally3LowerStr, &ib, work, &ib, info);

                magma_tally3_dsetmatrix( ib, ib,
                                  work,     ib,
                                  dA(i, i), ldda );

                if (i+ib < n) {
                    magma_tally3_dgemm( Magma_tally3ConjTrans, Magma_tally3NoTrans,
                                 ib, i, (n-i-ib), c_one, dA( i+ib,i),
                                 ldda, dA(i+ib, 0),ldda, c_one,
                                 dA(i,0), ldda);
                    magma_tally3_dsyrk( Magma_tally3Lower, Magma_tally3ConjTrans, ib, (n-i-ib),
                                 d_one, dA(i+ib, i), ldda,
                                 d_one, dA(i, i),    ldda);
                }
            }
        }
    }

    magma_tally3_queue_destroy( stream[0] );
    magma_tally3_queue_destroy( stream[1] );

    magma_tally3_free_pinned( work );

    return *info;
}
