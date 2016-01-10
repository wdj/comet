/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZLAUUM computes the product U * U' or L' * L, where the triangular
    factor U or L is stored in the upper or lower triangular part of
    the array A.

    If UPLO = Magma_tally4Upper then the upper triangle of the result is stored,
    overwriting the factor U in A.
    If UPLO = Magma_tally4Lower then the lower triangle of the result is stored,
    overwriting the factor L in A.
    This is the blocked form of the algorithm, calling Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
            Specifies whether the triangular factor stored in the array A
            is upper or lower triangular:
      -     = Magma_tally4Upper:  Upper triangular
      -     = Magma_tally4Lower:  Lower triangular

    @param[in]
    n       INTEGER
            The order of the triangular factor U or L.  N >= 0.

    @param[in,out]
    A       COPLEX_16 array, dimension (LDA,N)
            On entry, the triangular factor U or L.
            On exit, if UPLO = Magma_tally4Upper, the upper triangle of A is
            overwritten with the upper triangle of the product U * U';
            if UPLO = Magma_tally4Lower, the lower triangle of A is overwritten with
            the lower triangle of the product L' * L.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -k, the k-th argument had an illegal value

    @ingroup magma_tally4_zposv_aux
    ***************************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zlauum(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info)
{
#define A(i, j)  (A  + (j)*lda  + (i))
#define dA(i, j) (dA + (j)*ldda + (i))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const_tally4( uplo );
    magma_tally4_int_t     ldda, nb;
    magma_tally4_int_t i, ib;
    magma_tally4DoubleComplex c_one = MAGMA_tally4_Z_ONE;
    double             d_one = MAGMA_tally4_D_ONE;
    magma_tally4DoubleComplex    *dA;
    int upper = (uplo == Magma_tally4Upper);

    *info = 0;
    if (! upper && uplo != Magma_tally4Lower)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,n))
        *info = -4;

    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ( n == 0 )
        return *info;

    ldda = ((n+31)/32)*32;

    if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc( &dA, (n)*ldda )) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }

    magma_tally4_queue_t stream[2];
    magma_tally4_queue_create( &stream[0] );
    magma_tally4_queue_create( &stream[1] );

    nb = magma_tally4_get_zpotrf_nb(n);

    if (nb <= 1 || nb >= n)
        lapackf77_zlauum(uplo_, &n, A, &lda, info);
    else {
        if (upper) {
            /* Compute the product U * U'. */
            for (i=0; i < n; i += nb) {
                ib=min(nb,n-i);

                magma_tally4_zsetmatrix_async( ib, ib,
                                        A(i,i),   lda,
                                        dA(i, i), ldda, stream[1] );

                magma_tally4_zsetmatrix_async( ib, (n-i-ib),
                                        A(i,i+ib),  lda,
                                        dA(i,i+ib), ldda, stream[0] );

                magma_tally4_queue_sync( stream[1] );

                magma_tally4_ztrmm( Magma_tally4Right, Magma_tally4Upper,
                             Magma_tally4ConjTrans, Magma_tally4NonUnit, i, ib,
                             c_one, dA(i,i), ldda, dA(0, i),ldda);


                lapackf77_zlauum(Magma_tally4UpperStr, &ib, A(i,i), &lda, info);

                magma_tally4_zsetmatrix_async( ib, ib,
                                        A(i, i),  lda,
                                        dA(i, i), ldda, stream[0] );

                if (i+ib < n) {
                    magma_tally4_zgemm( Magma_tally4NoTrans, Magma_tally4ConjTrans,
                                 i, ib, (n-i-ib), c_one, dA(0,i+ib),
                                 ldda, dA(i, i+ib),ldda, c_one,
                                 dA(0,i), ldda);

                    magma_tally4_queue_sync( stream[0] );

                    magma_tally4_zherk( Magma_tally4Upper, Magma_tally4NoTrans, ib,(n-i-ib),
                                 d_one, dA(i, i+ib), ldda,
                                 d_one, dA(i, i), ldda);
                }

                magma_tally4_zgetmatrix( i+ib, ib,
                                  dA(0, i), ldda,
                                  A(0, i),  lda );
            }
        }
        else {
            /* Compute the product L' * L. */
            for (i=0; i < n; i += nb) {
                ib=min(nb,n-i);
                magma_tally4_zsetmatrix_async( ib, ib,
                                        A(i,i),   lda,
                                        dA(i, i), ldda, stream[1] );

                magma_tally4_zsetmatrix_async( (n-i-ib), ib,
                                        A(i+ib, i),  lda,
                                        dA(i+ib, i), ldda, stream[0] );

                magma_tally4_queue_sync( stream[1] );

                magma_tally4_ztrmm( Magma_tally4Left, Magma_tally4Lower,
                             Magma_tally4ConjTrans, Magma_tally4NonUnit, ib,
                             i, c_one, dA(i,i), ldda,
                             dA(i, 0),ldda);


                lapackf77_zlauum(Magma_tally4LowerStr, &ib, A(i,i), &lda, info);

                magma_tally4_zsetmatrix_async( ib, ib,
                                        A(i, i),  lda,
                                        dA(i, i), ldda, stream[0] );

                if (i+ib < n) {
                    magma_tally4_zgemm(Magma_tally4ConjTrans, Magma_tally4NoTrans,
                                    ib, i, (n-i-ib), c_one, dA( i+ib,i),
                                    ldda, dA(i+ib, 0),ldda, c_one,
                                    dA(i,0), ldda);

                    magma_tally4_queue_sync( stream[0] );

                    magma_tally4_zherk(Magma_tally4Lower, Magma_tally4ConjTrans, ib, (n-i-ib),
                                    d_one, dA(i+ib, i), ldda,
                                    d_one, dA(i, i), ldda);
                }
                magma_tally4_zgetmatrix( ib, i+ib,
                                  dA(i, 0), ldda,
                                  A(i, 0),  lda );
            }
        }
    }
    magma_tally4_queue_destroy( stream[0] );
    magma_tally4_queue_destroy( stream[1] );

    magma_tally4_free( dA );

    return *info;
}