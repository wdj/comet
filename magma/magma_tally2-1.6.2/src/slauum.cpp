/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlauum.cpp normal z -> s, Fri Jan 30 19:00:13 2015

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    SLAUUM computes the product U * U' or L' * L, where the triangular
    factor U or L is stored in the upper or lower triangular part of
    the array A.

    If UPLO = Magma_tally2Upper then the upper triangle of the result is stored,
    overwriting the factor U in A.
    If UPLO = Magma_tally2Lower then the lower triangle of the result is stored,
    overwriting the factor L in A.
    This is the blocked form of the algorithm, calling Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally2_uplo_t
            Specifies whether the triangular factor stored in the array A
            is upper or lower triangular:
      -     = Magma_tally2Upper:  Upper triangular
      -     = Magma_tally2Lower:  Lower triangular

    @param[in]
    n       INTEGER
            The order of the triangular factor U or L.  N >= 0.

    @param[in,out]
    A       COPLEX_16 array, dimension (LDA,N)
            On entry, the triangular factor U or L.
            On exit, if UPLO = Magma_tally2Upper, the upper triangle of A is
            overwritten with the upper triangle of the product U * U';
            if UPLO = Magma_tally2Lower, the lower triangle of A is overwritten with
            the lower triangle of the product L' * L.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -k, the k-th argument had an illegal value

    @ingroup magma_tally2_sposv_aux
    ***************************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_slauum(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info)
{
#define A(i, j)  (A  + (j)*lda  + (i))
#define dA(i, j) (dA + (j)*ldda + (i))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const_tally2( uplo );
    magma_tally2_int_t     ldda, nb;
    magma_tally2_int_t i, ib;
    float c_one = MAGMA_tally2_S_ONE;
    float             d_one = MAGMA_tally2_D_ONE;
    float    *dA;
    int upper = (uplo == Magma_tally2Upper);

    *info = 0;
    if (! upper && uplo != Magma_tally2Lower)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (lda < max(1,n))
        *info = -4;

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ( n == 0 )
        return *info;

    ldda = ((n+31)/32)*32;

    if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc( &dA, (n)*ldda )) {
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }

    magma_tally2_queue_t stream[2];
    magma_tally2_queue_create( &stream[0] );
    magma_tally2_queue_create( &stream[1] );

    nb = magma_tally2_get_spotrf_nb(n);

    if (nb <= 1 || nb >= n)
        lapackf77_slauum(uplo_, &n, A, &lda, info);
    else {
        if (upper) {
            /* Compute the product U * U'. */
            for (i=0; i < n; i += nb) {
                ib=min(nb,n-i);

                magma_tally2_ssetmatrix_async( ib, ib,
                                        A(i,i),   lda,
                                        dA(i, i), ldda, stream[1] );

                magma_tally2_ssetmatrix_async( ib, (n-i-ib),
                                        A(i,i+ib),  lda,
                                        dA(i,i+ib), ldda, stream[0] );

                magma_tally2_queue_sync( stream[1] );

                magma_tally2_strmm( Magma_tally2Right, Magma_tally2Upper,
                             Magma_tally2ConjTrans, Magma_tally2NonUnit, i, ib,
                             c_one, dA(i,i), ldda, dA(0, i),ldda);


                lapackf77_slauum(Magma_tally2UpperStr, &ib, A(i,i), &lda, info);

                magma_tally2_ssetmatrix_async( ib, ib,
                                        A(i, i),  lda,
                                        dA(i, i), ldda, stream[0] );

                if (i+ib < n) {
                    magma_tally2_sgemm( Magma_tally2NoTrans, Magma_tally2ConjTrans,
                                 i, ib, (n-i-ib), c_one, dA(0,i+ib),
                                 ldda, dA(i, i+ib),ldda, c_one,
                                 dA(0,i), ldda);

                    magma_tally2_queue_sync( stream[0] );

                    magma_tally2_ssyrk( Magma_tally2Upper, Magma_tally2NoTrans, ib,(n-i-ib),
                                 d_one, dA(i, i+ib), ldda,
                                 d_one, dA(i, i), ldda);
                }

                magma_tally2_sgetmatrix( i+ib, ib,
                                  dA(0, i), ldda,
                                  A(0, i),  lda );
            }
        }
        else {
            /* Compute the product L' * L. */
            for (i=0; i < n; i += nb) {
                ib=min(nb,n-i);
                magma_tally2_ssetmatrix_async( ib, ib,
                                        A(i,i),   lda,
                                        dA(i, i), ldda, stream[1] );

                magma_tally2_ssetmatrix_async( (n-i-ib), ib,
                                        A(i+ib, i),  lda,
                                        dA(i+ib, i), ldda, stream[0] );

                magma_tally2_queue_sync( stream[1] );

                magma_tally2_strmm( Magma_tally2Left, Magma_tally2Lower,
                             Magma_tally2ConjTrans, Magma_tally2NonUnit, ib,
                             i, c_one, dA(i,i), ldda,
                             dA(i, 0),ldda);


                lapackf77_slauum(Magma_tally2LowerStr, &ib, A(i,i), &lda, info);

                magma_tally2_ssetmatrix_async( ib, ib,
                                        A(i, i),  lda,
                                        dA(i, i), ldda, stream[0] );

                if (i+ib < n) {
                    magma_tally2_sgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans,
                                    ib, i, (n-i-ib), c_one, dA( i+ib,i),
                                    ldda, dA(i+ib, 0),ldda, c_one,
                                    dA(i,0), ldda);

                    magma_tally2_queue_sync( stream[0] );

                    magma_tally2_ssyrk(Magma_tally2Lower, Magma_tally2ConjTrans, ib, (n-i-ib),
                                    d_one, dA(i+ib, i), ldda,
                                    d_one, dA(i, i), ldda);
                }
                magma_tally2_sgetmatrix( ib, i+ib,
                                  dA(i, 0), ldda,
                                  A(i, 0),  lda );
            }
        }
    }
    magma_tally2_queue_destroy( stream[0] );
    magma_tally2_queue_destroy( stream[1] );

    magma_tally2_free( dA );

    return *info;
}
