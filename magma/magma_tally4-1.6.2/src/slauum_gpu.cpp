/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlauum_gpu.cpp normal z -> s, Fri Jan 30 19:00:12 2015

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    SLAUUM computes the product U * U' or L' * L, where the triangular
    factor U or L is stored in the upper or lower triangular part of
    the array dA.

    If UPLO = Magma_tally4Upper then the upper triangle of the result is stored,
    overwriting the factor U in dA.
    If UPLO = Magma_tally4Lower then the lower triangle of the result is stored,
    overwriting the factor L in dA.
    This is the blocked form of the algorithm, calling Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
            Specifies whether the triangular factor stored in the array dA
            is upper or lower triangular:
      -     = Magma_tally4Upper:  Upper triangular
      -     = Magma_tally4Lower:  Lower triangular

    @param[in]
    n       INTEGER
            The order of the triangular factor U or L.  N >= 0.

    @param[in,out]
    dA      REAL array on the GPU, dimension (LDDA,N)
            On entry, the triangular factor U or L.
            On exit, if UPLO = Magma_tally4Upper, the upper triangle of dA is
            overwritten with the upper triangle of the product U * U';
            if UPLO = Magma_tally4Lower, the lower triangle of dA is overwritten with
            the lower triangle of the product L' * L.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -k, the k-th argument had an illegal value

    @ingroup magma_tally4_sposv_aux
    ***************************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_slauum_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info)
{
#define dA(i, j) (dA + (j)*ldda + (i))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const_tally4( uplo );
    magma_tally4_int_t         nb, i, ib;
    float              d_one = MAGMA_tally4_D_ONE;
    float  c_one = MAGMA_tally4_S_ONE;
    float  *work;

    int upper  = (uplo == Magma_tally4Upper);

    *info = 0;

    if (! upper && uplo != Magma_tally4Lower)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,n))
        *info = -4;

    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    nb = magma_tally4_get_spotrf_nb(n);

    if (MAGMA_tally4_SUCCESS != magma_tally4_smalloc_pinned( &work, nb*nb )) {
        *info = MAGMA_tally4_ERR_HOST_ALLOC;
        return *info;
    }

    magma_tally4_queue_t stream[2];
    magma_tally4_queue_create( &stream[0] );
    magma_tally4_queue_create( &stream[1] );

    if (nb <= 1 || nb >= n) {
        magma_tally4_sgetmatrix( n, n, dA, ldda, work, n );
        lapackf77_slauum(uplo_, &n, work, &n, info);
        magma_tally4_ssetmatrix( n, n, work, n, dA, ldda );
    }
    else {
        if (upper) {
            /* Compute inverse of upper triangular matrix */
            for (i=0; i < n; i += nb) {
                ib = min(nb, (n-i));

                /* Compute the product U * U'. */
                magma_tally4_strmm( Magma_tally4Right, Magma_tally4Upper,
                             Magma_tally4ConjTrans, Magma_tally4NonUnit, i, ib,
                             c_one, dA(i,i), ldda, dA(0, i),ldda);

                magma_tally4_sgetmatrix( ib, ib,
                                  dA(i, i), ldda,
                                  work,     ib );

                lapackf77_slauum(Magma_tally4UpperStr, &ib, work, &ib, info);

                magma_tally4_ssetmatrix( ib, ib,
                                  work,     ib,
                                  dA(i, i), ldda );

                if (i+ib < n) {
                    magma_tally4_sgemm( Magma_tally4NoTrans, Magma_tally4ConjTrans,
                                 i, ib, (n-i-ib), c_one, dA(0,i+ib),
                                 ldda, dA(i, i+ib), ldda, c_one,
                                 dA(0,i), ldda);

                    magma_tally4_ssyrk( Magma_tally4Upper, Magma_tally4NoTrans, ib,(n-i-ib),
                                 d_one, dA(i, i+ib), ldda,
                                 d_one, dA(i, i),    ldda);
                }
            }
        }
        else {
            /* Compute the product L' * L. */
            for (i=0; i < n; i += nb) {
                ib=min(nb,(n-i));

                magma_tally4_strmm( Magma_tally4Left, Magma_tally4Lower,
                             Magma_tally4ConjTrans, Magma_tally4NonUnit, ib,
                             i, c_one, dA(i,i), ldda,
                             dA(i, 0),ldda);

                magma_tally4_sgetmatrix( ib, ib,
                                  dA(i, i), ldda,
                                  work,     ib );

                lapackf77_slauum(Magma_tally4LowerStr, &ib, work, &ib, info);

                magma_tally4_ssetmatrix( ib, ib,
                                  work,     ib,
                                  dA(i, i), ldda );

                if (i+ib < n) {
                    magma_tally4_sgemm( Magma_tally4ConjTrans, Magma_tally4NoTrans,
                                 ib, i, (n-i-ib), c_one, dA( i+ib,i),
                                 ldda, dA(i+ib, 0),ldda, c_one,
                                 dA(i,0), ldda);
                    magma_tally4_ssyrk( Magma_tally4Lower, Magma_tally4ConjTrans, ib, (n-i-ib),
                                 d_one, dA(i+ib, i), ldda,
                                 d_one, dA(i, i),    ldda);
                }
            }
        }
    }

    magma_tally4_queue_destroy( stream[0] );
    magma_tally4_queue_destroy( stream[1] );

    magma_tally4_free_pinned( work );

    return *info;
}