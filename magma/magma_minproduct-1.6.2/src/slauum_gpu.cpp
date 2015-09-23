/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlauum_gpu.cpp normal z -> s, Fri Jan 30 19:00:12 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    SLAUUM computes the product U * U' or L' * L, where the triangular
    factor U or L is stored in the upper or lower triangular part of
    the array dA.

    If UPLO = Magma_minproductUpper then the upper triangle of the result is stored,
    overwriting the factor U in dA.
    If UPLO = Magma_minproductLower then the lower triangle of the result is stored,
    overwriting the factor L in dA.
    This is the blocked form of the algorithm, calling Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
            Specifies whether the triangular factor stored in the array dA
            is upper or lower triangular:
      -     = Magma_minproductUpper:  Upper triangular
      -     = Magma_minproductLower:  Lower triangular

    @param[in]
    n       INTEGER
            The order of the triangular factor U or L.  N >= 0.

    @param[in,out]
    dA      REAL array on the GPU, dimension (LDDA,N)
            On entry, the triangular factor U or L.
            On exit, if UPLO = Magma_minproductUpper, the upper triangle of dA is
            overwritten with the upper triangle of the product U * U';
            if UPLO = Magma_minproductLower, the lower triangle of dA is overwritten with
            the lower triangle of the product L' * L.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -k, the k-th argument had an illegal value

    @ingroup magma_minproduct_sposv_aux
    ***************************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_slauum_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info)
{
#define dA(i, j) (dA + (j)*ldda + (i))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const( uplo );
    magma_minproduct_int_t         nb, i, ib;
    float              d_one = MAGMA_minproduct_D_ONE;
    float  c_one = MAGMA_minproduct_S_ONE;
    float  *work;

    int upper  = (uplo == Magma_minproductUpper);

    *info = 0;

    if (! upper && uplo != Magma_minproductLower)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,n))
        *info = -4;

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    nb = magma_minproduct_get_spotrf_nb(n);

    if (MAGMA_minproduct_SUCCESS != magma_minproduct_smalloc_pinned( &work, nb*nb )) {
        *info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return *info;
    }

    magma_minproduct_queue_t stream[2];
    magma_minproduct_queue_create( &stream[0] );
    magma_minproduct_queue_create( &stream[1] );

    if (nb <= 1 || nb >= n) {
        magma_minproduct_sgetmatrix( n, n, dA, ldda, work, n );
        lapackf77_slauum(uplo_, &n, work, &n, info);
        magma_minproduct_ssetmatrix( n, n, work, n, dA, ldda );
    }
    else {
        if (upper) {
            /* Compute inverse of upper triangular matrix */
            for (i=0; i < n; i += nb) {
                ib = min(nb, (n-i));

                /* Compute the product U * U'. */
                magma_minproduct_strmm( Magma_minproductRight, Magma_minproductUpper,
                             Magma_minproductConjTrans, Magma_minproductNonUnit, i, ib,
                             c_one, dA(i,i), ldda, dA(0, i),ldda);

                magma_minproduct_sgetmatrix( ib, ib,
                                  dA(i, i), ldda,
                                  work,     ib );

                lapackf77_slauum(Magma_minproductUpperStr, &ib, work, &ib, info);

                magma_minproduct_ssetmatrix( ib, ib,
                                  work,     ib,
                                  dA(i, i), ldda );

                if (i+ib < n) {
                    magma_minproduct_sgemm( Magma_minproductNoTrans, Magma_minproductConjTrans,
                                 i, ib, (n-i-ib), c_one, dA(0,i+ib),
                                 ldda, dA(i, i+ib), ldda, c_one,
                                 dA(0,i), ldda);

                    magma_minproduct_ssyrk( Magma_minproductUpper, Magma_minproductNoTrans, ib,(n-i-ib),
                                 d_one, dA(i, i+ib), ldda,
                                 d_one, dA(i, i),    ldda);
                }
            }
        }
        else {
            /* Compute the product L' * L. */
            for (i=0; i < n; i += nb) {
                ib=min(nb,(n-i));

                magma_minproduct_strmm( Magma_minproductLeft, Magma_minproductLower,
                             Magma_minproductConjTrans, Magma_minproductNonUnit, ib,
                             i, c_one, dA(i,i), ldda,
                             dA(i, 0),ldda);

                magma_minproduct_sgetmatrix( ib, ib,
                                  dA(i, i), ldda,
                                  work,     ib );

                lapackf77_slauum(Magma_minproductLowerStr, &ib, work, &ib, info);

                magma_minproduct_ssetmatrix( ib, ib,
                                  work,     ib,
                                  dA(i, i), ldda );

                if (i+ib < n) {
                    magma_minproduct_sgemm( Magma_minproductConjTrans, Magma_minproductNoTrans,
                                 ib, i, (n-i-ib), c_one, dA( i+ib,i),
                                 ldda, dA(i+ib, 0),ldda, c_one,
                                 dA(i,0), ldda);
                    magma_minproduct_ssyrk( Magma_minproductLower, Magma_minproductConjTrans, ib, (n-i-ib),
                                 d_one, dA(i+ib, i), ldda,
                                 d_one, dA(i, i),    ldda);
                }
            }
        }
    }

    magma_minproduct_queue_destroy( stream[0] );
    magma_minproduct_queue_destroy( stream[1] );

    magma_minproduct_free_pinned( work );

    return *info;
}
