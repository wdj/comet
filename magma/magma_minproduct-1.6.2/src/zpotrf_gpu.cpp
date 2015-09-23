/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @precisions normal z -> s d c
*/
#include "common_magma_minproduct.h"

#define PRECISION_z

// === Define what BLAS to use ============================================
//#if (defined(PRECISION_s) || defined(PRECISION_d))
    #define magma_minproduct_ztrsm magma_minproductblas_ztrsm
//#endif
// === End defining what BLAS to use =======================================

/**
    Purpose
    -------
    ZPOTRF computes the Cholesky factorization of a complex Hermitian
    positive definite matrix dA.

    The factorization has the form
        dA = U**H * U,   if UPLO = Magma_minproductUpper, or
        dA = L  * L**H,  if UPLO = Magma_minproductLower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  Upper triangle of dA is stored;
      -     = Magma_minproductLower:  Lower triangle of dA is stored.

    @param[in]
    n       INTEGER
            The order of the matrix dA.  N >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the Hermitian matrix dA.  If UPLO = Magma_minproductUpper, the leading
            N-by-N upper triangular part of dA contains the upper
            triangular part of the matrix dA, and the strictly lower
            triangular part of dA is not referenced.  If UPLO = Magma_minproductLower, the
            leading N-by-N lower triangular part of dA contains the lower
            triangular part of the matrix dA, and the strictly upper
            triangular part of dA is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization dA = U**H * U or dA = L * L**H.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).
            To benefit from coalescent memory accesses LDDA must be
            divisible by 16.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

    @ingroup magma_minproduct_zposv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_zpotrf_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info)
{
#define dA(i, j) (dA + (j)*ldda + (i))

    magma_minproduct_int_t     j, jb, nb;
    const char* uplo_ = lapack_uplo_const( uplo );
    magma_minproductDoubleComplex c_one     = MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex c_neg_one = MAGMA_minproduct_Z_NEG_ONE;
    magma_minproductDoubleComplex *work;
    double          d_one     =  1.0;
    double          d_neg_one = -1.0;
    int upper = (uplo == Magma_minproductUpper);

    *info = 0;
    if (! upper && uplo != Magma_minproductLower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    nb = magma_minproduct_get_zpotrf_nb(n);

    if (MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc_pinned( &work, nb*nb )) {
        *info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return *info;
    }

    /* Define user stream if current stream is NULL */
    magma_minproduct_queue_t stream[2];
    
    magma_minproduct_queue_t orig_stream;
    magma_minproductblasGetKernelStream( &orig_stream );
    
    magma_minproduct_queue_create( &stream[0] );
    if (orig_stream == NULL) {
        magma_minproduct_queue_create( &stream[1] );
        magma_minproductblasSetKernelStream(stream[1]);
    }
    else {
        stream[1] = orig_stream;
    }
    
    if ((nb <= 1) || (nb >= n)) {
        /* Use unblocked code. */
        magma_minproduct_zgetmatrix_async( n, n, dA, ldda, work, n, stream[1] );
        magma_minproduct_queue_sync( stream[1] );
        lapackf77_zpotrf(uplo_, &n, work, &n, info);
        magma_minproduct_zsetmatrix_async( n, n, work, n, dA, ldda, stream[1] );
    }
    else {
        /* Use blocked code. */
        if (upper) {
            
            /* Compute the Cholesky factorization A = U'*U. */
            for (j=0; j < n; j += nb) {
                
                /* Update and factorize the current diagonal block and test
                   for non-positive-definiteness. Computing MIN */
                jb = min(nb, (n-j));
                
                magma_minproduct_zherk(Magma_minproductUpper, Magma_minproductConjTrans, jb, j,
                            d_neg_one, dA(0, j), ldda,
                            d_one,     dA(j, j), ldda);

                magma_minproduct_queue_sync( stream[1] );
                magma_minproduct_zgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[0] );
                
                if ( (j+jb) < n) {
                    /* Compute the current block row. */
                    magma_minproduct_zgemm(Magma_minproductConjTrans, Magma_minproductNoTrans,
                                jb, (n-j-jb), j,
                                c_neg_one, dA(0, j   ), ldda,
                                           dA(0, j+jb), ldda,
                                c_one,     dA(j, j+jb), ldda);
                }
                
                magma_minproduct_queue_sync( stream[0] );
                lapackf77_zpotrf(Magma_minproductUpperStr, &jb, work, &jb, info);
                magma_minproduct_zsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[1] );
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }

                if ( (j+jb) < n) {
                    magma_minproduct_ztrsm( Magma_minproductLeft, Magma_minproductUpper, Magma_minproductConjTrans, Magma_minproductNonUnit,
                                 jb, (n-j-jb),
                                 c_one, dA(j, j   ), ldda,
                                        dA(j, j+jb), ldda);
                }
            }
        }
        else {
            //=========================================================
            // Compute the Cholesky factorization A = L*L'.
            for (j=0; j < n; j += nb) {
                //  Update and factorize the current diagonal block and test
                //  for non-positive-definiteness. Computing MIN
                jb = min(nb, (n-j));

                magma_minproduct_zherk(Magma_minproductLower, Magma_minproductNoTrans, jb, j,
                            d_neg_one, dA(j, 0), ldda,
                            d_one,     dA(j, j), ldda);
                
                magma_minproduct_queue_sync( stream[1] );
                magma_minproduct_zgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[0] );
                
                if ( (j+jb) < n) {
                    magma_minproduct_zgemm( Magma_minproductNoTrans, Magma_minproductConjTrans,
                                 (n-j-jb), jb, j,
                                 c_neg_one, dA(j+jb, 0), ldda,
                                            dA(j,    0), ldda,
                                 c_one,     dA(j+jb, j), ldda);
                }

                magma_minproduct_queue_sync( stream[0] );
                lapackf77_zpotrf(Magma_minproductLowerStr, &jb, work, &jb, info);
                magma_minproduct_zsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[1] );
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }
                
                if ( (j+jb) < n) {
                    magma_minproduct_ztrsm(Magma_minproductRight, Magma_minproductLower, Magma_minproductConjTrans, Magma_minproductNonUnit,
                                (n-j-jb), jb,
                                c_one, dA(j,    j), ldda,
                                       dA(j+jb, j), ldda);
                }
            }
        }
    }

    magma_minproduct_free_pinned( work );

    magma_minproduct_queue_destroy( stream[0] );
    if (orig_stream == NULL) {
        magma_minproduct_queue_destroy( stream[1] );
    }
    magma_minproductblasSetKernelStream( orig_stream );

    return *info;
} /* magma_minproduct_zpotrf_gpu */
