/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @precisions normal z -> s d c
*/
#include "common_magma_tally3.h"

#define PRECISION_z

// === Define what BLAS to use ============================================
//#if (defined(PRECISION_s) || defined(PRECISION_d))
    #define magma_tally3_ztrsm magma_tally3blas_ztrsm
//#endif
// === End defining what BLAS to use =======================================

/**
    Purpose
    -------
    ZPOTRF computes the Cholesky factorization of a complex Hermitian
    positive definite matrix dA.

    The factorization has the form
        dA = U**H * U,   if UPLO = Magma_tally3Upper, or
        dA = L  * L**H,  if UPLO = Magma_tally3Lower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  Upper triangle of dA is stored;
      -     = Magma_tally3Lower:  Lower triangle of dA is stored.

    @param[in]
    n       INTEGER
            The order of the matrix dA.  N >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the Hermitian matrix dA.  If UPLO = Magma_tally3Upper, the leading
            N-by-N upper triangular part of dA contains the upper
            triangular part of the matrix dA, and the strictly lower
            triangular part of dA is not referenced.  If UPLO = Magma_tally3Lower, the
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

    @ingroup magma_tally3_zposv_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zpotrf_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info)
{
#define dA(i, j) (dA + (j)*ldda + (i))

    magma_tally3_int_t     j, jb, nb;
    const char* uplo_ = lapack_uplo_const_tally3( uplo );
    magma_tally3DoubleComplex c_one     = MAGMA_tally3_Z_ONE;
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3DoubleComplex *work;
    double          d_one     =  1.0;
    double          d_neg_one = -1.0;
    int upper = (uplo == Magma_tally3Upper);

    *info = 0;
    if (! upper && uplo != Magma_tally3Lower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    nb = magma_tally3_get_zpotrf_nb(n);

    if (MAGMA_tally3_SUCCESS != magma_tally3_zmalloc_pinned( &work, nb*nb )) {
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        return *info;
    }

    /* Define user stream if current stream is NULL */
    magma_tally3_queue_t stream[2];
    
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );
    
    magma_tally3_queue_create( &stream[0] );
    if (orig_stream == NULL) {
        magma_tally3_queue_create( &stream[1] );
        magma_tally3blasSetKernelStream(stream[1]);
    }
    else {
        stream[1] = orig_stream;
    }
    
    if ((nb <= 1) || (nb >= n)) {
        /* Use unblocked code. */
        magma_tally3_zgetmatrix_async( n, n, dA, ldda, work, n, stream[1] );
        magma_tally3_queue_sync( stream[1] );
        lapackf77_zpotrf(uplo_, &n, work, &n, info);
        magma_tally3_zsetmatrix_async( n, n, work, n, dA, ldda, stream[1] );
    }
    else {
        /* Use blocked code. */
        if (upper) {
            
            /* Compute the Cholesky factorization A = U'*U. */
            for (j=0; j < n; j += nb) {
                
                /* Update and factorize the current diagonal block and test
                   for non-positive-definiteness. Computing MIN */
                jb = min(nb, (n-j));
                
                magma_tally3_zherk(Magma_tally3Upper, Magma_tally3ConjTrans, jb, j,
                            d_neg_one, dA(0, j), ldda,
                            d_one,     dA(j, j), ldda);

                magma_tally3_queue_sync( stream[1] );
                magma_tally3_zgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[0] );
                
                if ( (j+jb) < n) {
                    /* Compute the current block row. */
                    magma_tally3_zgemm(Magma_tally3ConjTrans, Magma_tally3NoTrans,
                                jb, (n-j-jb), j,
                                c_neg_one, dA(0, j   ), ldda,
                                           dA(0, j+jb), ldda,
                                c_one,     dA(j, j+jb), ldda);
                }
                
                magma_tally3_queue_sync( stream[0] );
                lapackf77_zpotrf(Magma_tally3UpperStr, &jb, work, &jb, info);
                magma_tally3_zsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[1] );
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }

                if ( (j+jb) < n) {
                    magma_tally3_ztrsm( Magma_tally3Left, Magma_tally3Upper, Magma_tally3ConjTrans, Magma_tally3NonUnit,
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

                magma_tally3_zherk(Magma_tally3Lower, Magma_tally3NoTrans, jb, j,
                            d_neg_one, dA(j, 0), ldda,
                            d_one,     dA(j, j), ldda);
                
                magma_tally3_queue_sync( stream[1] );
                magma_tally3_zgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[0] );
                
                if ( (j+jb) < n) {
                    magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans,
                                 (n-j-jb), jb, j,
                                 c_neg_one, dA(j+jb, 0), ldda,
                                            dA(j,    0), ldda,
                                 c_one,     dA(j+jb, j), ldda);
                }

                magma_tally3_queue_sync( stream[0] );
                lapackf77_zpotrf(Magma_tally3LowerStr, &jb, work, &jb, info);
                magma_tally3_zsetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[1] );
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }
                
                if ( (j+jb) < n) {
                    magma_tally3_ztrsm(Magma_tally3Right, Magma_tally3Lower, Magma_tally3ConjTrans, Magma_tally3NonUnit,
                                (n-j-jb), jb,
                                c_one, dA(j,    j), ldda,
                                       dA(j+jb, j), ldda);
                }
            }
        }
    }

    magma_tally3_free_pinned( work );

    magma_tally3_queue_destroy( stream[0] );
    if (orig_stream == NULL) {
        magma_tally3_queue_destroy( stream[1] );
    }
    magma_tally3blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally3_zpotrf_gpu */
