/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @generated from zpotrf_gpu.cpp normal z -> s, Fri Jan 30 19:00:12 2015
*/
#include "common_magma_tally4.h"

#define PRECISION_s

// === Define what BLAS to use ============================================
//#if (defined(PRECISION_s) || defined(PRECISION_d))
    #define magma_tally4_strsm magma_tally4blas_strsm
//#endif
// === End defining what BLAS to use =======================================

/**
    Purpose
    -------
    SPOTRF computes the Cholesky factorization of a real symmetric
    positive definite matrix dA.

    The factorization has the form
        dA = U**H * U,   if UPLO = Magma_tally4Upper, or
        dA = L  * L**H,  if UPLO = Magma_tally4Lower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of dA is stored;
      -     = Magma_tally4Lower:  Lower triangle of dA is stored.

    @param[in]
    n       INTEGER
            The order of the matrix dA.  N >= 0.

    @param[in,out]
    dA      REAL array on the GPU, dimension (LDDA,N)
            On entry, the symmetric matrix dA.  If UPLO = Magma_tally4Upper, the leading
            N-by-N upper triangular part of dA contains the upper
            triangular part of the matrix dA, and the strictly lower
            triangular part of dA is not referenced.  If UPLO = Magma_tally4Lower, the
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

    @ingroup magma_tally4_sposv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_spotrf_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info)
{
#define dA(i, j) (dA + (j)*ldda + (i))

    magma_tally4_int_t     j, jb, nb;
    const char* uplo_ = lapack_uplo_const_tally4( uplo );
    float c_one     = MAGMA_tally4_S_ONE;
    float c_neg_one = MAGMA_tally4_S_NEG_ONE;
    float *work;
    float          d_one     =  1.0;
    float          d_neg_one = -1.0;
    int upper = (uplo == Magma_tally4Upper);

    *info = 0;
    if (! upper && uplo != Magma_tally4Lower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    nb = magma_tally4_get_spotrf_nb(n);

    if (MAGMA_tally4_SUCCESS != magma_tally4_smalloc_pinned( &work, nb*nb )) {
        *info = MAGMA_tally4_ERR_HOST_ALLOC;
        return *info;
    }

    /* Define user stream if current stream is NULL */
    magma_tally4_queue_t stream[2];
    
    magma_tally4_queue_t orig_stream;
    magma_tally4blasGetKernelStream( &orig_stream );
    
    magma_tally4_queue_create( &stream[0] );
    if (orig_stream == NULL) {
        magma_tally4_queue_create( &stream[1] );
        magma_tally4blasSetKernelStream(stream[1]);
    }
    else {
        stream[1] = orig_stream;
    }
    
    if ((nb <= 1) || (nb >= n)) {
        /* Use unblocked code. */
        magma_tally4_sgetmatrix_async( n, n, dA, ldda, work, n, stream[1] );
        magma_tally4_queue_sync( stream[1] );
        lapackf77_spotrf(uplo_, &n, work, &n, info);
        magma_tally4_ssetmatrix_async( n, n, work, n, dA, ldda, stream[1] );
    }
    else {
        /* Use blocked code. */
        if (upper) {
            
            /* Compute the Cholesky factorization A = U'*U. */
            for (j=0; j < n; j += nb) {
                
                /* Update and factorize the current diagonal block and test
                   for non-positive-definiteness. Computing MIN */
                jb = min(nb, (n-j));
                
                magma_tally4_ssyrk(Magma_tally4Upper, Magma_tally4ConjTrans, jb, j,
                            d_neg_one, dA(0, j), ldda,
                            d_one,     dA(j, j), ldda);

                magma_tally4_queue_sync( stream[1] );
                magma_tally4_sgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[0] );
                
                if ( (j+jb) < n) {
                    /* Compute the current block row. */
                    magma_tally4_sgemm(Magma_tally4ConjTrans, Magma_tally4NoTrans,
                                jb, (n-j-jb), j,
                                c_neg_one, dA(0, j   ), ldda,
                                           dA(0, j+jb), ldda,
                                c_one,     dA(j, j+jb), ldda);
                }
                
                magma_tally4_queue_sync( stream[0] );
                lapackf77_spotrf(Magma_tally4UpperStr, &jb, work, &jb, info);
                magma_tally4_ssetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[1] );
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }

                if ( (j+jb) < n) {
                    magma_tally4_strsm( Magma_tally4Left, Magma_tally4Upper, Magma_tally4ConjTrans, Magma_tally4NonUnit,
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

                magma_tally4_ssyrk(Magma_tally4Lower, Magma_tally4NoTrans, jb, j,
                            d_neg_one, dA(j, 0), ldda,
                            d_one,     dA(j, j), ldda);
                
                magma_tally4_queue_sync( stream[1] );
                magma_tally4_sgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        work,     jb, stream[0] );
                
                if ( (j+jb) < n) {
                    magma_tally4_sgemm( Magma_tally4NoTrans, Magma_tally4ConjTrans,
                                 (n-j-jb), jb, j,
                                 c_neg_one, dA(j+jb, 0), ldda,
                                            dA(j,    0), ldda,
                                 c_one,     dA(j+jb, j), ldda);
                }

                magma_tally4_queue_sync( stream[0] );
                lapackf77_spotrf(Magma_tally4LowerStr, &jb, work, &jb, info);
                magma_tally4_ssetmatrix_async( jb, jb,
                                        work,     jb,
                                        dA(j, j), ldda, stream[1] );
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }
                
                if ( (j+jb) < n) {
                    magma_tally4_strsm(Magma_tally4Right, Magma_tally4Lower, Magma_tally4ConjTrans, Magma_tally4NonUnit,
                                (n-j-jb), jb,
                                c_one, dA(j,    j), ldda,
                                       dA(j+jb, j), ldda);
                }
            }
        }
    }

    magma_tally4_free_pinned( work );

    magma_tally4_queue_destroy( stream[0] );
    if (orig_stream == NULL) {
        magma_tally4_queue_destroy( stream[1] );
    }
    magma_tally4blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally4_spotrf_gpu */
