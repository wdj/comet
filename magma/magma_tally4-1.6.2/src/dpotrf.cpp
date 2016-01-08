/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @generated from zpotrf.cpp normal z -> d, Fri Jan 30 19:00:14 2015
*/
#include "common_magma_tally4.h"

#define PRECISION_d

// === Define what BLAS to use ============================================
//#if defined(PRECISION_s) || defined(PRECISION_d)
    #define magma_tally4_dtrsm magma_tally4blas_dtrsm
//#endif
// === End defining what BLAS to use ======================================

/**
    Purpose
    -------
    DPOTRF computes the Cholesky factorization of a real symmetric
    positive definite matrix A. This version does not require work
    space on the GPU passed as input. GPU memory is allocated in the
    routine.

    The factorization has the form
        A = U**H * U,  if uplo = Magma_tally4Upper, or
        A = L  * L**H, if uplo = Magma_tally4Lower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.
    
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of A is stored;
      -     = Magma_tally4Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If uplo = Magma_tally4Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If uplo = Magma_tally4Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization A = U**H * U or A = L * L**H.
    \n
            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using magma_tally4_malloc_pinned.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
      -     > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

    @ingroup magma_tally4_dposv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dpotrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info)
{
#define  A(i_, j_)  (A + (j_)*lda  + (i_))
#define dA(i_, j_) (dA + (j_)*ldda + (i_))

    /* Local variables */
    const char* uplo_ = lapack_uplo_const( uplo );
    magma_tally4_int_t        ldda, nb;
    magma_tally4_int_t j, jb;
    double    c_one     = MAGMA_tally4_D_ONE;
    double    c_neg_one = MAGMA_tally4_D_NEG_ONE;
    magma_tally4Double_ptr dA;
    double             d_one     =  1.0;
    double             d_neg_one = -1.0;
    int upper = (uplo == Magma_tally4Upper);

    *info = 0;
    if (! upper && uplo != Magma_tally4Lower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ( n == 0 )
        return *info;

    magma_tally4_int_t ngpu = magma_tally4_num_gpus();
    if ( ngpu > 1 ) {
        /* call multiple-GPU interface  */
        return magma_tally4_dpotrf_m(ngpu, uplo, n, A, lda, info);
    }

    ldda = ((n+31)/32)*32;
    
    if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc( &dA, (n)*ldda )) {
        /* alloc failed so call the non-GPU-resident version */
        return magma_tally4_dpotrf_m(ngpu, uplo, n, A, lda, info);
    }

    /* Define user stream if current stream is NULL */
    magma_tally4_queue_t stream[3];
    
    magma_tally4_queue_t orig_stream;
    magma_tally4blasGetKernelStream( &orig_stream );

    magma_tally4_queue_create( &stream[0] );
    magma_tally4_queue_create( &stream[2] );

    if (orig_stream == NULL) {
        magma_tally4_queue_create( &stream[1] );
        magma_tally4blasSetKernelStream(stream[1]);
    }
    else {
        stream[1] = orig_stream;
    }

    nb = magma_tally4_get_dpotrf_nb(n);

    if (nb <= 1 || nb >= n) {
        lapackf77_dpotrf(uplo_, &n, A, &lda, info);
    } else {
        /* Use hybrid blocked code. */
        if (upper) {
            /* Compute the Cholesky factorization A = U'*U. */
            for (j=0; j < n; j += nb) {
                /* Update and factorize the current diagonal block and test
                   for non-positive-definiteness. Computing MIN */
                jb = min(nb, (n-j));
                magma_tally4_dsetmatrix_async( jb, (n-j), A(j, j), lda, dA(j, j), ldda, stream[1]);
                
                magma_tally4_dsyrk(Magma_tally4Upper, Magma_tally4ConjTrans, jb, j,
                            d_neg_one, dA(0, j), ldda,
                            d_one,     dA(j, j), ldda);
                magma_tally4_queue_sync( stream[1] );

                magma_tally4_dgetmatrix_async( jb, jb,
                                        dA(j, j), ldda,
                                        A(j, j),  lda, stream[0] );
                
                if ( (j+jb) < n) {
                    magma_tally4_dgemm(Magma_tally4ConjTrans, Magma_tally4NoTrans,
                                jb, (n-j-jb), j,
                                c_neg_one, dA(0, j   ), ldda,
                                           dA(0, j+jb), ldda,
                                c_one,     dA(j, j+jb), ldda);
                }
                
                magma_tally4_dgetmatrix_async( j, jb,
                                        dA(0, j), ldda,
                                        A (0, j),  lda, stream[2] );

                magma_tally4_queue_sync( stream[0] );
                lapackf77_dpotrf(Magma_tally4UpperStr, &jb, A(j, j), &lda, info);
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }
                magma_tally4_dsetmatrix_async( jb, jb,
                                        A(j, j),  lda,
                                        dA(j, j), ldda, stream[0] );
                magma_tally4_queue_sync( stream[0] );

                if ( (j+jb) < n ) {
                    magma_tally4_dtrsm(Magma_tally4Left, Magma_tally4Upper, Magma_tally4ConjTrans, Magma_tally4NonUnit,
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
                magma_tally4_dsetmatrix_async( (n-j), jb, A(j, j), lda, dA(j, j), ldda, stream[1]);

                magma_tally4_dsyrk(Magma_tally4Lower, Magma_tally4NoTrans, jb, j,
                            d_neg_one, dA(j, 0), ldda,
                            d_one,     dA(j, j), ldda);
                magma_tally4_queue_sync( stream[1] );

                magma_tally4_dgetmatrix_async( jb, jb,
                                        dA(j,j), ldda,
                                        A(j,j),  lda, stream[0] );

                if ( (j+jb) < n) {
                    magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4ConjTrans,
                                 (n-j-jb), jb, j,
                                 c_neg_one, dA(j+jb, 0), ldda,
                                            dA(j,    0), ldda,
                                 c_one,     dA(j+jb, j), ldda);
                }
                
                magma_tally4_dgetmatrix_async( jb, j,
                                        dA(j, 0), ldda,
                                        A(j, 0),  lda, stream[2] );

                magma_tally4_queue_sync( stream[0] );
                lapackf77_dpotrf(Magma_tally4LowerStr, &jb, A(j, j), &lda, info);
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }
                magma_tally4_dsetmatrix_async( jb, jb,
                                        A(j, j),  lda,
                                        dA(j, j), ldda, stream[0] );
                magma_tally4_queue_sync( stream[0] );

                if ( (j+jb) < n) {
                    magma_tally4_dtrsm(Magma_tally4Right, Magma_tally4Lower, Magma_tally4ConjTrans, Magma_tally4NonUnit,
                                (n-j-jb), jb,
                                c_one, dA(j,    j), ldda,
                                       dA(j+jb, j), ldda);
                }
            }
        }
    }
    
    magma_tally4_queue_destroy( stream[0] );
    magma_tally4_queue_destroy( stream[2] );
    if (orig_stream == NULL) {
        magma_tally4_queue_destroy( stream[1] );
    }
    magma_tally4blasSetKernelStream( orig_stream );

    magma_tally4_free( dA );
    
    return *info;
} /* magma_tally4_dpotrf */
