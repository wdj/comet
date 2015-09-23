/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Ichitaro Yamazaki                                                                   
       @author Stan Tomov

       @generated from zhetrf_nopiv.cpp normal z -> d, Fri Jan 30 19:00:16 2015
*/

#include "common_magma_minproduct.h"
#include "trace.h"

/**
    Purpose   
    =======   

    DSYTRF_nopiv computes the LDLt factorization of a real symmetric   
    matrix A. This version does not require work space on the GPU passed 
    as input. GPU memory is allocated in the routine.

    The factorization has the form   
       A = U^H * D * U  , if UPLO = 'U', or   
       A = L   * D * L^H, if UPLO = 'L',   
    where U is an upper triangular matrix, L is lower triangular, and
    D is a diagonal matrix.

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    ---------
    @param[in]
    UPLO    CHARACTER*1   
      -     = 'U':  Upper triangle of A is stored;   
      -     = 'L':  Lower triangle of A is stored.   

    @param[in]
    N       INTEGER   
            The order of the matrix A.  N >= 0.   

    @param[in,out]
    A       DOUBLE_PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U^H D U or A = L D L^H.   
    \n
            Higher performance is achieved if A is in pinned memory.

    @param[in]
    LDA     INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    @param[out]
    INFO    INTEGER   
      -     = 0:  successful exit   
      -     < 0:  if INFO = -i, the i-th argument had an illegal value 
                  if INFO = -6, the GPU memory allocation failed 
      -     > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.

    @ingroup magma_minproduct_dsysv_comp
    ******************************************************************* */
extern "C" magma_minproduct_int_t
magma_minproduct_dsytrf_nopiv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info)
{
    #define  A(i, j)  ( A +(j)*lda  + (i))
    #define dA(i, j)  (dA +(j)*ldda + (i))
    #define dW(i, j)  (dW +(j)*ldda + (i))
    #define dWt(i, j) (dW +(j)*nb   + (i))

    double zone  = MAGMA_minproduct_D_ONE;
    double mzone = MAGMA_minproduct_D_NEG_ONE;
    int                upper = (uplo == Magma_minproductUpper);
    magma_minproduct_int_t j, k, jb, ldda, nb, ib, iinfo;
    magma_minproductDouble_ptr dA;
    magma_minproductDouble_ptr dW;

    *info = 0;
    if (! upper && uplo != Magma_minproductLower) {
      *info = -1;
    } else if (n < 0) {
      *info = -2;
    } else if (lda < max(1,n)) {
      *info = -4;
    }
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return MAGMA_minproduct_ERR_ILLEGAL_VALUE;
    }

    /* Quick return */
    if ( n == 0 )
      return MAGMA_minproduct_SUCCESS;

    ldda = ((n+31)/32)*32;
    nb = magma_minproduct_get_dsytrf_nopiv_nb(n);
    ib = min(32, nb); // inner-block for diagonal factorization

    if ((MAGMA_minproduct_SUCCESS != magma_minproduct_dmalloc(&dA, n *ldda)) ||
        (MAGMA_minproduct_SUCCESS != magma_minproduct_dmalloc(&dW, nb*ldda))) {
        /* alloc failed so call the non-GPU-resident version */
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        return *info;
    }

    magma_minproduct_queue_t stream[2];
    magma_minproduct_event_t event;
    magma_minproduct_queue_create(&stream[0]);
    magma_minproduct_queue_create(&stream[1]);
    magma_minproduct_event_create( &event );
    trace_init( 1, 1, 2, stream );

    /* Use hybrid blocked code. */
    if (upper) {
        //=========================================================
        // Compute the LDLt factorization A = U'*D*U without pivoting.
        // copy matrix to GPU
        for (j=0; j<n; j+=nb) {
            jb = min(nb, (n-j));
            trace_gpu_start( 0, 0, "set", "set" );
            magma_minproduct_dsetmatrix_async(j+jb, jb, A(0, j), lda, dA(0, j), ldda, stream[0]);
            trace_gpu_end( 0, 0 );
        }
        
        // main loop
        for (j=0; j<n; j += nb) {
            jb = min(nb, (n-j));
            
            // copy A(j,j) back to CPU
            trace_gpu_start( 0, 0, "get", "get" );
            if ( j!=0) {
                //magma_minproduct_event_sync(event);
                magma_minproduct_dgetmatrix_async(jb, jb, dA(j, j), ldda, A(j,j), lda, stream[1]);
            }
            trace_gpu_end( 0, 0 );
            
            // factorize the diagonal block
            magma_minproduct_queue_sync(stream[1]);
            trace_cpu_start( 0, "potrf", "potrf" );
            dsytrf_nopiv_cpu(Magma_minproductUpper, jb, ib, A(j, j), lda, info);
            trace_cpu_end( 0 );
            if (*info != 0){
                *info = *info + j;
                break;
            }
            
            // copy A(j,j) back to GPU
            trace_gpu_start( 0, 0, "set", "set" );
            magma_minproduct_dsetmatrix_async(jb, jb, A(j, j), lda, dA(j, j), ldda, stream[0]);
            trace_gpu_end( 0, 0 );
            
            // copy j-th column of U back to CPU
            trace_gpu_start( 0, 1, "get", "get" );                                                                         
            magma_minproduct_dgetmatrix_async(j, jb, dA(0, j), ldda, A(0, j), lda, stream[1]);                                        
            trace_gpu_end( 0, 1 );

            if ( (j+jb) < n) {
                // compute the off-diagonal blocks of current block column
                magma_minproductblasSetKernelStream( stream[0] );
                trace_gpu_start( 0, 0, "trsm", "trsm" );
                magma_minproduct_dtrsm(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductConjTrans, Magma_minproductUnit, 
                            jb, (n-j-jb), 
                            zone, dA(j, j),    ldda, 
                            dA(j, j+jb), ldda);
                magma_minproduct_dcopymatrix( jb, n-j-jb, dA( j, j+jb ), ldda, dWt( 0, j+jb ), nb );
                
                // update the trailing submatrix with D
                magma_minproductblas_dlascl_diag(Magma_minproductUpper, jb, n-j-jb,
                                      dA(j,    j), ldda,
                                      dA(j, j+jb), ldda,
                                      &iinfo);
                trace_gpu_end( 0, 0 );
                
                // update the trailing submatrix with U and W
                trace_gpu_start( 0, 0, "gemm", "gemm" );
                for (k=j+jb; k<n; k+=nb) {
                    magma_minproduct_int_t kb = min(nb,n-k);
                    magma_minproduct_dgemm(Magma_minproductConjTrans, Magma_minproductNoTrans, kb, n-k, jb,
                                mzone, dWt(0, k), nb, 
                                       dA(j, k), ldda,
                                zone,  dA(k, k), ldda);
                    if (k==j+jb) {
                        // magma_minproduct_event_record( event, stream[0] );
                        magma_minproduct_queue_sync( stream[0] );
                    }
                }
                trace_gpu_end( 0, 0 );
            }
        }
    } else {
        //=========================================================
        // Compute the LDLt factorization A = L*D*L' without pivoting.
        // copy the matrix to GPU
        for (j=0; j<n; j+=nb) {
            jb = min(nb, (n-j));
            trace_gpu_start( 0, 0, "set", "set" );
            magma_minproduct_dsetmatrix_async((n-j), jb, A(j, j), lda, dA(j, j), ldda, stream[0]);
            trace_gpu_end( 0, 0 );
        }
        
        // main loop
        for (j=0; j<n; j+=nb) {
            jb = min(nb, (n-j));
            
            // copy A(j,j) back to CPU
            trace_gpu_start( 0, 0, "get", "get" );
            if (j!=0) {
                //magma_minproduct_event_sync(event);
                magma_minproduct_dgetmatrix_async(jb, jb, dA(j, j), ldda, A(j,j), lda, stream[1]);
            }
            trace_gpu_end( 0, 0 );
            
            // factorize the diagonal block
            magma_minproduct_queue_sync(stream[1]);
            trace_cpu_start( 0, "potrf", "potrf" );
            dsytrf_nopiv_cpu(Magma_minproductLower, jb, ib, A(j, j), lda, info);
            trace_cpu_end( 0 );
            if (*info != 0){
                *info = *info + j;
                break;
            }

            // copy A(j,j) back to GPU
            trace_gpu_start( 0, 0, "set", "set" );
            magma_minproduct_dsetmatrix_async(jb, jb, A(j, j), lda, dA(j, j), ldda, stream[0]);
            trace_gpu_end( 0, 0 );
            
            // copy j-th row of L back to CPU
            trace_gpu_start( 0, 1, "get", "get" );
            magma_minproduct_dgetmatrix_async(jb, j, dA(j, 0), ldda, A(j, 0), lda, stream[1]);
            trace_gpu_end( 0, 1 ); 
            
            if ( (j+jb) < n) {
                // compute the off-diagonal blocks of current block column
                magma_minproductblasSetKernelStream( stream[0] );
                trace_gpu_start( 0, 0, "trsm", "trsm" );
                magma_minproduct_dtrsm(Magma_minproductRight, Magma_minproductLower, Magma_minproductConjTrans, Magma_minproductUnit, 
                            (n-j-jb), jb, 
                            zone, dA(j,    j), ldda, 
                            dA(j+jb, j), ldda);
                magma_minproduct_dcopymatrix( n-j-jb,jb, dA( j+jb, j ), ldda, dW( j+jb, 0 ), ldda );
                
                // update the trailing submatrix with D
                magma_minproductblas_dlascl_diag(Magma_minproductLower, n-j-jb, jb,
                                      dA(j,    j), ldda,
                                      dA(j+jb, j), ldda,
                                      &iinfo);
                trace_gpu_end( 0, 0 );
                
                // update the trailing submatrix with L and W
                trace_gpu_start( 0, 0, "gemm", "gemm" );
                for (k=j+jb; k<n; k+=nb) {
                    magma_minproduct_int_t kb = min(nb,n-k);
                    magma_minproduct_dgemm(Magma_minproductNoTrans, Magma_minproductConjTrans, n-k, kb, jb,
                                mzone, dA(k, j), ldda, 
                                       dW(k, 0), ldda,
                                zone,  dA(k, k), ldda);
                    if (k==j+jb) {
                        //magma_minproduct_event_record( event, stream[0] );
                        magma_minproduct_queue_sync(stream[0]);
                    }
                }
                trace_gpu_end( 0, 0 );
            }
        }
    }
    
    trace_finalize( "dsytrf.svg","trace.css" );
    magma_minproduct_queue_destroy(stream[0]);
    magma_minproduct_queue_destroy(stream[1]);
    magma_minproduct_event_destroy( event );
    magma_minproduct_free(dW);
    magma_minproduct_free(dA);
    
    return MAGMA_minproduct_SUCCESS;
} /* magma_minproduct_dsytrf_nopiv */

