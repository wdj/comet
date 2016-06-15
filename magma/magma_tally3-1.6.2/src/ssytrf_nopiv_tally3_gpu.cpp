/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Ichitaro Yamazaki                                                                   
       @author Stan Tomov

       @generated from zhetrf_nopiv_tally3_gpu.cpp normal z -> s, Fri Jan 30 19:00:17 2015
*/
#include "common_magma_tally3.h"
#include "trace.h"

/**
    Purpose   
    =======   

    SSYTRF_nopiv_gpu computes the LDLt factorization of a real symmetric   
    matrix A.

    The factorization has the form   
       A = U^H * D * U , if UPLO = 'U', or   
       A = L  * D * L^H, if UPLO = 'L',   
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
    dA      REAL array on the GPU, dimension (LDA,N)   
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
            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using cudaMallocHost.

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
    
    @ingroup magma_tally3_ssysv_comp
    ******************************************************************* */
extern "C" magma_tally3_int_t
magma_tally3_ssytrf_nopiv_tally3_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info)
{
    #define  A(i, j)  (A)
    #define dA(i, j)  (dA +(j)*ldda + (i))
    #define dW(i, j)  (dW +(j)*ldda + (i))
    #define dWt(i, j) (dW +(j)*nb   + (i))

    /* Local variables */
    float zone  = MAGMA_tally3_S_ONE;
    float mzone = MAGMA_tally3_S_NEG_ONE;
    int                upper = (uplo == Magma_tally3Upper);
    magma_tally3_int_t j, k, jb, nb, ib, iinfo;

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
        return MAGMA_tally3_ERR_ILLEGAL_VALUE;
    }

    /* Quick return */
    if ( n == 0 )
      return MAGMA_tally3_SUCCESS;

    nb = magma_tally3_get_ssytrf_nopiv_tally3_nb(n);
    ib = min(32, nb); // inner-block for diagonal factorization

    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );


    magma_tally3_queue_t stream[2];
    magma_tally3_event_t event;
    magma_tally3_queue_create(&stream[0]);
    magma_tally3_queue_create(&stream[1]);
    magma_tally3_event_create( &event );
    trace_init( 1, 1, 2, stream );

    // CPU workspace
    float *A;
    if (MAGMA_tally3_SUCCESS != magma_tally3_smalloc_pinned( &A, nb*nb )) {
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        return *info;
    }

    // GPU workspace
    magma_tally3Float_ptr dW;
    if (MAGMA_tally3_SUCCESS != magma_tally3_smalloc( &dW, (1+nb)*ldda )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }

    /* Use hybrid blocked code. */
    if (upper) {
        //=========================================================
        // Compute the LDLt factorization A = U'*D*U without pivoting.
        // main loop
        for (j=0; j<n; j += nb) {
            jb = min(nb, (n-j));
            
            // copy A(j,j) back to CPU
            trace_gpu_start( 0, 0, "get", "get" );
            //magma_tally3_queue_wait_event( stream[1], event );                                                                
            magma_tally3_event_sync(event);
            magma_tally3_sgetmatrix_async(jb, jb, dA(j, j), ldda, A(j,j), nb, stream[1]);
            trace_gpu_end( 0, 0 );

            // factorize the diagonal block
            magma_tally3_queue_sync(stream[1]);
            trace_cpu_start( 0, "potrf", "potrf" );
            ssytrf_nopiv_tally3_cpu(Magma_tally3Upper, jb, ib, A(j, j), nb, info);
            trace_cpu_end( 0 );
            if (*info != 0){
                *info = *info + j;
                break;
            }
            
            // copy A(j,j) back to GPU
            trace_gpu_start( 0, 0, "set", "set" );
            magma_tally3_ssetmatrix_async(jb, jb, A(j, j), nb, dA(j, j), ldda, stream[0]);
            trace_gpu_end( 0, 0 );
                
            if ( (j+jb) < n) {
                // compute the off-diagonal blocks of current block column
                magma_tally3blasSetKernelStream( stream[0] );
                trace_gpu_start( 0, 0, "trsm", "trsm" );
                magma_tally3_strsm(Magma_tally3Left, Magma_tally3Upper, Magma_tally3ConjTrans, Magma_tally3Unit, 
                            jb, (n-j-jb), 
                            zone, dA(j, j),    ldda, 
                            dA(j, j+jb), ldda);
                magma_tally3_scopymatrix( jb, n-j-jb, dA( j, j+jb ), ldda, dWt( 0, j+jb ), nb );
                
                // update the trailing submatrix with D
                magma_tally3blas_slascl_diag(Magma_tally3Upper, jb, n-j-jb,
                                      dA(j,    j), ldda,
                                      dA(j, j+jb), ldda,
                                      &iinfo);
                trace_gpu_end( 0, 0 );
                
                // update the trailing submatrix with U and W
                trace_gpu_start( 0, 0, "gemm", "gemm" );
                for (k=j+jb; k<n; k+=nb) {
                    magma_tally3_int_t kb = min(nb,n-k);
                    magma_tally3_sgemm(Magma_tally3ConjTrans, Magma_tally3NoTrans, kb, n-k, jb,
                                mzone, dWt(0, k), nb, 
                                       dA(j, k), ldda,
                                zone,  dA(k, k), ldda);
                    if (k==j+jb)
                        magma_tally3_event_record( event, stream[0] );
                }
                trace_gpu_end( 0, 0 );
            }
        }
    } else {
        //=========================================================
        // Compute the LDLt factorization A = L*D*L' without pivoting.
        // main loop
        for (j=0; j<n; j+=nb) {
            jb = min(nb, (n-j));
            
            // copy A(j,j) back to CPU
            trace_gpu_start( 0, 0, "get", "get" );
            //magma_tally3_queue_wait_event( stream[0], event );                                                                
            magma_tally3_event_sync(event);
            magma_tally3_sgetmatrix_async(jb, jb, dA(j, j), ldda, A(j,j), nb, stream[1]);
            trace_gpu_end( 0, 0 );
            
            // factorize the diagonal block
            magma_tally3_queue_sync(stream[1]);
            trace_cpu_start( 0, "potrf", "potrf" );
            ssytrf_nopiv_tally3_cpu(Magma_tally3Lower, jb, ib, A(j, j), nb, info);
            trace_cpu_end( 0 );
            if (*info != 0){
                *info = *info + j;
                break;
            }

            // copy A(j,j) back to GPU
            trace_gpu_start( 0, 0, "set", "set" );
            magma_tally3_ssetmatrix_async(jb, jb, A(j, j), nb, dA(j, j), ldda, stream[0]);
            trace_gpu_end( 0, 0 );
            
            if ( (j+jb) < n) {
                // compute the off-diagonal blocks of current block column
                magma_tally3blasSetKernelStream( stream[0] );
                trace_gpu_start( 0, 0, "trsm", "trsm" );
                magma_tally3_strsm(Magma_tally3Right, Magma_tally3Lower, Magma_tally3ConjTrans, Magma_tally3Unit, 
                            (n-j-jb), jb, 
                            zone, dA(j,    j), ldda, 
                            dA(j+jb, j), ldda);
                magma_tally3_scopymatrix( n-j-jb,jb, dA( j+jb, j ), ldda, dW( j+jb, 0 ), ldda );
                
                // update the trailing submatrix with D
                magma_tally3blas_slascl_diag(Magma_tally3Lower, n-j-jb, jb,
                                      dA(j,    j), ldda,
                                      dA(j+jb, j), ldda,
                                      &iinfo);
                trace_gpu_end( 0, 0 );
                
                // update the trailing submatrix with L and W
                trace_gpu_start( 0, 0, "gemm", "gemm" );
                for (k=j+jb; k<n; k+=nb) {
                    magma_tally3_int_t kb = min(nb,n-k);
                    magma_tally3_sgemm(Magma_tally3NoTrans, Magma_tally3ConjTrans, n-k, kb, jb,
                                mzone, dA(k, j), ldda, 
                                       dW(k, 0), ldda,
                                zone,  dA(k, k), ldda);
                    if (k==j+jb)
                        magma_tally3_event_record( event, stream[0] );
                }
                trace_gpu_end( 0, 0 );
            }
        }
    }
    
    trace_finalize( "ssytrf.svg","trace.css" );
    magma_tally3_queue_destroy(stream[0]);
    magma_tally3_queue_destroy(stream[1]);
    magma_tally3_event_destroy( event );
    magma_tally3_free( dW );
    magma_tally3_free_pinned( A );
    
    magma_tally3blasSetKernelStream( orig_stream );
    return MAGMA_tally3_SUCCESS;
} /* magma_tally3_ssytrf_nopiv_tally3 */
