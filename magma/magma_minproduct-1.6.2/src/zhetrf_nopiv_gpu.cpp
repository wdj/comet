/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Ichitaro Yamazaki                                                                   
       @author Stan Tomov

       @precisions normal z -> s d c
*/
#include "common_magma_minproduct.h"
#include "trace.h"

/**
    Purpose   
    =======   

    ZHETRF_nopiv_gpu computes the LDLt factorization of a complex Hermitian   
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
    dA      COMPLEX_16 array on the GPU, dimension (LDA,N)   
            On entry, the Hermitian matrix A.  If UPLO = 'U', the leading   
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
    
    @ingroup magma_minproduct_zhesv_comp
    ******************************************************************* */
extern "C" magma_minproduct_int_t
magma_minproduct_zhetrf_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info)
{
    #define  A(i, j)  (A)
    #define dA(i, j)  (dA +(j)*ldda + (i))
    #define dW(i, j)  (dW +(j)*ldda + (i))
    #define dWt(i, j) (dW +(j)*nb   + (i))

    /* Local variables */
    magma_minproductDoubleComplex zone  = MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex mzone = MAGMA_minproduct_Z_NEG_ONE;
    int                upper = (uplo == Magma_minproductUpper);
    magma_minproduct_int_t j, k, jb, nb, ib, iinfo;

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
        return MAGMA_minproduct_ERR_ILLEGAL_VALUE;
    }

    /* Quick return */
    if ( n == 0 )
      return MAGMA_minproduct_SUCCESS;

    nb = magma_minproduct_get_zhetrf_nopiv_nb(n);
    ib = min(32, nb); // inner-block for diagonal factorization

    magma_minproduct_queue_t orig_stream;
    magma_minproductblasGetKernelStream( &orig_stream );


    magma_minproduct_queue_t stream[2];
    magma_minproduct_event_t event;
    magma_minproduct_queue_create(&stream[0]);
    magma_minproduct_queue_create(&stream[1]);
    magma_minproduct_event_create( &event );
    trace_init( 1, 1, 2, stream );

    // CPU workspace
    magma_minproductDoubleComplex *A;
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc_pinned( &A, nb*nb )) {
        *info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return *info;
    }

    // GPU workspace
    magma_minproductDoubleComplex_ptr dW;
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc( &dW, (1+nb)*ldda )) {
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
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
            //magma_minproduct_queue_wait_event( stream[1], event );                                                                
            magma_minproduct_event_sync(event);
            magma_minproduct_zgetmatrix_async(jb, jb, dA(j, j), ldda, A(j,j), nb, stream[1]);
            trace_gpu_end( 0, 0 );

            // factorize the diagonal block
            magma_minproduct_queue_sync(stream[1]);
            trace_cpu_start( 0, "potrf", "potrf" );
            zhetrf_nopiv_cpu(Magma_minproductUpper, jb, ib, A(j, j), nb, info);
            trace_cpu_end( 0 );
            if (*info != 0){
                *info = *info + j;
                break;
            }
            
            // copy A(j,j) back to GPU
            trace_gpu_start( 0, 0, "set", "set" );
            magma_minproduct_zsetmatrix_async(jb, jb, A(j, j), nb, dA(j, j), ldda, stream[0]);
            trace_gpu_end( 0, 0 );
                
            if ( (j+jb) < n) {
                // compute the off-diagonal blocks of current block column
                magma_minproductblasSetKernelStream( stream[0] );
                trace_gpu_start( 0, 0, "trsm", "trsm" );
                magma_minproduct_ztrsm(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductConjTrans, Magma_minproductUnit, 
                            jb, (n-j-jb), 
                            zone, dA(j, j),    ldda, 
                            dA(j, j+jb), ldda);
                magma_minproduct_zcopymatrix( jb, n-j-jb, dA( j, j+jb ), ldda, dWt( 0, j+jb ), nb );
                
                // update the trailing submatrix with D
                magma_minproductblas_zlascl_diag(Magma_minproductUpper, jb, n-j-jb,
                                      dA(j,    j), ldda,
                                      dA(j, j+jb), ldda,
                                      &iinfo);
                trace_gpu_end( 0, 0 );
                
                // update the trailing submatrix with U and W
                trace_gpu_start( 0, 0, "gemm", "gemm" );
                for (k=j+jb; k<n; k+=nb) {
                    magma_minproduct_int_t kb = min(nb,n-k);
                    magma_minproduct_zgemm(Magma_minproductConjTrans, Magma_minproductNoTrans, kb, n-k, jb,
                                mzone, dWt(0, k), nb, 
                                       dA(j, k), ldda,
                                zone,  dA(k, k), ldda);
                    if (k==j+jb)
                        magma_minproduct_event_record( event, stream[0] );
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
            //magma_minproduct_queue_wait_event( stream[0], event );                                                                
            magma_minproduct_event_sync(event);
            magma_minproduct_zgetmatrix_async(jb, jb, dA(j, j), ldda, A(j,j), nb, stream[1]);
            trace_gpu_end( 0, 0 );
            
            // factorize the diagonal block
            magma_minproduct_queue_sync(stream[1]);
            trace_cpu_start( 0, "potrf", "potrf" );
            zhetrf_nopiv_cpu(Magma_minproductLower, jb, ib, A(j, j), nb, info);
            trace_cpu_end( 0 );
            if (*info != 0){
                *info = *info + j;
                break;
            }

            // copy A(j,j) back to GPU
            trace_gpu_start( 0, 0, "set", "set" );
            magma_minproduct_zsetmatrix_async(jb, jb, A(j, j), nb, dA(j, j), ldda, stream[0]);
            trace_gpu_end( 0, 0 );
            
            if ( (j+jb) < n) {
                // compute the off-diagonal blocks of current block column
                magma_minproductblasSetKernelStream( stream[0] );
                trace_gpu_start( 0, 0, "trsm", "trsm" );
                magma_minproduct_ztrsm(Magma_minproductRight, Magma_minproductLower, Magma_minproductConjTrans, Magma_minproductUnit, 
                            (n-j-jb), jb, 
                            zone, dA(j,    j), ldda, 
                            dA(j+jb, j), ldda);
                magma_minproduct_zcopymatrix( n-j-jb,jb, dA( j+jb, j ), ldda, dW( j+jb, 0 ), ldda );
                
                // update the trailing submatrix with D
                magma_minproductblas_zlascl_diag(Magma_minproductLower, n-j-jb, jb,
                                      dA(j,    j), ldda,
                                      dA(j+jb, j), ldda,
                                      &iinfo);
                trace_gpu_end( 0, 0 );
                
                // update the trailing submatrix with L and W
                trace_gpu_start( 0, 0, "gemm", "gemm" );
                for (k=j+jb; k<n; k+=nb) {
                    magma_minproduct_int_t kb = min(nb,n-k);
                    magma_minproduct_zgemm(Magma_minproductNoTrans, Magma_minproductConjTrans, n-k, kb, jb,
                                mzone, dA(k, j), ldda, 
                                       dW(k, 0), ldda,
                                zone,  dA(k, k), ldda);
                    if (k==j+jb)
                        magma_minproduct_event_record( event, stream[0] );
                }
                trace_gpu_end( 0, 0 );
            }
        }
    }
    
    trace_finalize( "zhetrf.svg","trace.css" );
    magma_minproduct_queue_destroy(stream[0]);
    magma_minproduct_queue_destroy(stream[1]);
    magma_minproduct_event_destroy( event );
    magma_minproduct_free( dW );
    magma_minproduct_free_pinned( A );
    
    magma_minproductblasSetKernelStream( orig_stream );
    return MAGMA_minproduct_SUCCESS;
} /* magma_minproduct_zhetrf_nopiv */

