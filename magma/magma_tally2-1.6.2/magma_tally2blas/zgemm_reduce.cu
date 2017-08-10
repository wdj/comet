/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"
#include "magma_tally2_templates.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
// size of work for a thread block
#define BLK_M 16
#define BLK_N 16

// BLK_K gets defined in magma_tally2blas_zgemm_reduce,
// because it depends on the CUDA architecture at runtime.


//==============================================================================
// BLK_K size is templated, as it depends on CUDA architecture at runtime.
// Hmm... how to compile for both CUDA arch 1.x and 2.x?

template< int BLK_K >
__global__
void zgemm_reduce_kernel(
    int m, int n, int k,
    magma_tally2DoubleComplex alpha,
    const magma_tally2DoubleComplex* __restrict__ dA, int lda,
    const magma_tally2DoubleComplex* __restrict__ dB, int ldb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex      * __restrict__ dC, int ldc)
{
#if (__CUDA_ARCH__ >= 200)
    const int tx = threadIdx.x;
    
    if (blockIdx.x*BLK_M + threadIdx.y < m && blockIdx.y*BLK_N + threadIdx.z < n){
    
        dA += (blockIdx.x*BLK_M + threadIdx.y) * lda;
        dB += (blockIdx.y*BLK_N + threadIdx.z) * ldb;
        dC +=  blockIdx.x*BLK_M + blockIdx.y*BLK_N * ldc;
        
        // was: sum[BLK_M][BLK_N+1][BLK_K+1];
        // moved 3rd dimension to 1st dimension to make magma_tally2_sum_reduce_3d interface nicer.
        __shared__ magma_tally2DoubleComplex sum[BLK_K][BLK_M+1][BLK_N+1];
        magma_tally2DoubleComplex lsum;
        
        /*  w := v**H * C  */
        lsum = MAGMA_tally2_Z_ZERO;
        for( int j = tx; j < k; j += BLK_K )
            lsum += MAGMA_tally2_Z_CNJG( dA[j] )* dB[j];
        
        sum[tx][threadIdx.y][threadIdx.z] = lsum;
        magma_tally2_sum_reduce_3d< BLK_K, BLK_M+1, BLK_N+1 >( tx, threadIdx.y, threadIdx.z, sum );
        
        /*  C := C - v * w  */
        __syncthreads();
        if (threadIdx.x == 0) {
            if (MAGMA_tally2_Z_EQUAL(beta, MAGMA_tally2_Z_ZERO))
                dC[threadIdx.y + threadIdx.z*ldc] = alpha*sum[0][threadIdx.y][threadIdx.z];
            else
                dC[threadIdx.y + threadIdx.z*ldc] = beta* dC[threadIdx.y + threadIdx.z*ldc] +
                                                    alpha*sum[0][threadIdx.y][threadIdx.z];
        }
    }
#endif
}

//==============================================================================

/**
    Purpose
    -------
    ZGEMM_REDUCE  performs one of the matrix-matrix operations
    
        C := alpha*A^T*B + beta*C,
    
    where alpha and beta are scalars, and A, B and C are matrices, with A
    a k-by-m matrix, B a k-by-n matrix, and C an m-by-n matrix.
    
    This routine is tuned for m, n << k. Typically, m and n are expected
    to be less than 128.

    @ingroup magma_tally2_zblas3
    ********************************************************************/
extern "C" void
magma_tally2blas_zgemm_reduce(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc )
{
    magma_tally2_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( k < 0 )
        info = -3;
    else if ( ldda < m )
        info = -6;
    else if ( lddb < k )
        info = -8;
    else if ( lddc < m )
        info = -11;
    
    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally2_int_t arch = magma_tally2_getdevice_arch();
    if ( arch < 200  ) {
        // --------------------
        // call CUDA ARCH 1.x -- maximum 512 threads
        const int NUM_THREADS = 512;
        const int BLK_K = (NUM_THREADS / (BLK_M * BLK_N)); // == 2
        dim3 blocks( (m-1)/BLK_M + 1, (n-1)/BLK_N + 1 );
        dim3 threads( BLK_K, BLK_M, BLK_N );
        zgemm_reduce_kernel<BLK_K> <<< blocks, threads, 0, magma_tally2_stream >>>
            ( m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc );
    }
    else {
        // --------------------
        // call CUDA ARCH 2.x -- maximum 1024 threads
        const int NUM_THREADS = 1024;
        const int BLK_K = (NUM_THREADS / (BLK_M * BLK_N)); // == 4
        dim3 blocks( (m-1)/BLK_M + 1, (n-1)/BLK_N + 1 );
        dim3 threads( BLK_K, BLK_M, BLK_N );
        zgemm_reduce_kernel<BLK_K> <<< blocks, threads, 0, magma_tally2_stream >>>
            ( m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc );
    }
}

//==============================================================================