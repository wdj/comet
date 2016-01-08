/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

#define BLOCK_SIZE 64

/*********************************************************/
/*
 *  Blocked version: swap several pairs of lines
 */
typedef struct {
    magma_tally4DoubleComplex *A;
    magma_tally4DoubleComplex *B;
    int n, ldda, lddb, npivots;
    short ipiv[BLOCK_SIZE];
} magma_tally4gpu_zswapblk_params_t;

__global__ void magma_tally4gpu_zswapblkrm( magma_tally4gpu_zswapblk_params_t params )
{
    unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
    if( y < params.n )
    {
        magma_tally4DoubleComplex *A = params.A + y - params.ldda;
        magma_tally4DoubleComplex *B = params.B + y;
      
        for( int i = 0; i < params.npivots; i++ )
        {
            A += params.ldda;
            if ( params.ipiv[i] == -1 )
                continue;
            magma_tally4DoubleComplex  tmp1 = *A;
            magma_tally4DoubleComplex *tmp2 = B + params.ipiv[i]*params.lddb;
            *A    = *tmp2;
            *tmp2 =  tmp1;
        }
    }
}

__global__ void magma_tally4gpu_zswapblkcm( magma_tally4gpu_zswapblk_params_t params )
{
    unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned int offset1 = y*params.ldda;
    unsigned int offset2 = y*params.lddb;
    if( y < params.n )
    {
        magma_tally4DoubleComplex *A = params.A + offset1 - 1;
        magma_tally4DoubleComplex *B = params.B + offset2;
      
        for( int i = 0; i < params.npivots; i++ )
        {
            A++;
            if ( params.ipiv[i] == -1 )
                continue;
            magma_tally4DoubleComplex  tmp1 = *A;
            magma_tally4DoubleComplex *tmp2 = B + params.ipiv[i];
            *A    = *tmp2;
            *tmp2 =  tmp1;
        }
    }
    __syncthreads();
}


/**
    @ingroup magma_tally4_zblas2
    ********************************************************************/
extern "C" void 
magma_tally4blas_zswapblk_q(
    magma_tally4_order_t order, magma_tally4_int_t n, 
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t i1, magma_tally4_int_t i2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci, magma_tally4_int_t offset,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t  blocksize = 64;
    dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
    magma_tally4_int_t  k, im;
    
    /* Quick return */
    if ( n == 0 )
        return;
    
    if ( order == Magma_tally4ColMajor ) {
        for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
        {
            magma_tally4_int_t sb = min(BLOCK_SIZE, i2-k);
            magma_tally4gpu_zswapblk_params_t params = { dA+k, dB, n, ldda, lddb, sb };
            for( magma_tally4_int_t j = 0; j < sb; j++ )
            {
                im = ipiv[(k+j)*inci] - 1;
                if ( (k+j) == im )
                    params.ipiv[j] = -1;
                else
                    params.ipiv[j] = im - offset;
            }
            magma_tally4gpu_zswapblkcm<<< blocks, blocksize, 0, queue >>>( params );
        }
    }
    else {
        for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
        {
            magma_tally4_int_t sb = min(BLOCK_SIZE, i2-k);
            magma_tally4gpu_zswapblk_params_t params = { dA+k*ldda, dB, n, ldda, lddb, sb };
            for( magma_tally4_int_t j = 0; j < sb; j++ )
            {
                im = ipiv[(k+j)*inci] - 1;
                if ( (k+j) == im )
                    params.ipiv[j] = -1;
                else
                    params.ipiv[j] = im - offset;
            }
            magma_tally4gpu_zswapblkrm<<< blocks, blocksize, 0, queue >>>( params );
        }
    }
}


/**
    @see magma_tally4blas_zswapblk_q
    @ingroup magma_tally4_zblas2
    ********************************************************************/
extern "C" void 
magma_tally4blas_zswapblk(
    magma_tally4_order_t order, magma_tally4_int_t n, 
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t i1, magma_tally4_int_t i2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci, magma_tally4_int_t offset )
{
    magma_tally4blas_zswapblk_q(
        order, n, dA, ldda, dB, lddb, i1, i2, ipiv, inci, offset, magma_tally4_stream );
}
