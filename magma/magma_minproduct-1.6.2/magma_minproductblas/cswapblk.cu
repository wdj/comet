/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zswapblk.cu normal z -> c, Fri Jan 30 19:00:09 2015

*/
#include "common_magma_minproduct.h"

#define BLOCK_SIZE 64

/*********************************************************/
/*
 *  Blocked version: swap several pairs of lines
 */
typedef struct {
    magma_minproductFloatComplex *A;
    magma_minproductFloatComplex *B;
    int n, ldda, lddb, npivots;
    short ipiv[BLOCK_SIZE];
} magma_minproductgpu_cswapblk_params_t;

__global__ void magma_minproductgpu_cswapblkrm( magma_minproductgpu_cswapblk_params_t params )
{
    unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
    if( y < params.n )
    {
        magma_minproductFloatComplex *A = params.A + y - params.ldda;
        magma_minproductFloatComplex *B = params.B + y;
      
        for( int i = 0; i < params.npivots; i++ )
        {
            A += params.ldda;
            if ( params.ipiv[i] == -1 )
                continue;
            magma_minproductFloatComplex  tmp1 = *A;
            magma_minproductFloatComplex *tmp2 = B + params.ipiv[i]*params.lddb;
            *A    = *tmp2;
            *tmp2 =  tmp1;
        }
    }
}

__global__ void magma_minproductgpu_cswapblkcm( magma_minproductgpu_cswapblk_params_t params )
{
    unsigned int y = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned int offset1 = y*params.ldda;
    unsigned int offset2 = y*params.lddb;
    if( y < params.n )
    {
        magma_minproductFloatComplex *A = params.A + offset1 - 1;
        magma_minproductFloatComplex *B = params.B + offset2;
      
        for( int i = 0; i < params.npivots; i++ )
        {
            A++;
            if ( params.ipiv[i] == -1 )
                continue;
            magma_minproductFloatComplex  tmp1 = *A;
            magma_minproductFloatComplex *tmp2 = B + params.ipiv[i];
            *A    = *tmp2;
            *tmp2 =  tmp1;
        }
    }
    __syncthreads();
}


/**
    @ingroup magma_minproduct_cblas2
    ********************************************************************/
extern "C" void 
magma_minproductblas_cswapblk_q(
    magma_minproduct_order_t order, magma_minproduct_int_t n, 
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci, magma_minproduct_int_t offset,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t  blocksize = 64;
    dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
    magma_minproduct_int_t  k, im;
    
    /* Quick return */
    if ( n == 0 )
        return;
    
    if ( order == Magma_minproductColMajor ) {
        for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
        {
            magma_minproduct_int_t sb = min(BLOCK_SIZE, i2-k);
            magma_minproductgpu_cswapblk_params_t params = { dA+k, dB, n, ldda, lddb, sb };
            for( magma_minproduct_int_t j = 0; j < sb; j++ )
            {
                im = ipiv[(k+j)*inci] - 1;
                if ( (k+j) == im )
                    params.ipiv[j] = -1;
                else
                    params.ipiv[j] = im - offset;
            }
            magma_minproductgpu_cswapblkcm<<< blocks, blocksize, 0, queue >>>( params );
        }
    }
    else {
        for( k=(i1-1); k<i2; k+=BLOCK_SIZE )
        {
            magma_minproduct_int_t sb = min(BLOCK_SIZE, i2-k);
            magma_minproductgpu_cswapblk_params_t params = { dA+k*ldda, dB, n, ldda, lddb, sb };
            for( magma_minproduct_int_t j = 0; j < sb; j++ )
            {
                im = ipiv[(k+j)*inci] - 1;
                if ( (k+j) == im )
                    params.ipiv[j] = -1;
                else
                    params.ipiv[j] = im - offset;
            }
            magma_minproductgpu_cswapblkrm<<< blocks, blocksize, 0, queue >>>( params );
        }
    }
}


/**
    @see magma_minproductblas_cswapblk_q
    @ingroup magma_minproduct_cblas2
    ********************************************************************/
extern "C" void 
magma_minproductblas_cswapblk(
    magma_minproduct_order_t order, magma_minproduct_int_t n, 
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci, magma_minproduct_int_t offset )
{
    magma_minproductblas_cswapblk_q(
        order, n, dA, ldda, dB, lddb, i1, i2, ipiv, inci, offset, magma_minproduct_stream );
}
