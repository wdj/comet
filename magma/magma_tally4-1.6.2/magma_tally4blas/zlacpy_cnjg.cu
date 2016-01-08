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

/*********************************************************
 *
 * SWAP BLAS: permute to set of N elements
 *
 ********************************************************/
/*
 *  First version: line per line
 */
typedef struct {
    magma_tally4DoubleComplex *A1;
    magma_tally4DoubleComplex *A2;
    int n, lda1, lda2;
} magma_tally4gpu_zlacpy_cnjg_params_t;

__global__ void magma_tally4gpu_zlacpy_cnjg( magma_tally4gpu_zlacpy_cnjg_params_t params )
{
    unsigned int x = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned int offset1 = x*params.lda1;
    unsigned int offset2 = x*params.lda2;
    if( x < params.n )
    {
        magma_tally4DoubleComplex *A1  = params.A1 + offset1;
        magma_tally4DoubleComplex *A2  = params.A2 + offset2;
        *A2 = MAGMA_tally4_Z_CNJG(*A1);
    }
}


extern "C" void 
magma_tally4blas_zlacpy_cnjg_q(
    magma_tally4_int_t n, magma_tally4DoubleComplex *dA1, magma_tally4_int_t lda1, 
    magma_tally4DoubleComplex *dA2, magma_tally4_int_t lda2,
    magma_tally4_queue_t queue )
{
    int blocksize = 64;
    dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
    magma_tally4gpu_zlacpy_cnjg_params_t params = { dA1, dA2, n, lda1, lda2 };
    magma_tally4gpu_zlacpy_cnjg<<< blocks, blocksize, 0, queue >>>( params );
}


extern "C" void 
magma_tally4blas_zlacpy_cnjg(
    magma_tally4_int_t n, magma_tally4DoubleComplex *dA1, magma_tally4_int_t lda1, 
    magma_tally4DoubleComplex *dA2, magma_tally4_int_t lda2)
{
    magma_tally4blas_zlacpy_cnjg_q( n, dA1, lda1, dA2, lda2, magma_tally4_stream );
}
