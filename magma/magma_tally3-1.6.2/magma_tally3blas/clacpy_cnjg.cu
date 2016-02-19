/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlacpy_cnjg.cu normal z -> c, Fri Jan 30 19:00:08 2015

*/
#include "common_magma_tally3.h"

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
    magma_tally3FloatComplex *A1;
    magma_tally3FloatComplex *A2;
    int n, lda1, lda2;
} magma_tally3gpu_clacpy_cnjg_params_t;

__global__ void magma_tally3gpu_clacpy_cnjg( magma_tally3gpu_clacpy_cnjg_params_t params )
{
    unsigned int x = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned int offset1 = x*params.lda1;
    unsigned int offset2 = x*params.lda2;
    if( x < params.n )
    {
        magma_tally3FloatComplex *A1  = params.A1 + offset1;
        magma_tally3FloatComplex *A2  = params.A2 + offset2;
        *A2 = MAGMA_tally3_C_CNJG(*A1);
    }
}


extern "C" void 
magma_tally3blas_clacpy_cnjg_q(
    magma_tally3_int_t n, magma_tally3FloatComplex *dA1, magma_tally3_int_t lda1, 
    magma_tally3FloatComplex *dA2, magma_tally3_int_t lda2,
    magma_tally3_queue_t queue )
{
    int blocksize = 64;
    dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
    magma_tally3gpu_clacpy_cnjg_params_t params = { dA1, dA2, n, lda1, lda2 };
    magma_tally3gpu_clacpy_cnjg<<< blocks, blocksize, 0, queue >>>( params );
}


extern "C" void 
magma_tally3blas_clacpy_cnjg(
    magma_tally3_int_t n, magma_tally3FloatComplex *dA1, magma_tally3_int_t lda1, 
    magma_tally3FloatComplex *dA2, magma_tally3_int_t lda2)
{
    magma_tally3blas_clacpy_cnjg_q( n, dA1, lda1, dA2, lda2, magma_tally3_stream );
}
