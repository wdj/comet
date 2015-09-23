/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_minproduct.h"

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
    magma_minproductDoubleComplex *A1;
    magma_minproductDoubleComplex *A2;
    int n, lda1, lda2;
} magma_minproductgpu_zlacpy_cnjg_params_t;

__global__ void magma_minproductgpu_zlacpy_cnjg( magma_minproductgpu_zlacpy_cnjg_params_t params )
{
    unsigned int x = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned int offset1 = x*params.lda1;
    unsigned int offset2 = x*params.lda2;
    if( x < params.n )
    {
        magma_minproductDoubleComplex *A1  = params.A1 + offset1;
        magma_minproductDoubleComplex *A2  = params.A2 + offset2;
        *A2 = MAGMA_minproduct_Z_CNJG(*A1);
    }
}


extern "C" void 
magma_minproductblas_zlacpy_cnjg_q(
    magma_minproduct_int_t n, magma_minproductDoubleComplex *dA1, magma_minproduct_int_t lda1, 
    magma_minproductDoubleComplex *dA2, magma_minproduct_int_t lda2,
    magma_minproduct_queue_t queue )
{
    int blocksize = 64;
    dim3 blocks( (n+blocksize-1) / blocksize, 1, 1);
    magma_minproductgpu_zlacpy_cnjg_params_t params = { dA1, dA2, n, lda1, lda2 };
    magma_minproductgpu_zlacpy_cnjg<<< blocks, blocksize, 0, queue >>>( params );
}


extern "C" void 
magma_minproductblas_zlacpy_cnjg(
    magma_minproduct_int_t n, magma_minproductDoubleComplex *dA1, magma_minproduct_int_t lda1, 
    magma_minproductDoubleComplex *dA2, magma_minproduct_int_t lda2)
{
    magma_minproductblas_zlacpy_cnjg_q( n, dA1, lda1, dA2, lda2, magma_minproduct_stream );
}
