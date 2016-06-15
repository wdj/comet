/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally3.h"
#include "commonblas_z.h"
#include "magma_tally3_templates.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512




//==============================================================================
extern "C"
__global__ void 
magma_tally3_zgemv_kernel1(int m, const magma_tally3DoubleComplex * __restrict__ V, int ldv, 
                    const magma_tally3DoubleComplex * __restrict__ c, 
                    magma_tally3DoubleComplex *dwork)
{
    const int i = threadIdx.x;
    const magma_tally3DoubleComplex *dV = V + (blockIdx.x) * ldv;

    __shared__ magma_tally3DoubleComplex sum[ BLOCK_SIZE ];
    magma_tally3DoubleComplex lsum;

    /*  lsum := v**H * C  */
    lsum = MAGMA_tally3_Z_ZERO;
    for( int j = i; j < m; j += BLOCK_SIZE )
       lsum += MAGMA_tally3_Z_MUL( MAGMA_tally3_Z_CNJG( dV[j] ), c[j] );
    
    sum[i] = lsum;
    magma_tally3_sum_reduce< BLOCK_SIZE >( i, sum );

    __syncthreads();
    if (i==0)
       dwork [blockIdx.x] = sum[0];
}

//==============================================================================
/*  ----------------------------------------------------------------------------- 
    Call 
        magma_tally3_zgemv_kernel3<<< n, BLOCK_SIZE>>>(m, V, ldv, c, dwork, tau)
    to compute
        ZGEMV( "Conjugate transpose", m, n, -tau[0], V, ldv, c, 1, zero, dwork, 1)
        and to set c[0] to 1.
    i.e., 
        work = -tau[0] V**H c
    ----------------------------------------------------------------------------- */
extern "C"
__global__ void
magma_tally3_zgemv_kernel3(int m, const magma_tally3DoubleComplex * __restrict__ V, int ldv, magma_tally3DoubleComplex *c,
                    magma_tally3DoubleComplex *dwork, magma_tally3DoubleComplex *tau)
{
    const int i = threadIdx.x;
    const magma_tally3DoubleComplex *dV = V + (blockIdx.x) * ldv;

    __shared__ magma_tally3DoubleComplex sum[ BLOCK_SIZE ];
    magma_tally3DoubleComplex lsum;

    if (i==0)
       c[0] = MAGMA_tally3_Z_ONE;           

    /*  lsum := v**H * C  */
    lsum = MAGMA_tally3_Z_ZERO;
    for( int j = i; j < m; j += BLOCK_SIZE )
       lsum += MAGMA_tally3_Z_MUL( MAGMA_tally3_Z_CNJG( dV[j] ), c[j] );

    sum[i] = lsum;
    magma_tally3_sum_reduce< BLOCK_SIZE >( i, sum );

    __syncthreads();
    if (i==0)
       dwork [blockIdx.x] = -tau[0]*sum[0];
}

//==============================================================================
extern "C"
__global__ void
magma_tally3_zgemv_kernel2(int m, int n, const magma_tally3DoubleComplex * __restrict__ V, int ldv, 
                    const magma_tally3DoubleComplex * __restrict__ x, magma_tally3DoubleComplex *c)
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    magma_tally3DoubleComplex lsum;

    V += j;

    lsum = MAGMA_tally3_Z_ZERO;
    if (j < m){
       for(int k=0; k<n; k++)
          lsum += MAGMA_tally3_Z_MUL( V[k*ldv], x[k]);
       
       c[j] -= lsum;
    }
}

//==============================================================================

/*
    Apply a complex block reflector H to a complex vector C from the left
    (i.e., C = H C). H is represented in the form
          H = I - V T V**H
    where T is the complex k-by-k upper triangular matrix in the 
    representation of the block reflector, and V is a complex block of
    k elementary reflectors. 
*/
extern "C" void
magma_tally3_zlarfbx_gpu(
    magma_tally3_int_t m, magma_tally3_int_t k,
    magma_tally3DoubleComplex_ptr V,  magma_tally3_int_t ldv,
    magma_tally3DoubleComplex_ptr dT, magma_tally3_int_t ldt,
    magma_tally3DoubleComplex_ptr c,
    magma_tally3DoubleComplex_ptr dwork)
{
    /* dwork = V**H c     */
    magma_tally3_zgemv_kernel1<<< k, BLOCK_SIZE, 0, magma_tally3_stream >>>(m, V, ldv, c, dwork); 

    /* dwork = T**H dwork */
    magma_tally3_ztrmv_tkernel<<< k, k, 0, magma_tally3_stream >>>( dT, ldt, dwork, dwork+k);
 
    /* c = c - V dwork    */
    dim3  blocks3( (m + BLOCK_SIZE-1) / BLOCK_SIZE );
    dim3 threads3( BLOCK_SIZE );     
    magma_tally3_zgemv_kernel2<<< blocks3, threads3, 0, magma_tally3_stream >>>( m, k, V, ldv, dwork+k, c);
}

//==============================================================================