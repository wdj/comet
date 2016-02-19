/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlarfbx.cu normal z -> s, Fri Jan 30 19:00:08 2015

*/
#include "common_magma_tally3.h"
#include "commonblas_s.h"
#include "magma_tally3_templates.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512




//==============================================================================
extern "C"
__global__ void 
magma_tally3_sgemv_kernel1(int m, const float * __restrict__ V, int ldv, 
                    const float * __restrict__ c, 
                    float *dwork)
{
    const int i = threadIdx.x;
    const float *dV = V + (blockIdx.x) * ldv;

    __shared__ float sum[ BLOCK_SIZE ];
    float lsum;

    /*  lsum := v**H * C  */
    lsum = MAGMA_tally3_S_ZERO;
    for( int j = i; j < m; j += BLOCK_SIZE )
       lsum += MAGMA_tally3_S_MUL( MAGMA_tally3_S_CNJG( dV[j] ), c[j] );
    
    sum[i] = lsum;
    magma_tally3_sum_reduce< BLOCK_SIZE >( i, sum );

    __syncthreads();
    if (i==0)
       dwork [blockIdx.x] = sum[0];
}

//==============================================================================
/*  ----------------------------------------------------------------------------- 
    Call 
        magma_tally3_sgemv_kernel3<<< n, BLOCK_SIZE>>>(m, V, ldv, c, dwork, tau)
    to compute
        SGEMV( "Conjugate transpose", m, n, -tau[0], V, ldv, c, 1, zero, dwork, 1)
        and to set c[0] to 1.
    i.e., 
        work = -tau[0] V**H c
    ----------------------------------------------------------------------------- */
extern "C"
__global__ void
magma_tally3_sgemv_kernel3(int m, const float * __restrict__ V, int ldv, float *c,
                    float *dwork, float *tau)
{
    const int i = threadIdx.x;
    const float *dV = V + (blockIdx.x) * ldv;

    __shared__ float sum[ BLOCK_SIZE ];
    float lsum;

    if (i==0)
       c[0] = MAGMA_tally3_S_ONE;           

    /*  lsum := v**H * C  */
    lsum = MAGMA_tally3_S_ZERO;
    for( int j = i; j < m; j += BLOCK_SIZE )
       lsum += MAGMA_tally3_S_MUL( MAGMA_tally3_S_CNJG( dV[j] ), c[j] );

    sum[i] = lsum;
    magma_tally3_sum_reduce< BLOCK_SIZE >( i, sum );

    __syncthreads();
    if (i==0)
       dwork [blockIdx.x] = -tau[0]*sum[0];
}

//==============================================================================
extern "C"
__global__ void
magma_tally3_sgemv_kernel2(int m, int n, const float * __restrict__ V, int ldv, 
                    const float * __restrict__ x, float *c)
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    float lsum;

    V += j;

    lsum = MAGMA_tally3_S_ZERO;
    if (j < m){
       for(int k=0; k<n; k++)
          lsum += MAGMA_tally3_S_MUL( V[k*ldv], x[k]);
       
       c[j] -= lsum;
    }
}

//==============================================================================

/*
    Apply a real block reflector H to a real vector C from the left
    (i.e., C = H C). H is represented in the form
          H = I - V T V**H
    where T is the real k-by-k upper triangular matrix in the 
    representation of the block reflector, and V is a real block of
    k elementary reflectors. 
*/
extern "C" void
magma_tally3_slarfbx_gpu(
    magma_tally3_int_t m, magma_tally3_int_t k,
    magma_tally3Float_ptr V,  magma_tally3_int_t ldv,
    magma_tally3Float_ptr dT, magma_tally3_int_t ldt,
    magma_tally3Float_ptr c,
    magma_tally3Float_ptr dwork)
{
    /* dwork = V**H c     */
    magma_tally3_sgemv_kernel1<<< k, BLOCK_SIZE, 0, magma_tally3_stream >>>(m, V, ldv, c, dwork); 

    /* dwork = T**H dwork */
    magma_tally3_strmv_tkernel<<< k, k, 0, magma_tally3_stream >>>( dT, ldt, dwork, dwork+k);
 
    /* c = c - V dwork    */
    dim3  blocks3( (m + BLOCK_SIZE-1) / BLOCK_SIZE );
    dim3 threads3( BLOCK_SIZE );     
    magma_tally3_sgemv_kernel2<<< blocks3, threads3, 0, magma_tally3_stream >>>( m, k, V, ldv, dwork+k, c);
}

//==============================================================================
