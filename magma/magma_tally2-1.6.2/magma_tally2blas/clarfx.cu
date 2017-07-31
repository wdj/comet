/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlarfx.cu normal z -> c, Fri Jan 30 19:00:09 2015

*/
#include "common_magma_tally2.h"
#include "commonblas_c.h"
#include "magma_tally2_templates.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define BLOCK_SIZEx  32
#define BLOCK_SIZEy  16


//==============================================================================

__global__
void magma_tally2_clarfx_kernel( int m, magma_tally2FloatComplex *v, magma_tally2FloatComplex *tau,
                         magma_tally2FloatComplex *c, int ldc, float *xnorm,
                         magma_tally2FloatComplex *T, int it )
{
    if ( !MAGMA_tally2_C_EQUAL(*tau, MAGMA_tally2_C_ZERO) ) {
        const int tx = threadIdx.x;
        //magma_tally2FloatComplex *dc = c + (blockIdx.x-it-1) * ldc;
        magma_tally2FloatComplex *dc = c + (blockIdx.x) * ldc;

        __shared__ magma_tally2FloatComplex sum[ BLOCK_SIZE ];
        magma_tally2FloatComplex lsum;

        /* NOTE HERE C is the C at position C(i, 0) 
         * if blockIdx.x<it it performs the V(i:n,i)' * V(i:n,1:i-1)' used for computing T
         * if blockIdx.x>it it perform  w := v**H * C  */
        lsum = MAGMA_tally2_C_ZERO;
        for( int j = tx; j < m; j += BLOCK_SIZE ){
            if (j==0){
               lsum += MAGMA_tally2_C_MUL( MAGMA_tally2_C_ONE, dc[j] );
               v[j] = MAGMA_tally2_C_ONE;
            }
            else
               lsum += MAGMA_tally2_C_MUL( MAGMA_tally2_C_CNJG( v[j] ), dc[j] );
        }
        sum[tx] = lsum;
        magma_tally2_sum_reduce< BLOCK_SIZE >( tx, sum );

        /*  C := C - v * w  */
        __syncthreads();
        magma_tally2FloatComplex z__1 = - MAGMA_tally2_C_CNJG(*tau) * sum[0];
        if (blockIdx.x>it){
           for( int j = m-tx-1; j>=0 ; j -= BLOCK_SIZE )
                 dc[j] += z__1 * v[j];
           __syncthreads();

           /* Adjust the rest of the column norms */
           /*
           if (tx==0){
             float temp = MAGMA_tally2_C_ABS( dc[0] ) / xnorm[blockIdx.x-it-1];
             temp = (temp + 1.) * (1. - temp);
             xnorm[blockIdx.x-it-1] = xnorm[blockIdx.x-it-1] * sqrt(temp); 
           }
           */
        }
        else
        {
           if (blockIdx.x==it)
              *(T+it) = *tau;
           else
              *(T+blockIdx.x) = MAGMA_tally2_C_CNJG(z__1);
        }
    }
    else if (blockIdx.x<=it)// in case tau is zero put the corresponding column of T to zero
    {
        *(T+blockIdx.x) = MAGMA_tally2_C_ZERO;
    }

}

//==============================================================================
extern "C"
__global__
void magma_tally2_ctrmv_kernel(const magma_tally2FloatComplex *T, int ldt, magma_tally2FloatComplex *t)
{
   const int tx = threadIdx.x;
   T += tx;

   __shared__ magma_tally2FloatComplex tlocal[ BLOCK_SIZE ];
   magma_tally2FloatComplex res = MAGMA_tally2_C_MAKE(0., 0.);

   tlocal[tx] = t[tx];
   __syncthreads();

   #pragma unroll
   for(int j=0; j<blockDim.x; j++)
      res +=  T[j*ldt]*tlocal[j];

   t[tx] = res;
}

extern "C"
__global__
void magma_tally2_ctrmv_kernel2(const magma_tally2FloatComplex *T, int ldt, magma_tally2FloatComplex *t, 
                         magma_tally2FloatComplex *y, magma_tally2FloatComplex *tau)
{
   const int tx = threadIdx.x;
   T += blockIdx.x;

   __shared__ magma_tally2FloatComplex sum[ 128 ];

   sum[tx] = T[tx*ldt]*t[tx];
   magma_tally2_sum_reduce_n(blockDim.x, tx, sum);

   __syncthreads();

   if (tx==0){
      y[blockIdx.x] = sum[0];
      if (blockIdx.x==0)
         y[gridDim.x] = tau[0];
   }
}

//==============================================================================
extern "C"
__global__
void magma_tally2_ctrmv_tkernel(magma_tally2FloatComplex *T, int ldt, magma_tally2FloatComplex *t, magma_tally2FloatComplex *y)
{
   const int tx = threadIdx.x;
   T += blockIdx.x*ldt;

   __shared__ magma_tally2FloatComplex sum[ 128 ];

   sum[tx] = MAGMA_tally2_C_CNJG(T[tx])*t[tx];
   magma_tally2_sum_reduce_n(blockDim.x, tx, sum);

   __syncthreads();

   if (tx==0)
      y[blockIdx.x] = sum[0];
}

//==============================================================================

/*
    Apply a complex elementary reflector H to a complex M-by-N
    matrix C from the left. H is represented in the form
          H = I - tau * v * v**H
    where tau is a complex scalar and v is a complex vector.
    If tau = 0, then H is taken to be the unit matrix.

    To apply H**H (the conjugate transpose of H), supply conjg(tau) 
    instead tau.

    The norms of v(:, 1:n) are given as input in xnorm(1:n). On exit, the norms
    are adjusted to hold the norms of v(2:m,2:n). This is a difference with the 
    LAPACK's clarf routine. 
 */
extern "C" void
magma_tally2_clarfx_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr v,
    magma_tally2FloatComplex_ptr tau,
    magma_tally2FloatComplex_ptr C, magma_tally2_int_t ldc,
    magma_tally2Float_ptr        xnorm, 
    magma_tally2FloatComplex_ptr dT, magma_tally2_int_t iter,
    magma_tally2FloatComplex_ptr work )
{
    magma_tally2_int_t N = n + iter + 1;

    if (iter==0)
        magma_tally2_clarfx_kernel<<< N, BLOCK_SIZE, 0, magma_tally2_stream >>>( m, v, tau, C, ldc, xnorm, dT+iter*N, iter);
    else
        magma_tally2_clarfx_kernel<<< N, BLOCK_SIZE, 0, magma_tally2_stream >>>( m, v, tau, C, ldc, xnorm, work, iter);

    if (iter > 0){
        //magma_tally2_ctrmv_kernel<<< 1, iter, 0, magma_tally2_stream >>>( dT, N, dT+iter*N);
        magma_tally2_ctrmv_kernel2<<< iter, iter, 0, magma_tally2_stream  >>>( dT, N, work, dT+iter*N, tau);
    }
}

//==============================================================================
