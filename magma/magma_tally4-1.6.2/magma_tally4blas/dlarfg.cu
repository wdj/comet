/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlarfg.cu normal z -> d, Fri Jan 30 19:00:08 2015
       
       @author Mark Gates
*/
#include "common_magma_tally4.h"
#include "magma_tally4_templates.h"

#define REAL

// 512 is maximum number of threads for CUDA capability 1.x
#define NB 512


// ----------------------------------------
// CUDA kernel for magma_tally4_dlarfg.
// Uses one block of NB (currently 512) threads.
// Each thread sums dx[ tx + k*NB ]^2 for k = 0, 1, ...,
// then does parallel sum reduction to get norm-squared.
// 
// Currently setup to use NB threads, no matter how small dx is.
// This was slightly faster (5%) than passing n to magma_tally4_sum_reduce.
// To use number of threads = min( NB, max( 1, n-1 )), pass n as
// argument to magma_tally4_sum_reduce, rather than as template parameter.
__global__ void
dlarfg_kernel(
    int n,
    double* dalpha, double* dx, int incx,
    double* dtau )
{
    const int tx = threadIdx.x;
    __shared__ double swork[ NB ];
    // TODO is it faster for each thread to have its own scale (register)?
    // if so, communicate it via swork[0]
    __shared__ double sscale;
    __shared__ double sscale2;
    double tmp;
    
    // find max of [dalpha, dx], to use as scaling to avoid unnecesary under- and overflow
    if ( tx == 0 ) {
        tmp = *dalpha;
        #ifdef COMPLEX
        swork[tx] = max( fabs(real(tmp)), fabs(imag(tmp)) );
        #else
        swork[tx] = fabs(tmp);
        #endif
    }
    else {
        swork[tx] = 0;
    }
    for( int j = tx; j < n-1; j += NB ) {
        tmp = dx[j*incx];
        #ifdef COMPLEX
        swork[tx] = max( swork[tx], max( fabs(real(tmp)), fabs(imag(tmp)) ));
        #else
        swork[tx] = max( swork[tx], fabs(tmp) );
        #endif
    }
    magma_tally4_max_reduce< NB >( tx, swork );
    if ( tx == 0 )
        sscale = swork[0];
    __syncthreads();
    
    // sum norm^2 of dx/sscale
    // dx has length n-1
    swork[tx] = 0;
    if ( sscale > 0 ) {
        for( int j = tx; j < n-1; j += NB ) {
            tmp = dx[j*incx] / sscale;
            swork[tx] += real(tmp)*real(tmp) + imag(tmp)*imag(tmp);
        }
        magma_tally4_sum_reduce< NB >( tx, swork );
        //magma_tally4_sum_reduce( blockDim.x, tx, swork );
    }
    
    if ( tx == 0 ) {
        double alpha = *dalpha;
        if ( swork[0] == 0 && imag(alpha) == 0 ) {
            // H = I
            *dtau = MAGMA_tally4_D_ZERO;
        }
        else {
            // beta = norm( [dalpha, dx] )
            double beta;
            tmp  = alpha / sscale;
            beta = sscale * sqrt( real(tmp)*real(tmp) + imag(tmp)*imag(tmp) + swork[0] );
            beta = -copysign( beta, real(alpha) );
            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau   = MAGMA_tally4_D_MAKE( (beta - real(alpha)) / beta, -imag(alpha) / beta );
            *dalpha = MAGMA_tally4_D_MAKE( beta, 0 );
            sscale2 = 1 / (alpha - beta);
        }
    }
    
    // scale x (if norm was not 0)
    __syncthreads();
    if ( swork[0] != 0 ) {
        for( int j = tx; j < n-1; j += NB ) {
            dx[j*incx] *= sscale2;
        }
    }
}


/**
    Purpose
    -------
    DLARFG generates a real elementary reflector (Householder matrix)
    H of order n, such that

         H * ( alpha ) = ( beta ),   H**H * H = I.
             (   x   )   (   0  )

    where alpha and beta are scalars, with beta real and beta = ±norm([alpha, x]),
    and x is an (n-1)-element real vector. H is represented in the form

         H = I - tau * ( 1 ) * ( 1 v**H ),
                       ( v )

    where tau is a real scalar and v is a real (n-1)-element vector.
    Note that H is not symmetric.

    If the elements of x are all zero and dalpha is real, then tau = 0
    and H is taken to be the unit matrix.

    Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the elementary reflector.

    @param[in,out]
    dalpha  DOUBLE_PRECISION* on the GPU.
            On entry, pointer to the value alpha, i.e., the first entry of the vector.
            On exit, it is overwritten with the value beta.

    @param[in,out]
    dx      DOUBLE_PRECISION array, dimension (1+(N-2)*abs(INCX)), on the GPU
            On entry, the (n-1)-element vector x.
            On exit, it is overwritten with the vector v.

    @param[in]
    incx    INTEGER
            The increment between elements of X. INCX > 0.

    @param[out]
    dtau    DOUBLE_PRECISION* on the GPU.
            Pointer to the value tau.

    @param[in]
    queue   magma_tally4_queue_t
            Queue to execute in.

    @ingroup magma_tally4_daux1
    ********************************************************************/
extern "C"
void magma_tally4blas_dlarfg_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dalpha,
    magma_tally4Double_ptr dx, magma_tally4_int_t incx,
    magma_tally4Double_ptr dtau,
    magma_tally4_queue_t queue )
{
    dim3 blocks( 1 );
    dim3 threads( NB );
    //dim3 threads( min( NB, max( n-1, 1 )));
    dlarfg_kernel<<< blocks, threads, 0, queue >>>( n, dalpha, dx, incx, dtau );
}


/**
    @see magma_tally4blas_dlarfg_q
    @ingroup magma_tally4_daux1
    ********************************************************************/
extern "C"
void magma_tally4blas_dlarfg(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dalpha,
    magma_tally4Double_ptr dx, magma_tally4_int_t incx,
    magma_tally4Double_ptr dtau )
{
    dim3 blocks( 1 );
    dim3 threads( NB );
    //dim3 threads( min( NB, max( n-1, 1 )));
    dlarfg_kernel<<< blocks, threads >>>( n, dalpha, dx, incx, dtau );
}