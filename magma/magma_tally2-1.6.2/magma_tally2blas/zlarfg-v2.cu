/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define PRECISION_z


__global__
void magma_tally2_zlarfg_gpu_kernel( int n, magma_tally2DoubleComplex* dx0, magma_tally2DoubleComplex* dx,
                              magma_tally2DoubleComplex *dtau, double *dxnorm, magma_tally2DoubleComplex* dAkk)
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    __shared__ magma_tally2DoubleComplex scale;
    double xnorm;

    magma_tally2DoubleComplex dxi;

#if (defined(PRECISION_s) || defined(PRECISION_d))
    if( n <= 1 ) {
#else
    if( n <= 0 ) {
#endif
        *dtau = MAGMA_tally2_Z_ZERO;
        *dAkk = *dx0;
        return;
    }

    if ( j < n-1)
        dxi = dx[j];

    xnorm = *dxnorm;
    magma_tally2DoubleComplex alpha = *dx0;

#if (defined(PRECISION_s) || defined(PRECISION_d))
    if ( xnorm != 0 ) {
        if (i == 0) {  
            double beta  = sqrt( alpha*alpha + xnorm*xnorm );
            beta  = -copysign( beta, alpha );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau = (beta - alpha) / beta;
            *dAkk  = beta;

            scale = 1. / (alpha - beta);
        }
#else
    double alphar = MAGMA_tally2_Z_REAL(alpha);
    double alphai = MAGMA_tally2_Z_IMAG(alpha);
    if ( xnorm != 0 || alphai != 0) {
        if (i == 0) {
            double beta  = sqrt( alphar*alphar + alphai*alphai + xnorm*xnorm );
            beta  = -copysign( beta, alphar );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau = MAGMA_tally2_Z_MAKE((beta - alphar)/beta, -alphai/beta);
            *dAkk = MAGMA_tally2_Z_MAKE(beta, 0.);

            alpha = MAGMA_tally2_Z_MAKE( MAGMA_tally2_Z_REAL(alpha) - beta, MAGMA_tally2_Z_IMAG(alpha));
            scale = MAGMA_tally2_Z_DIV( MAGMA_tally2_Z_ONE, alpha);
        }
#endif

        // scale x
        __syncthreads();
        if ( xnorm != 0 && j < n-1)
            dx[j] = MAGMA_tally2_Z_MUL(dxi, scale);

    } else {
        *dtau = MAGMA_tally2_Z_ZERO;
        *dAkk = *dx0; 
    }
}


/*
    Generates Householder elementary reflector H = I - tau v v^T to reduce
        H [ dx0 ] = [ beta ]
          [ dx  ]   [ 0    ]
    with beta = ±norm( [dx0, dx] ) = ±dxnorm[0].
    Stores v over dx; first element of v is 1 and is not stored.
    Stores beta over dx0.
    Stores tau.  
    
    The difference with LAPACK's zlarfg is that the norm of dx, and hence beta,
    are computed outside the routine and passed to it in dxnorm (array on the GPU).
*/
extern "C" void
magma_tally2_zlarfg_gpu(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx0,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2Double_ptr        dxnorm,
    magma_tally2DoubleComplex_ptr dAkk)
{
    dim3 blocks((n+BLOCK_SIZE-1) / BLOCK_SIZE);
    dim3 threads( BLOCK_SIZE );

    /* recomputing the norm */
    //magma_tally2blas_dznrm2_cols(n, 1, dx0, n, dxnorm);
    magma_tally2blas_dznrm2_cols(n-1, 1, dx0+1, n, dxnorm);

    magma_tally2_zlarfg_gpu_kernel<<< blocks, threads,
                               0, magma_tally2_stream >>>(n, dx0, dx, dtau, dxnorm, dAkk);
}
