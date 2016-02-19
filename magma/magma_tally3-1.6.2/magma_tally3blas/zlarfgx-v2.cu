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

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define PRECISION_z


//==============================================================================

__global__
void magma_tally3_zlarfgx_gpu_kernel( int n, magma_tally3DoubleComplex* dx0, magma_tally3DoubleComplex* dx,
                               magma_tally3DoubleComplex *dtau, double *dxnorm,
                               magma_tally3DoubleComplex *dA, int it)
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    __shared__ magma_tally3DoubleComplex scale;
    __shared__ double xnorm;
  
    magma_tally3DoubleComplex dxi;

    if ( j < n-1 )
        dxi = dx[j];
  
    if ( i == 0 ) {
        xnorm = *dxnorm;
#if (defined(PRECISION_s) || defined(PRECISION_d))
        double alpha = *dx0;
        double alphai = MAGMA_tally3_Z_ZERO;
        if ( (xnorm == 0 && alphai == MAGMA_tally3_Z_ZERO ) || n == 1 )
#else
        magma_tally3DoubleComplex alpha = *dx0;
        double alphar =  MAGMA_tally3_Z_REAL(alpha), alphai = MAGMA_tally3_Z_IMAG(alpha);
        if ( (xnorm == 0 && alphai == MAGMA_tally3_Z_ZERO ) || n == 0 )
#endif
        {
            *dtau = MAGMA_tally3_Z_ZERO;
            *dA   = *dx0;
        }
        else {

#if (defined(PRECISION_s) || defined(PRECISION_d))
            // no need to compute the norm as it is passed as input
            double beta  = xnorm; // sqrt( alpha*alpha + xnorm*xnorm );
            beta  = -copysign( beta, alpha );
 
            // todo: deal with badly scaled vectors (see lapack's larfg)
            if (j==0){
                *dtau = (beta - alpha) / beta;
                //*dx0  = 1.; //cannot be done here because raise condition all threadblock need to read it for alpha
                *dA   = beta;
            }

            scale = 1. / (alpha - beta);
#else
            // no need to compute the norm as it is passed as input
            double beta  = xnorm; // sqrt( alphar*alphar + alphai*alphai + xnorm*xnorm );
            beta  = -copysign( beta, alphar );

            // todo: deal with badly scaled vectors (see lapack's larfg)
            if (j==0){
                *dtau = MAGMA_tally3_Z_MAKE((beta - alphar)/beta, -alphai/beta);
                //*dx0  = MAGMA_tally3_Z_MAKE(  1., 0.); //cannot be done here because raise condition all threadblock need to read it for alpha
                *dA   = MAGMA_tally3_Z_MAKE(beta, 0.);
            }

            alpha = MAGMA_tally3_Z_MAKE( MAGMA_tally3_Z_REAL(alpha) - beta, MAGMA_tally3_Z_IMAG(alpha));
            scale = MAGMA_tally3_Z_DIV( MAGMA_tally3_Z_ONE, alpha);
#endif
        }
    }

    // scale x
    __syncthreads();
    if ( xnorm != 0 && j < n-1)
        dx[j] = MAGMA_tally3_Z_MUL(dxi, scale);

    if (j<it){
        *( dA-it+j) = *(dx0-it+j);
        *(dx0-it+j) = MAGMA_tally3_Z_MAKE(0., 0.);
    }
}

//==============================================================================

/*
    Generates Householder elementary reflector H = I - tau v v^T to reduce
        H [ dx0 ] = [ beta ]
          [ dx  ]   [ 0    ]
    with beta = ±norm( [dx0, dx] ) = ±dxnorm[0].
    Stores v over dx; first element of v is 1 and is not stored.
    Stores beta over dx0.
    Stores tau.
    
    The difference with LAPACK's zlarfg is that the norm of dx, and hance beta,
    are computed outside the routine and passed to it in dxnorm (array on the GPU).
*/
extern "C" void
magma_tally3_zlarfgx_gpu(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dx0,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3Double_ptr        dxnorm,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t iter)
{
    dim3 blocks((n+BLOCK_SIZE-1) / BLOCK_SIZE);
    dim3 threads( BLOCK_SIZE );
 
    magma_tally3_zlarfgx_gpu_kernel<<< blocks, threads, 0, magma_tally3_stream >>>( n, dx0, dx, dtau, dxnorm, dA, iter);
}

//==============================================================================

/*
    Generates Householder elementary reflector H = I - tau v v^T to reduce
        H [ dx0 ] = [ beta ]
          [ dx  ]   [ 0    ]
    with beta = ±norm( [dx0, dx] ) = ±dxnorm[0].
    Stores v over dx; first element of v is 1 and is not stored.
    Stores beta over dx0.
    Stores tau.
    
    The difference with LAPACK's zlarfg is that the norm of dx, and hance beta,
    are computed outside the routine and passed to it in dxnorm (array on the GPU).
*/
extern "C" void
magma_tally3_zlarfgtx_gpu(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dx0,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3Double_ptr        dxnorm,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t iter,
    magma_tally3DoubleComplex_ptr V,  magma_tally3_int_t ldv,
    magma_tally3DoubleComplex_ptr T,  magma_tally3_int_t ldt,
    magma_tally3DoubleComplex_ptr dwork)
{
    /*  Generate the elementary reflector H(iter)  */
    magma_tally3_zlarfgx_gpu(n, dx0, dx, dtau, dxnorm, dA, iter);
    
    if (iter==0) {
        magma_tally3DoubleComplex tt = MAGMA_tally3_Z_ONE;
        magma_tally3blas_zlacpy(Magma_tally3UpperLower, 1, 1, dtau, 1, T+iter+iter*ldt, 1);
        magma_tally3_zsetmatrix(1,1, &tt,1, dx0,1);
    }
    else {
        /* Compute the iter-th column of T */
        magma_tally3_zgemv_kernel3<<< iter, BLOCK_SIZE, 0, magma_tally3_stream >>>( n, V, ldv, dx0, dwork, dtau );
        magma_tally3_ztrmv_kernel2<<< iter, iter,       0, magma_tally3_stream >>>( T, ldt, dwork, T+iter*ldt, dtau );
    }
}

//==============================================================================
