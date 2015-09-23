/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlarfgx-v2.cu normal z -> d, Fri Jan 30 19:00:09 2015

*/
#include "common_magma_minproduct.h"
#include "commonblas_d.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define PRECISION_d


//==============================================================================

__global__
void magma_minproduct_dlarfgx_gpu_kernel( int n, double* dx0, double* dx,
                               double *dtau, double *dxnorm,
                               double *dA, int it)
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    __shared__ double scale;
    __shared__ double xnorm;
  
    double dxi;

    if ( j < n-1 )
        dxi = dx[j];
  
    if ( i == 0 ) {
        xnorm = *dxnorm;
#if (defined(PRECISION_s) || defined(PRECISION_d))
        double alpha = *dx0;
        double alphai = MAGMA_minproduct_D_ZERO;
        if ( (xnorm == 0 && alphai == MAGMA_minproduct_D_ZERO ) || n == 1 )
#else
        double alpha = *dx0;
        double alphar =  MAGMA_minproduct_D_REAL(alpha), alphai = MAGMA_minproduct_D_IMAG(alpha);
        if ( (xnorm == 0 && alphai == MAGMA_minproduct_D_ZERO ) || n == 0 )
#endif
        {
            *dtau = MAGMA_minproduct_D_ZERO;
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
                *dtau = MAGMA_minproduct_D_MAKE((beta - alphar)/beta, -alphai/beta);
                //*dx0  = MAGMA_minproduct_D_MAKE(  1., 0.); //cannot be done here because raise condition all threadblock need to read it for alpha
                *dA   = MAGMA_minproduct_D_MAKE(beta, 0.);
            }

            alpha = MAGMA_minproduct_D_MAKE( MAGMA_minproduct_D_REAL(alpha) - beta, MAGMA_minproduct_D_IMAG(alpha));
            scale = MAGMA_minproduct_D_DIV( MAGMA_minproduct_D_ONE, alpha);
#endif
        }
    }

    // scale x
    __syncthreads();
    if ( xnorm != 0 && j < n-1)
        dx[j] = MAGMA_minproduct_D_MUL(dxi, scale);

    if (j<it){
        *( dA-it+j) = *(dx0-it+j);
        *(dx0-it+j) = MAGMA_minproduct_D_MAKE(0., 0.);
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
    
    The difference with LAPACK's dlarfg is that the norm of dx, and hance beta,
    are computed outside the routine and passed to it in dxnorm (array on the GPU).
*/
extern "C" void
magma_minproduct_dlarfgx_gpu(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx0,
    magma_minproductDouble_ptr dx,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr        dxnorm,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t iter)
{
    dim3 blocks((n+BLOCK_SIZE-1) / BLOCK_SIZE);
    dim3 threads( BLOCK_SIZE );
 
    magma_minproduct_dlarfgx_gpu_kernel<<< blocks, threads, 0, magma_minproduct_stream >>>( n, dx0, dx, dtau, dxnorm, dA, iter);
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
    
    The difference with LAPACK's dlarfg is that the norm of dx, and hance beta,
    are computed outside the routine and passed to it in dxnorm (array on the GPU).
*/
extern "C" void
magma_minproduct_dlarfgtx_gpu(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx0,
    magma_minproductDouble_ptr dx,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr        dxnorm,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t iter,
    magma_minproductDouble_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductDouble_ptr T,  magma_minproduct_int_t ldt,
    magma_minproductDouble_ptr dwork)
{
    /*  Generate the elementary reflector H(iter)  */
    magma_minproduct_dlarfgx_gpu(n, dx0, dx, dtau, dxnorm, dA, iter);
    
    if (iter==0) {
        double tt = MAGMA_minproduct_D_ONE;
        magma_minproductblas_dlacpy(Magma_minproductUpperLower, 1, 1, dtau, 1, T+iter+iter*ldt, 1);
        magma_minproduct_dsetmatrix(1,1, &tt,1, dx0,1);
    }
    else {
        /* Compute the iter-th column of T */
        magma_minproduct_dgemv_kernel3<<< iter, BLOCK_SIZE, 0, magma_minproduct_stream >>>( n, V, ldv, dx0, dwork, dtau );
        magma_minproduct_dtrmv_kernel2<<< iter, iter,       0, magma_minproduct_stream >>>( T, ldt, dwork, T+iter*ldt, dtau );
    }
}

//==============================================================================
