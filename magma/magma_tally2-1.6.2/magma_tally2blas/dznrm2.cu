/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"
#include "commonblas_z.h"
#include "magma_tally2_templates.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE  512
#define BLOCK_SIZEx  32
#define BLOCK_SIZEy  16

#define PRECISION_z


//==============================================================================

__global__ void
magma_tally2blas_dznrm2_kernel( int m, magma_tally2DoubleComplex *dA, int ldda, double *dxnorm )
{
    const int tx = threadIdx.x;
    magma_tally2DoubleComplex *dx = dA + blockIdx.x * ldda;

    __shared__ double sum[ BLOCK_SIZE ];
    double re, lsum;

    // get norm of dx
    lsum = 0;
    for( int j = tx; j < m; j += BLOCK_SIZE ) {
        #if (defined(PRECISION_s) || defined(PRECISION_d))
        re = dx[j];
        lsum += re*re;
        #else
        re = MAGMA_tally2_Z_REAL( dx[j] );
        double im = MAGMA_tally2_Z_IMAG( dx[j] );
        lsum += re*re + im*im;
        #endif
    }
    sum[tx] = lsum;
    magma_tally2_sum_reduce< BLOCK_SIZE >( tx, sum );
    
    if (tx==0)
        dxnorm[blockIdx.x] = sqrt(sum[0]);
}


//==============================================================================
__global__ void
magma_tally2blas_dznrm2_check_kernel( int m, magma_tally2DoubleComplex *dA, int ldda, double *dxnorm, 
                               double *lsticc )
{
    const int tx = threadIdx.x;
    magma_tally2DoubleComplex *dx = dA + blockIdx.x * ldda;

    __shared__ double sum[ BLOCK_SIZE ];
    double re, lsum;

    // get norm of dx only if lsticc[blockIdx+1] != 0
    if ( lsticc[blockIdx.x + 1] == 0 )
        return;

    lsum = 0;
    for( int j = tx; j < m; j += BLOCK_SIZE ) {
        #if (defined(PRECISION_s) || defined(PRECISION_d))
        re = dx[j];
        lsum += re*re;
        #else
        re = MAGMA_tally2_Z_REAL( dx[j] );
        double im = MAGMA_tally2_Z_IMAG( dx[j] );
        lsum += re*re + im*im;
        #endif
    }
    sum[tx] = lsum;
    magma_tally2_sum_reduce< BLOCK_SIZE >( tx, sum );
    
    if (tx==0)
        dxnorm[blockIdx.x] = sqrt(sum[0]);
}

extern "C" void
magma_tally2blas_dznrm2_check(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2Double_ptr dxnorm,
    magma_tally2Double_ptr dlsticc) 
{
    dim3  blocks( n );
    dim3 threads( BLOCK_SIZE );
    
    magma_tally2blas_dznrm2_check_kernel<<< blocks, threads >>>( m, dA, ldda, dxnorm, dlsticc );
}


//==============================================================================
__global__ void
magma_tally2blas_dznrm2_smkernel( int m, int n, magma_tally2DoubleComplex *dA, int ldda,
                           double *dxnorm )
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    __shared__ double sum[ BLOCK_SIZEx ][ BLOCK_SIZEy + 1];
    double re, lsum;

    for( int k = ty; k < n; k += BLOCK_SIZEy ) {
        magma_tally2DoubleComplex *dx = dA + k * ldda;

        // get norm of dx
        lsum = 0;
        for( int j = tx; j < m; j += BLOCK_SIZEx ) {
            #if (defined(PRECISION_s) || defined(PRECISION_d))
            re = dx[j];
            lsum += re*re;
            #else
            re = MAGMA_tally2_Z_REAL( dx[j] );
            double im = MAGMA_tally2_Z_IMAG( dx[j] );
            lsum += re*re + im*im;
            #endif
        }
        sum[tx][ty] = lsum;
        magma_tally2_sum_reduce_2d< BLOCK_SIZEx, BLOCK_SIZEy+1 >( tx, ty, sum );

        if (tx == 0)
            dxnorm[k] = sqrt(sum[0][ty]);
        __syncthreads();
    }
}


//==============================================================================
/*
    Compute the dznrm2 of each column of m-by-n matrix dA.
    The resulting norms are written in the dxnorm array.
    This routine uses only one SM (block).
*/
extern "C" void
magma_tally2blas_dznrm2_sm(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double *dxnorm)
{
    dim3  blocks( 1 );
    dim3 threads( BLOCK_SIZEx, BLOCK_SIZEy );

    magma_tally2blas_dznrm2_smkernel<<< blocks, threads, 0, magma_tally2_stream >>>( m, n, dA, ldda, dxnorm );
}

//==============================================================================
extern "C"
__global__ void
magma_tally2_dznrm2_adjust_kernel(double *xnorm, magma_tally2DoubleComplex *c)
{
    const int tx = threadIdx.x;

    __shared__ double sum[ BLOCK_SIZE ];
    double temp;

    temp = MAGMA_tally2_Z_ABS( c[tx] ) / xnorm[0];
    sum[tx] = -temp * temp;
    magma_tally2_sum_reduce_n( blockDim.x, tx, sum );

    __syncthreads();
    if (tx == 0)
        xnorm[0] = xnorm[0] * sqrt(1+sum[0]);
}


/*
    Adjust the norm of c to give the norm of c[k+1:], assuming that
    c was changed with orthogonal transformations.
*/
extern "C" void
magma_tally2blas_dznrm2_adjust(magma_tally2_int_t k, magma_tally2Double_ptr dxnorm, magma_tally2DoubleComplex_ptr dc)
{
    magma_tally2_dznrm2_adjust_kernel<<< 1, k, 0, magma_tally2_stream >>> (dxnorm, dc);
}

//==============================================================================

#define BS 256

__global__ void
magma_tally2_dznrm2_row_check_adjust_kernel(
    int n, double tol, double *xnorm, double *xnorm2, 
    magma_tally2DoubleComplex *C, int ldc, double *lsticc)
{
    const int tx = threadIdx.x + blockIdx.x*BS;
    lsticc[tx+1] = 0;

    if (tx < n) {
        double temp = MAGMA_tally2_Z_ABS( C[tx*ldc] ) / xnorm[tx];
        temp = max( 0.0, ((1.0 + temp) * (1.0 - temp)) );
        
        
        double temp2 = xnorm[tx] / xnorm2[tx];
        temp2 = temp * (temp2 * temp2);
        
        if (temp2 <= tol) {
            lsticc[tx+1] = 1;
        } else {
            xnorm[tx] *= sqrt(temp);
        }
    }
    if (tx == 0)
        lsticc[0] = 0;
    magma_tally2_sum_reduce_n( blockDim.x, tx, lsticc );
}

/*
    Adjust the norm of C[,1:k] to give the norm of C[k+1:,1:k], assuming that
    C was changed with orthogonal transformations.
    It also do checks for QP3
*/
extern "C" void
magma_tally2blas_dznrm2_row_check_adjust(
    magma_tally2_int_t k, double tol,
    magma_tally2Double_ptr dxnorm,
    magma_tally2Double_ptr dxnorm2, 
    magma_tally2DoubleComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2Double_ptr dlsticc)
{
    int nblocks = (k+BS-1)/BS;
    magma_tally2_dznrm2_row_check_adjust_kernel<<< nblocks, BS >>> (k, tol, dxnorm, dxnorm2, dC, lddc, dlsticc);
}

//==============================================================================

/*
    Compute the dznrm2 of each column of m-by-n matrix dA.
    The resulting norms are written in the dxnorm array. 
    The computation can be done using n blocks (default) or on one SM (commented).
*/
extern "C" void
magma_tally2blas_dznrm2_cols(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2Double_ptr dxnorm) 
{
    dim3  blocks( n );
    dim3 threads( BLOCK_SIZE );
    
    magma_tally2blas_dznrm2_kernel<<< blocks, threads, 0, magma_tally2_stream >>>( m, dA, ldda, dxnorm );

    // The following would do the computation on one SM
    // magma_tally2blas_dznrm2_sm(m, n, dA, ldda, dxnorm);
}

//==============================================================================
