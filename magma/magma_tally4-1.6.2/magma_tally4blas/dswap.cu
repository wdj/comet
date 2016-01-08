/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Mark Gates

       @generated from zswap.cu normal z -> d, Fri Jan 30 19:00:09 2015

*/
#include "common_magma_tally4.h"

#define NB 64


/* Vector is divided into ceil(n/nb) blocks.
   Each thread swaps one element, x[tid] <---> y[tid].
*/
__global__ void dswap_kernel(
    int n,
    double *x, int incx,
    double *y, int incy )
{
    double tmp;
    int ind = threadIdx.x + blockDim.x*blockIdx.x;
    if ( ind < n ) {
        x += ind*incx;
        y += ind*incy;
        tmp = *x;
        *x  = *y;
        *y  = tmp;
    }
}


/**
    Purpose:
    =============
    Swap vector x and y; \f$ x <-> y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally4_dblas1
    ********************************************************************/
extern "C" void 
magma_tally4blas_dswap_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dx, magma_tally4_int_t incx, 
    magma_tally4Double_ptr dy, magma_tally4_int_t incy,
    magma_tally4_queue_t queue )
{
    dim3 grid( (n+NB-1) / NB );
    dim3 threads( NB );
    dswap_kernel<<< grid, threads, 0, queue >>>( n, dx, incx, dy, incy );
}


/**
    @see magma_tally4blas_dswap_q
    @ingroup magma_tally4_dblas1
    ********************************************************************/
extern "C" void 
magma_tally4blas_dswap(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dx, magma_tally4_int_t incx, 
    magma_tally4Double_ptr dy, magma_tally4_int_t incy)
{
    magma_tally4blas_dswap_q( n, dx, incx, dy, incy, magma_tally4_stream );
}
