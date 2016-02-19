/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zcaxpycp.cu mixed zc -> ds, Fri Jan 30 19:00:07 2015

*/
#include "common_magma_tally3.h"

#define NB 64

// adds   x += r (including conversion to double)  --and--
// copies w = b
// each thread does one index, x[i] and w[i]
__global__ void
dsaxpycp_kernel(
    int m, float *r, double *x,
    const double *b, double *w )
{
    const int i = threadIdx.x + blockIdx.x*NB;
    if ( i < m ) {
        x[i] = MAGMA_tally3_D_ADD( x[i], (double)( r[i] ) );
        w[i] = b[i];
    }
}


// adds   x += r  --and--
// copies r = b
// each thread does one index, x[i] and r[i]
__global__ void
daxpycp_kernel(
    int m, double *r, double *x,
    const double *b)
{
    const int i = threadIdx.x + blockIdx.x*NB;
    if ( i < m ) {
        x[i] = MAGMA_tally3_D_ADD( x[i], r[i] );
        r[i] = b[i];
    }
}


// ----------------------------------------------------------------------
// adds   x += r (including conversion to double)  --and--
// copies w = b
extern "C" void
magma_tally3blas_dsaxpycp_q(
    magma_tally3_int_t m,
    magma_tally3Float_ptr r,
    magma_tally3Double_ptr x,
    magma_tally3Double_const_ptr b,
    magma_tally3Double_ptr w,
    magma_tally3_queue_t queue )
{
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    dsaxpycp_kernel <<< grid, threads, 0, queue >>> ( m, r, x, b, w );
}


extern "C" void
magma_tally3blas_dsaxpycp(
    magma_tally3_int_t m,
    magma_tally3Float_ptr r,
    magma_tally3Double_ptr x,
    magma_tally3Double_const_ptr b,
    magma_tally3Double_ptr w)
{
    magma_tally3blas_dsaxpycp_q( m, r, x, b, w, magma_tally3_stream );
}


// ----------------------------------------------------------------------
// adds   x += r  --and--
// copies r = b
extern "C" void
magma_tally3blas_daxpycp_q(
    magma_tally3_int_t m,
    magma_tally3Double_ptr r,
    magma_tally3Double_ptr x,
    magma_tally3Double_const_ptr b,
    magma_tally3_queue_t queue )
{
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    daxpycp_kernel <<< grid, threads, 0, queue >>> ( m, r, x, b );
}

extern "C" void
magma_tally3blas_daxpycp(
    magma_tally3_int_t m,
    magma_tally3Double_ptr r,
    magma_tally3Double_ptr x,
    magma_tally3Double_const_ptr b)
{
    magma_tally3blas_daxpycp_q( m, r, x, b, magma_tally3_stream );
}
