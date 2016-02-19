/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds

*/
#include "common_magma_tally3.h"

#define NB 64

// adds   x += r (including conversion to double)  --and--
// copies w = b
// each thread does one index, x[i] and w[i]
__global__ void
zcaxpycp_kernel(
    int m, magma_tally3FloatComplex *r, magma_tally3DoubleComplex *x,
    const magma_tally3DoubleComplex *b, magma_tally3DoubleComplex *w )
{
    const int i = threadIdx.x + blockIdx.x*NB;
    if ( i < m ) {
        x[i] = MAGMA_tally3_Z_ADD( x[i], cuComplexFloatToDouble( r[i] ) );
        w[i] = b[i];
    }
}


// adds   x += r  --and--
// copies r = b
// each thread does one index, x[i] and r[i]
__global__ void
zaxpycp_kernel(
    int m, magma_tally3DoubleComplex *r, magma_tally3DoubleComplex *x,
    const magma_tally3DoubleComplex *b)
{
    const int i = threadIdx.x + blockIdx.x*NB;
    if ( i < m ) {
        x[i] = MAGMA_tally3_Z_ADD( x[i], r[i] );
        r[i] = b[i];
    }
}


// ----------------------------------------------------------------------
// adds   x += r (including conversion to double)  --and--
// copies w = b
extern "C" void
magma_tally3blas_zcaxpycp_q(
    magma_tally3_int_t m,
    magma_tally3FloatComplex_ptr r,
    magma_tally3DoubleComplex_ptr x,
    magma_tally3DoubleComplex_const_ptr b,
    magma_tally3DoubleComplex_ptr w,
    magma_tally3_queue_t queue )
{
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    zcaxpycp_kernel <<< grid, threads, 0, queue >>> ( m, r, x, b, w );
}


extern "C" void
magma_tally3blas_zcaxpycp(
    magma_tally3_int_t m,
    magma_tally3FloatComplex_ptr r,
    magma_tally3DoubleComplex_ptr x,
    magma_tally3DoubleComplex_const_ptr b,
    magma_tally3DoubleComplex_ptr w)
{
    magma_tally3blas_zcaxpycp_q( m, r, x, b, w, magma_tally3_stream );
}


// ----------------------------------------------------------------------
// adds   x += r  --and--
// copies r = b
extern "C" void
magma_tally3blas_zaxpycp_q(
    magma_tally3_int_t m,
    magma_tally3DoubleComplex_ptr r,
    magma_tally3DoubleComplex_ptr x,
    magma_tally3DoubleComplex_const_ptr b,
    magma_tally3_queue_t queue )
{
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    zaxpycp_kernel <<< grid, threads, 0, queue >>> ( m, r, x, b );
}

extern "C" void
magma_tally3blas_zaxpycp(
    magma_tally3_int_t m,
    magma_tally3DoubleComplex_ptr r,
    magma_tally3DoubleComplex_ptr x,
    magma_tally3DoubleComplex_const_ptr b)
{
    magma_tally3blas_zaxpycp_q( m, r, x, b, magma_tally3_stream );
}
