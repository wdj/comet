/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds

*/
#include "common_magma_minproduct.h"

#define NB 64

// adds   x += r (including conversion to double)  --and--
// copies w = b
// each thread does one index, x[i] and w[i]
__global__ void
zcaxpycp_kernel(
    int m, magma_minproductFloatComplex *r, magma_minproductDoubleComplex *x,
    const magma_minproductDoubleComplex *b, magma_minproductDoubleComplex *w )
{
    const int i = threadIdx.x + blockIdx.x*NB;
    if ( i < m ) {
        x[i] = MAGMA_minproduct_Z_ADD( x[i], cuComplexFloatToDouble( r[i] ) );
        w[i] = b[i];
    }
}


// adds   x += r  --and--
// copies r = b
// each thread does one index, x[i] and r[i]
__global__ void
zaxpycp_kernel(
    int m, magma_minproductDoubleComplex *r, magma_minproductDoubleComplex *x,
    const magma_minproductDoubleComplex *b)
{
    const int i = threadIdx.x + blockIdx.x*NB;
    if ( i < m ) {
        x[i] = MAGMA_minproduct_Z_ADD( x[i], r[i] );
        r[i] = b[i];
    }
}


// ----------------------------------------------------------------------
// adds   x += r (including conversion to double)  --and--
// copies w = b
extern "C" void
magma_minproductblas_zcaxpycp_q(
    magma_minproduct_int_t m,
    magma_minproductFloatComplex_ptr r,
    magma_minproductDoubleComplex_ptr x,
    magma_minproductDoubleComplex_const_ptr b,
    magma_minproductDoubleComplex_ptr w,
    magma_minproduct_queue_t queue )
{
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    zcaxpycp_kernel <<< grid, threads, 0, queue >>> ( m, r, x, b, w );
}


extern "C" void
magma_minproductblas_zcaxpycp(
    magma_minproduct_int_t m,
    magma_minproductFloatComplex_ptr r,
    magma_minproductDoubleComplex_ptr x,
    magma_minproductDoubleComplex_const_ptr b,
    magma_minproductDoubleComplex_ptr w)
{
    magma_minproductblas_zcaxpycp_q( m, r, x, b, w, magma_minproduct_stream );
}


// ----------------------------------------------------------------------
// adds   x += r  --and--
// copies r = b
extern "C" void
magma_minproductblas_zaxpycp_q(
    magma_minproduct_int_t m,
    magma_minproductDoubleComplex_ptr r,
    magma_minproductDoubleComplex_ptr x,
    magma_minproductDoubleComplex_const_ptr b,
    magma_minproduct_queue_t queue )
{
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    zaxpycp_kernel <<< grid, threads, 0, queue >>> ( m, r, x, b );
}

extern "C" void
magma_minproductblas_zaxpycp(
    magma_minproduct_int_t m,
    magma_minproductDoubleComplex_ptr r,
    magma_minproductDoubleComplex_ptr x,
    magma_minproductDoubleComplex_const_ptr b)
{
    magma_minproductblas_zaxpycp_q( m, r, x, b, magma_minproduct_stream );
}
