/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Stan Tomov

       @generated from zgemv_conjv.cu normal z -> c, Fri Jan 30 19:00:08 2015
*/
#include "common_magma_minproduct.h"
#include "commonblas_c.h"

#define PRECISION_c

#define num_threads 256


__global__ void
cgemv_conjv_kernel(
    int m, int n, magma_minproductFloatComplex alpha,
    const magma_minproductFloatComplex * __restrict__ A, int lda,
    const magma_minproductFloatComplex * __restrict__ x, int incx, magma_minproductFloatComplex beta,
    magma_minproductFloatComplex *       __restrict__ y, int incy)
{
    int ind = blockIdx.x*num_threads + threadIdx.x;
    
    A += ind;

    if ( ind < m ) {
        magma_minproductFloatComplex res = MAGMA_minproduct_C_ZERO;
        
        #pragma unroll
        for( int i=0; i < n; i ++ ) {
            res += A[0] * MAGMA_minproduct_C_CNJG(x[0]);
            A += lda;
            x += incx;
        }
        
        y[ind*incy] = alpha * res + beta * y[ind*incy];
    }
}


/**
    Purpose
    -------
    CGEMV_CONJV performs the matrix-vector operation
    
        y := alpha*A*conj(x)    + beta*y, 
    
    where alpha and beta are scalars, x and y are vectors and A is an
    m by n matrix.

    Arguments
    ----------
    @param[in]
    m       INTEGER
            On entry, m specifies the number of rows of the matrix A.

    @param[in]
    n       INTEGER
            On entry, n specifies the number of columns of the matrix A

    @param[in]
    alpha   COMPLEX
            On entry, ALPHA specifies the scalar alpha.

    @param[in]
    dA      COMPLEX array of dimension ( LDA, n ) on the GPU.

    @param[in]
    lda     INTEGER
            LDA specifies the leading dimension of A.

    @param[in]
    dx      COMPLEX array of dimension n

    @param[in]
    incx    Specifies the increment for the elements of X.
            INCX must not be zero.

    @param[in]
    beta    DOUBLE REAL
            On entry, BETA specifies the scalar beta. When BETA is
            supplied as zero then Y need not be set on input.

    @param[out]
    dy      REAL array of dimension m

    @param[in]
    incy    Specifies the increment for the elements of Y.
            INCY must not be zero.

    @ingroup magma_minproduct_cblas2
    ********************************************************************/
extern "C" void
magma_minproductblas_cgemv_conjv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy, magma_minproduct_int_t incy)
{
    magma_minproduct_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < m )
        info = -5;
    else if ( incx == 0 )
        info = -7;
    else if ( incy == 0 )
        info = -10;
    
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_minproduct_int_t blocks = (m - 1)/num_threads + 1;
    dim3 grid(blocks, 1, 1);
    dim3 threads(num_threads, 1, 1);

    cgemv_conjv_kernel<<< grid, threads, 0, magma_minproduct_stream >>>
            (m, n, alpha, dA, ldda, dx, incx, beta, dy, incy);

}

#undef num_threads
