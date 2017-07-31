/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Stan Tomov

       @precisions normal z -> s d c
*/
#include "common_magma_tally2.h"
#include "commonblas_z.h"

#define PRECISION_z

#define num_threads 256


__global__ void
zgemv_conjv_kernel(
    int m, int n, magma_tally2DoubleComplex alpha,
    const magma_tally2DoubleComplex * __restrict__ A, int lda,
    const magma_tally2DoubleComplex * __restrict__ x, int incx, magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex *       __restrict__ y, int incy)
{
    int ind = blockIdx.x*num_threads + threadIdx.x;
    
    A += ind;

    if ( ind < m ) {
        magma_tally2DoubleComplex res = MAGMA_tally2_Z_ZERO;
        
        #pragma unroll
        for( int i=0; i < n; i ++ ) {
            res += A[0] * MAGMA_tally2_Z_CNJG(x[0]);
            A += lda;
            x += incx;
        }
        
        y[ind*incy] = alpha * res + beta * y[ind*incy];
    }
}


/**
    Purpose
    -------
    ZGEMV_CONJV performs the matrix-vector operation
    
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
    alpha   COMPLEX_16
            On entry, ALPHA specifies the scalar alpha.

    @param[in]
    dA      COMPLEX_16 array of dimension ( LDA, n ) on the GPU.

    @param[in]
    lda     INTEGER
            LDA specifies the leading dimension of A.

    @param[in]
    dx      COMPLEX_16 array of dimension n

    @param[in]
    incx    Specifies the increment for the elements of X.
            INCX must not be zero.

    @param[in]
    beta    DOUBLE REAL
            On entry, BETA specifies the scalar beta. When BETA is
            supplied as zero then Y need not be set on input.

    @param[out]
    dy      DOUBLE PRECISION array of dimension m

    @param[in]
    incy    Specifies the increment for the elements of Y.
            INCY must not be zero.

    @ingroup magma_tally2_zblas2
    ********************************************************************/
extern "C" void
magma_tally2blas_zgemv_conjv(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy, magma_tally2_int_t incy)
{
    magma_tally2_int_t info = 0;
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
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally2_int_t blocks = (m - 1)/num_threads + 1;
    dim3 grid(blocks, 1, 1);
    dim3 threads(num_threads, 1, 1);

    zgemv_conjv_kernel<<< grid, threads, 0, magma_tally2_stream >>>
            (m, n, alpha, dA, ldda, dx, incx, beta, dy, incy);

}

#undef num_threads
