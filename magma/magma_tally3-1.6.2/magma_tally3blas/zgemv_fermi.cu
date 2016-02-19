/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Mark Gates

       @precisions normal z -> s d c
*/
#include "common_magma_tally3.h"
#include "commonblas_z.h"
#include "magma_tally3_templates.h"

#define PRECISION_z

#define BLK_X 128
#define BLK_Y 128

/* Compute y = alpha*A*x + beta*y.
 * Each thread block does a BLK_X x N block row of A.
 * Each thread goes across one row, accumulating dot product of row ind and x into res.
 * This simple implementation loads x directly, relying on the cache,
 * without using shared memory.
 */
__global__ void
zgemvn_kernel1_fermi(
    int m, int n, magma_tally3DoubleComplex alpha,
    const magma_tally3DoubleComplex * __restrict__ A, int lda,
    const magma_tally3DoubleComplex * __restrict__ x, int incx, magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex       * __restrict__ y, int incy)
{
#if (__CUDA_ARCH__ >= 200)
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    if ( ind < m ) {
        A += ind;
        
        magma_tally3DoubleComplex res = MAGMA_tally3_Z_ZERO;
        
        #pragma unroll
        for(int j=0; j < n; j++) {
            res += A[j*lda] * x[j*incx];
        }
        
        y[ind*incy] = alpha*res + beta*y[ind*incy];
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}


/* Compute y = alpha*A*x + beta*y.
 * Each thread block does a BLK_X x N block row of A.
 * Each thread goes across one row, accumulating dot product of row ind and x into res.
 * This implementation loads BLK_Y elements into sx, then multiplies
 * BLK_Y columns of A*sx.
 */
__global__ void
zgemvn_kernel2_fermi(
    int m, int n, magma_tally3DoubleComplex alpha,
    const magma_tally3DoubleComplex * __restrict__ A, int lda,
    const magma_tally3DoubleComplex * __restrict__ x, int incx, magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex       * __restrict__ y, int incy)
{
#if (__CUDA_ARCH__ >= 200)
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    
    // threads past last row redundantly work on last row
    A += min( ind, m-1 );
    x += threadIdx.x*incx;
    
    magma_tally3DoubleComplex res = MAGMA_tally3_Z_ZERO;
    
    __shared__ magma_tally3DoubleComplex sx[BLK_Y];
    
    // full block-columns
    int nfull = (n / BLK_Y) * BLK_Y;
    for( int j=0; j < nfull; j += BLK_Y ) {
        // load BLK_Y elements of x into sx
        sx[threadIdx.x] = x[0];
        x += BLK_Y*incx;
        __syncthreads();
        
        // multiply A*sx
        #pragma unroll
        for(int j2=0; j2 < BLK_Y; j2++) {
            res += A[0] * sx[j2];
            A += lda;
        }
        __syncthreads();
    }
    
    // last, partial block-column
    // load remaining npart elements of x into sx
    int npart = n % BLK_Y;
    if ( threadIdx.x < npart ) {
        sx[threadIdx.x] = x[0];
    }
    else {
        sx[threadIdx.x] = MAGMA_tally3_Z_ZERO;
    }
    __syncthreads();
        
    // multiply A*sx
    #pragma unroll
    for(int j2=0; j2 < npart; j2++) {
        res += A[0]*sx[j2];
        A += lda;
    }
    
    if ( ind < m ) {
        y[ind*incy] = alpha*res + beta*y[ind*incy];
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}


/* Compute y = alpha * A^T * x + beta*y.
 * Each thread block does one column of A (i.e., one row of A^T).
 * Each thread does a partial sum, then collectively they do a reduction.
 */
__global__ void
zgemvt_kernel_fermi(
    int m, int n, magma_tally3DoubleComplex alpha,
    const magma_tally3DoubleComplex * __restrict__ A, int lda,
    const magma_tally3DoubleComplex * __restrict__ x, int incx, magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex       * __restrict__ y, int incy)
{
#if (__CUDA_ARCH__ >= 200)
    int tx = threadIdx.x;

    __shared__ magma_tally3DoubleComplex sdata[BLK_X];

    magma_tally3DoubleComplex res = MAGMA_tally3_Z_ZERO;
    
    A += blockIdx.y*lda + threadIdx.x;
 
    // partial sums
    int mfull = (m / BLK_X) * BLK_X;
    for(int i=0; i < mfull; i += BLK_X) {
        res += A[i] * x[tx + i];
    }
    if ( tx + mfull < m ) {
        res += A[mfull] * x[tx + mfull];
    }
    sdata[tx] = res;

    // tree reduction of partial sums,
    // from BLK_X sums to ... 128 to 64 to 32 ... to 1 sum in sdata[0]
    magma_tally3_sum_reduce< BLK_X >( tx, sdata );

    if ( tx == 0 ) {
        y[blockIdx.y*incy] = alpha*sdata[0] + beta*y[blockIdx.y*incy];
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}


/* Compute y = alpha * A^H * x + beta*y.
 * Same as zgemvt_kernel_fermi but conjugates entries of A.
 */
__global__ void
zgemvc_kernel_fermi(
    int m, int n, magma_tally3DoubleComplex alpha,
    const magma_tally3DoubleComplex * __restrict__ A, int lda,
    const magma_tally3DoubleComplex * __restrict__ x, int incx, magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex       * __restrict__ y, int incy)
{
#if (__CUDA_ARCH__ >= 200)
    int tx = threadIdx.x;

    __shared__ magma_tally3DoubleComplex sdata[BLK_X];

    magma_tally3DoubleComplex res = MAGMA_tally3_Z_ZERO;
    
    A += blockIdx.y*lda + threadIdx.x;
 
    // partial sums
    int mfull = (m / BLK_X) * BLK_X;
    for(int i=0; i < mfull; i += BLK_X) {
        res += conj(A[i]) * x[tx + i];
    }
    if ( tx + mfull < m ) {
        res += conj(A[mfull]) * x[tx + mfull];
    }
    sdata[tx] = res;

    // tree reduction of partial sums,
    // from BLK_X sums to ... 128 to 64 to 32 ... to 1 sum in sdata[0]
    magma_tally3_sum_reduce< BLK_X >( tx, sdata );

    if ( tx == 0 ) {
        y[blockIdx.y*incy] = alpha*sdata[0] + beta*y[blockIdx.y*incy];
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}


/**
    Purpose
    -------
    ZGEMV performs one of the matrix-vector operations
    
        y := alpha*A*x    + beta*y,   or
        y := alpha*A**T*x + beta*y,   or
        y := alpha*A**H*x + beta*y,
    
    where alpha and beta are scalars, x and y are vectors and A is an
    m by n matrix.

    Arguments
    ----------
    @param[in]
    trans   magma_tally3_trans_t
            On entry, TRANS specifies the operation to be performed as
            follows:
      -     = Magma_tally3NoTrans:    y := alpha*A  *x + beta*y
      -     = Magma_tally3Trans:      y := alpha*A^T*x + beta*y
      -     = Magma_tally3ConjTrans:  y := alpha*A^H*x + beta*y

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
    dx      COMPLEX_16 array of dimension
            n if trans == Magma_tally3NoTrans
            m if trans == Magma_tally3Trans or Magma_tally3ConjTrans
     
    @param[in]
    incx    Specifies the increment for the elements of X.
            INCX must not be zero.
  
    @param[in]
    beta    DOUBLE REAL
            On entry, BETA specifies the scalar beta. When BETA is
            supplied as zero then Y need not be set on input.

    @param[out]
    dy      DOUBLE PRECISION array of dimension
            m if trans == Magma_tally3NoTrans
            n if trans == Magma_tally3Trans or Magma_tally3ConjTrans

    @param[in]
    incy    Specifies the increment for the elements of Y.
            INCY must not be zero.

    @ingroup magma_tally3_dblas2
    ********************************************************************/
extern "C" void
magma_tally3blas_zgemv(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy, magma_tally3_int_t incy)
{
    magma_tally3_int_t info = 0;
    if ( trans != Magma_tally3NoTrans && trans != Magma_tally3Trans && trans != Magma_tally3ConjTrans )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < m )
        info = -6;
    else if ( incx == 0 )
        info = -8;
    else if ( incy == 0 )
        info = -11;
    
    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally3_int_t arch = magma_tally3_getdevice_arch();
    if ( arch < 200  ) {
        // --------------------
        // call CUDA ARCH 1.x version
        // magma_tally3blas for [sd] precisions, cublas for [zc] precisions.
        #if defined(PRECISION_z) || defined(PRECISION_c)
        magma_tally3_zgemv( trans, m, n, alpha, dA, ldda, dx, incx, beta, dy, incy );
        #else
        magma_tally3blas_zgemv_tesla( trans, m, n, alpha, dA, ldda, dx, incx, beta, dy, incy );
        #endif
        return;
    }
    
    // --------------------
    // CUDA ARCH 2.x (Fermi) version
    if ( trans == Magma_tally3NoTrans ) {
        dim3 grid( (m - 1)/BLK_X + 1 );
        dim3 threads( BLK_X, 1, 1 );
        zgemvn_kernel1_fermi<<< grid, threads, 0, magma_tally3_stream >>>
            ( m, n, alpha, dA, ldda, dx, incx, beta, dy, incy );
    }
    else if ( trans == Magma_tally3Trans ) {
        dim3 grid    ( 1, n, 1 );
        dim3 threads ( BLK_X, 1, 1 );
        zgemvt_kernel_fermi<<< grid, threads, 0, magma_tally3_stream >>>
            ( m, n, alpha, dA, ldda, dx, incx, beta, dy, incy );
    }
    else if ( trans == Magma_tally3ConjTrans ) {
        dim3 grid    ( 1, n, 1 );
        dim3 threads ( BLK_X, 1, 1 );
        zgemvc_kernel_fermi<<< grid, threads, 0, magma_tally3_stream >>>
            ( m, n, alpha, dA, ldda, dx, incx, beta, dy, incy );
    }
}
