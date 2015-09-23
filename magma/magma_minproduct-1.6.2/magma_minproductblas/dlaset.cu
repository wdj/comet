/*
    -- MAGMA_minproduct (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Mark Gates
       @author Azzam Haidar
       
       @generated from zlaset.cu normal z -> d, Fri Mar 13 15:22:19 2015

*/
#include "common_magma_minproduct.h"
#include "batched_kernel_param.h"

// To deal with really large matrices, this launchs multiple super blocks,
// each with up to 64K-1 x 64K-1 thread blocks, which is up to 4194240 x 4194240 matrix with BLK=64.
// CUDA architecture 2.0 limits each grid dimension to 64K-1.
// Instances arose for vectors used by sparse matrices with M > 4194240, though N is small.
const magma_minproduct_int_t max_blocks = 65535;

// BLK_X and BLK_Y need to be equal for dlaset_q to deal with diag & offdiag
// when looping over super blocks.
// Formerly, BLK_X and BLK_Y could be different.
#define BLK_X 64
#define BLK_Y BLK_X

/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.

    Code similar to dlaset, dlacpy, dlag2s, clag2z, dgeadd.
*/
static __device__
void dlaset_full_device(
    int m, int n,
    double offdiag, double diag,
    double *A, int lda )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (below diag || above diag || offdiag == diag) */
    bool full = (iby + BLK_Y <= n && (ind >= iby + BLK_Y || ind + BLK_X <= iby || MAGMA_minproduct_D_EQUAL( offdiag, diag )));
    /* do only rows inside matrix */
    if ( ind < m ) {
        A += ind + iby*lda;
        if ( full ) {
            // full block-column, off-diagonal block or offdiag == diag
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                A[j*lda] = offdiag;
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                if ( iby+j == ind )
                    A[j*lda] = diag;
                else
                    A[j*lda] = offdiag;
            }
        }
    }
}


/*
    Similar to dlaset_full, but updates only the diagonal and below.
    Blocks that are fully above the diagonal exit immediately.

    Code similar to dlaset, dlacpy, zlat2c, clat2z.
*/
static __device__
void dlaset_lower_device(
    int m, int n,
    double offdiag, double diag,
    double *A, int lda )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (below diag) */
    bool full = (iby + BLK_Y <= n && (ind >= iby + BLK_Y));
    /* do only rows inside matrix, and blocks not above diag */
    if ( ind < m && ind + BLK_X > iby ) {
        A += ind + iby*lda;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                A[j*lda] = offdiag;
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                if ( iby+j == ind )
                    A[j*lda] = diag;
                else if ( ind > iby+j )
                    A[j*lda] = offdiag;
            }
        }
    }
}


/*
    Similar to dlaset_full, but updates only the diagonal and above.
    Blocks that are fully below the diagonal exit immediately.

    Code similar to dlaset, dlacpy, zlat2c, clat2z.
*/
static __device__
void dlaset_upper_device(
    int m, int n,
    double offdiag, double diag,
    double *A, int lda )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (above diag) */
    bool full = (iby + BLK_Y <= n && (ind + BLK_X <= iby));
    /* do only rows inside matrix, and blocks not below diag */
    if ( ind < m && ind < iby + BLK_Y ) {
        A += ind + iby*lda;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                A[j*lda] = offdiag;
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                if ( iby+j == ind )
                    A[j*lda] = diag;
                else if ( ind < iby+j )
                    A[j*lda] = offdiag;
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
/*
    kernel wrappers to call the device functions.
*/
__global__
void dlaset_full_kernel(
    int m, int n,
    double offdiag, double diag,
    double *dA, int ldda )
{
    dlaset_full_device(m, n, offdiag, diag, dA, ldda);
}

__global__
void dlaset_lower_kernel(
    int m, int n,
    double offdiag, double diag,
    double *dA, int ldda )
{
    dlaset_lower_device(m, n, offdiag, diag, dA, ldda);
}

__global__
void dlaset_upper_kernel(
    int m, int n,
    double offdiag, double diag,
    double *dA, int ldda )
{
    dlaset_upper_device(m, n, offdiag, diag, dA, ldda);
}


//////////////////////////////////////////////////////////////////////////////////////
/*
    kernel wrappers to call the device functions for the batched routine.
*/
__global__
void dlaset_full_kernel_batched(
    int m, int n,
    double offdiag, double diag,
    double **dAarray, int ldda )
{
    int batchid = blockIdx.z;
    dlaset_full_device(m, n, offdiag, diag, dAarray[batchid], ldda);
}

__global__
void dlaset_lower_kernel_batched(
    int m, int n,
    double offdiag, double diag,
    double **dAarray, int ldda )
{
    int batchid = blockIdx.z;
    dlaset_lower_device(m, n, offdiag, diag, dAarray[batchid], ldda);
}

__global__
void dlaset_upper_kernel_batched(
    int m, int n,
    double offdiag, double diag,
    double **dAarray, int ldda )
{
    int batchid = blockIdx.z;
    dlaset_upper_device(m, n, offdiag, diag, dAarray[batchid], ldda);
}


//////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    DLASET_Q initializes a 2-D array A to DIAG on the diagonal and
    OFFDIAG on the off-diagonals.
    
    This is the same as DLASET, but adds queue argument.
    
    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
            Specifies the part of the matrix dA to be set.
      -     = Magma_minproductUpper:      Upper triangular part
      -     = Magma_minproductLower:      Lower triangular part
      -     = Magma_minproductFull:       All of the matrix dA
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    offdiag DOUBLE_PRECISION
            The scalar OFFDIAG. (In LAPACK this is called ALPHA.)
    
    @param[in]
    diag    DOUBLE_PRECISION
            The scalar DIAG. (In LAPACK this is called BETA.)
    
    @param[in]
    dA      DOUBLE_PRECISION array, dimension (LDDA,N)
            The M-by-N matrix dA.
            If UPLO = Magma_minproductUpper, only the upper triangle or trapezoid is accessed;
            if UPLO = Magma_minproductLower, only the lower triangle or trapezoid is accessed.
            On exit, A(i,j) = OFFDIAG, 1 <= i <= m, 1 <= j <= n, i != j;
                     A(i,i) = DIAG,    1 <= i <= min(m,n)
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[in]
    queue   magma_minproduct_queue_t
            Queue to execute in.
    
    @ingroup magma_minproduct_daux2
    ********************************************************************/
extern "C"
void magma_minproductblas_dlaset_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double offdiag, double diag,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue)
{
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda)
    
    magma_minproduct_int_t info = 0;
    if ( uplo != Magma_minproductLower && uplo != Magma_minproductUpper && uplo != Magma_minproductFull )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < max(1,m) )
        info = -7;
    
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    if ( m == 0 || n == 0 ) {
        return;
    }
    
    assert( BLK_X == BLK_Y );
    const magma_minproduct_int_t super_NB = max_blocks*BLK_X;
    dim3 super_grid( magma_minproduct_ceildiv( m, super_NB ), magma_minproduct_ceildiv( n, super_NB ) );
    
    dim3 threads( BLK_X, 1 );
    dim3 grid;
    
    magma_minproduct_int_t mm, nn;
    if (uplo == Magma_minproductLower) {
        for( int i=0; i < super_grid.x; ++i ) {
            mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
            grid.x = magma_minproduct_ceildiv( mm, BLK_X );
            for( int j=0; j < super_grid.y && j <= i; ++j ) {  // from left to diagonal
                nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
                grid.y = magma_minproduct_ceildiv( nn, BLK_Y );
                if ( i == j ) {  // diagonal super block
                    dlaset_lower_kernel<<< grid, threads, 0, queue >>>
                        ( mm, nn, offdiag, diag, dA(i*super_NB, j*super_NB), ldda );
                }
                else {           // off diagonal super block
                    dlaset_full_kernel<<< grid, threads, 0, queue >>>
                        ( mm, nn, offdiag, offdiag, dA(i*super_NB, j*super_NB), ldda );
                }
            }
        }
    }
    else if (uplo == Magma_minproductUpper) {
        for( int i=0; i < super_grid.x; ++i ) {
            mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
            grid.x = magma_minproduct_ceildiv( mm, BLK_X );
            for( int j=i; j < super_grid.y; ++j ) {  // from diagonal to right
                nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
                grid.y = magma_minproduct_ceildiv( nn, BLK_Y );
                if ( i == j ) {  // diagonal super block
                    dlaset_upper_kernel<<< grid, threads, 0, queue >>>
                        ( mm, nn, offdiag, diag, dA(i*super_NB, j*super_NB), ldda );
                }
                else {           // off diagonal super block
                    dlaset_full_kernel<<< grid, threads, 0, queue >>>
                        ( mm, nn, offdiag, offdiag, dA(i*super_NB, j*super_NB), ldda );
                }
            }
        }
    }
    else {
        // if continuous in memory & set to zero, cudaMemset is faster.
        // TODO: use cudaMemset2D ?
        if ( m == ldda &&
             MAGMA_minproduct_D_EQUAL( offdiag, MAGMA_minproduct_D_ZERO ) &&
             MAGMA_minproduct_D_EQUAL( diag,    MAGMA_minproduct_D_ZERO ) )
        {
            size_t size = m*n;
            cudaError_t err = cudaMemsetAsync( dA, 0, size*sizeof(double), queue );
            assert( err == cudaSuccess );
        }
        else {
            for( int i=0; i < super_grid.x; ++i ) {
                mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
                grid.x = magma_minproduct_ceildiv( mm, BLK_X );
                for( int j=0; j < super_grid.y; ++j ) {  // full row
                    nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
                    grid.y = magma_minproduct_ceildiv( nn, BLK_Y );
                    if ( i == j ) {  // diagonal super block
                        dlaset_full_kernel<<< grid, threads, 0, queue >>>
                            ( mm, nn, offdiag, diag, dA(i*super_NB, j*super_NB), ldda );
                    }
                    else {           // off diagonal super block
                        dlaset_full_kernel<<< grid, threads, 0, queue >>>
                            ( mm, nn, offdiag, offdiag, dA(i*super_NB, j*super_NB), ldda );
                    }
                }
            }
        }
    }
}


/**
    @see magma_minproductblas_dlaset_q
    @ingroup magma_minproduct_daux2
    ********************************************************************/
extern "C"
void magma_minproductblas_dlaset(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double offdiag, double diag,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda )
{
    magma_minproductblas_dlaset_q( uplo, m, n, offdiag, diag, dA, ldda, magma_minproduct_stream );
}


////////////////////////////////////////////////////////////////////////////////////////
extern "C"
void magma_minproductblas_dlaset_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double offdiag, double diag,
    magma_minproductDouble_ptr dAarray[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{
    magma_minproduct_int_t info = 0;
    if ( uplo != Magma_minproductLower && uplo != Magma_minproductUpper && uplo != Magma_minproductFull )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < max(1,m) )
        info = -7;
    
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    if ( m == 0 || n == 0 ) {
        return;
    }
    
    dim3 threads( BLK_X, 1, 1 );
    dim3 grid( magma_minproduct_ceildiv( m, BLK_X ), magma_minproduct_ceildiv( n, BLK_Y ), batchCount );
    
    if (uplo == Magma_minproductLower) {
        dlaset_lower_kernel_batched<<< grid, threads, 0, queue >>> (m, n, offdiag, diag, dAarray, ldda);
    }
    else if (uplo == Magma_minproductUpper) {
        dlaset_upper_kernel_batched<<< grid, threads, 0, queue >>> (m, n, offdiag, diag, dAarray, ldda);
    }
    else {
        dlaset_full_kernel_batched<<< grid, threads, 0, queue >>> (m, n, offdiag, diag, dAarray, ldda);
    }
}
