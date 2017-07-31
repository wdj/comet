/*
    -- MAGMA_tally2 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Mark Gates
       @author Azzam Haidar
       
       @generated from zlacpy.cu normal z -> c, Fri Mar 13 15:22:41 2015

*/
#include "common_magma_tally2.h"

// To deal with really large matrices, this launchs multiple super blocks,
// each with up to 64K-1 x 64K-1 thread blocks, which is up to 4194240 x 4194240 matrix with BLK=64.
// CUDA architecture 2.0 limits each grid dimension to 64K-1.
// Instances arose for vectors used by sparse matrices with M > 4194240, though N is small.
const magma_tally2_int_t max_blocks = 65535;

// BLK_X and BLK_Y need to be equal for claset_q to deal with diag & offdiag
// when looping over super blocks.
// Formerly, BLK_X and BLK_Y could be different.
#define BLK_X 64
#define BLK_Y BLK_X

/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.

    Code similar to claset, clacpy, clag2z, clag2z, cgeadd.
*/
static __device__
void clacpy_full_device(
    int m, int n,
    const magma_tally2FloatComplex *dA, int ldda,
    magma_tally2FloatComplex       *dB, int lddb )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);
    /* do only rows inside matrix */
    if ( ind < m ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // partial block-column
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
    }
}


/*
    Similar to clacpy_full, but updates only the diagonal and below.
    Blocks that are fully above the diagonal exit immediately.

    Code similar to claset, clacpy, zlat2c, clat2z.
*/
static __device__
void clacpy_lower_device(
    int m, int n,
    const magma_tally2FloatComplex *dA, int ldda,
    magma_tally2FloatComplex       *dB, int lddb )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (below diag) */
    bool full = (iby + BLK_Y <= n && (ind >= iby + BLK_Y));
    /* do only rows inside matrix, and blocks not above diag */
    if ( ind < m && ind + BLK_X > iby ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n && ind >= iby+j; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
    }
}


/*
    Similar to clacpy_full, but updates only the diagonal and above.
    Blocks that are fully below the diagonal exit immediately.

    Code similar to claset, clacpy, zlat2c, clat2z.
*/
static __device__
void clacpy_upper_device(
    int m, int n,
    const magma_tally2FloatComplex *dA, int ldda,
    magma_tally2FloatComplex       *dB, int lddb )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (above diag) */
    bool full = (iby + BLK_Y <= n && (ind + BLK_X <= iby));
    /* do only rows inside matrix, and blocks not below diag */
    if ( ind < m && ind < iby + BLK_Y ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = dA[j*ldda];
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                if ( ind <= iby+j ) {
                    dB[j*lddb] = dA[j*ldda];
                }
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
/*
    kernel wrappers to call the device functions.
*/
__global__
void clacpy_full_kernel(
    int m, int n,
    const magma_tally2FloatComplex *dA, int ldda,
    magma_tally2FloatComplex       *dB, int lddb )
{
    clacpy_full_device(m, n, dA, ldda, dB, lddb);
}

__global__
void clacpy_lower_kernel(
    int m, int n,
    const magma_tally2FloatComplex *dA, int ldda,
    magma_tally2FloatComplex       *dB, int lddb )
{
    clacpy_lower_device(m, n, dA, ldda, dB, lddb);
}

__global__
void clacpy_upper_kernel(
    int m, int n,
    const magma_tally2FloatComplex *dA, int ldda,
    magma_tally2FloatComplex       *dB, int lddb )
{
    clacpy_upper_device(m, n, dA, ldda, dB, lddb);
}


//////////////////////////////////////////////////////////////////////////////////////
/*
    kernel wrappers to call the device functions for the batched routine.
*/
__global__
void clacpy_full_kernel_batched(
    int m, int n,
    magma_tally2FloatComplex const * const *dAarray, int ldda,
    magma_tally2FloatComplex **dBarray, int lddb )
{
    int batchid = blockIdx.z;
    clacpy_full_device(m, n, dAarray[batchid], ldda, dBarray[batchid], lddb);
}

__global__
void clacpy_lower_kernel_batched(
    int m, int n,
    magma_tally2FloatComplex const * const *dAarray, int ldda,
    magma_tally2FloatComplex **dBarray, int lddb )
{
    int batchid = blockIdx.z;
    clacpy_lower_device(m, n, dAarray[batchid], ldda, dBarray[batchid], lddb);
}

__global__
void clacpy_upper_kernel_batched(
    int m, int n,
    magma_tally2FloatComplex const * const *dAarray, int ldda,
    magma_tally2FloatComplex **dBarray, int lddb )
{
    int batchid = blockIdx.z;
    clacpy_upper_device(m, n, dAarray[batchid], ldda, dBarray[batchid], lddb);
}


//////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    CLACPY_Q copies all or part of a two-dimensional matrix dA to another
    matrix dB.
    
    This is the same as CLACPY, but adds queue argument.
    
    Arguments
    ---------
    @param[in]
    uplo    magma_tally2_uplo_t
            Specifies the part of the matrix dA to be copied to dB.
      -     = Magma_tally2Upper:      Upper triangular part
      -     = Magma_tally2Lower:      Lower triangular part
      -     = Magma_tally2Full:       All of the matrix dA
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    dA      COMPLEX array, dimension (LDDA,N)
            The M-by-N matrix dA.
            If UPLO = Magma_tally2Upper, only the upper triangle or trapezoid is accessed;
            if UPLO = Magma_tally2Lower, only the lower triangle or trapezoid is accessed.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[out]
    dB      COMPLEX array, dimension (LDDB,N)
            The M-by-N matrix dB.
            On exit, dB = dA in the locations specified by UPLO.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  LDDB >= max(1,M).
    
    @param[in]
    queue   magma_tally2_queue_t
            Queue to execute in.

    @ingroup magma_tally2_caux2
    ********************************************************************/
extern "C" void
magma_tally2blas_clacpy_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue )
{
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda)
    #define dB(i_, j_) (dB + (i_) + (j_)*lddb)
    
    magma_tally2_int_t info = 0;
    if ( uplo != Magma_tally2Lower && uplo != Magma_tally2Upper && uplo != Magma_tally2Full )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < max(1,m))
        info = -5;
    else if ( lddb < max(1,m))
        info = -7;
    
    if ( info != 0 ) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    if ( m == 0 || n == 0 ) {
        return;
    }
    
    assert( BLK_X == BLK_Y );
    const magma_tally2_int_t super_NB = max_blocks*BLK_X;
    dim3 super_grid( magma_tally2_ceildiv( m, super_NB ), magma_tally2_ceildiv( n, super_NB ) );
    
    dim3 threads( BLK_X, 1 );
    dim3 grid;
    
    magma_tally2_int_t mm, nn;
    if ( uplo == Magma_tally2Lower ) {
        for( int i=0; i < super_grid.x; ++i ) {
            mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
            grid.x = magma_tally2_ceildiv( mm, BLK_X );
            for( int j=0; j < super_grid.y && j <= i; ++j ) {  // from left to diagonal
                nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
                grid.y = magma_tally2_ceildiv( nn, BLK_Y );
                if ( i == j ) {  // diagonal super block
                    clacpy_lower_kernel<<< grid, threads, 0, queue >>>
                        ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
                }
                else {           // off diagonal super block
                    clacpy_full_kernel <<< grid, threads, 0, queue >>>
                        ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
                }
            }
        }
    }
    else if ( uplo == Magma_tally2Upper ) {
        for( int i=0; i < super_grid.x; ++i ) {
            mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
            grid.x = magma_tally2_ceildiv( mm, BLK_X );
            for( int j=i; j < super_grid.y; ++j ) {  // from diagonal to right
                nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
                grid.y = magma_tally2_ceildiv( nn, BLK_Y );
                if ( i == j ) {  // diagonal super block
                    clacpy_upper_kernel<<< grid, threads, 0, queue >>>
                        ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
                }
                else {           // off diagonal super block
                    clacpy_full_kernel <<< grid, threads, 0, queue >>>
                        ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
                }
            }
        }
    }
    else {
        // TODO: use cudaMemcpy or cudaMemcpy2D ?
        for( int i=0; i < super_grid.x; ++i ) {
            mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
            grid.x = magma_tally2_ceildiv( mm, BLK_X );
            for( int j=0; j < super_grid.y; ++j ) {  // full row
                nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
                grid.y = magma_tally2_ceildiv( nn, BLK_Y );
                clacpy_full_kernel <<< grid, threads, 0, queue >>>
                    ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
            }
        }
    }
}


/**
    @see magma_tally2blas_clacpy_q
    @ingroup magma_tally2_caux2
    ********************************************************************/
extern "C" void
magma_tally2blas_clacpy(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr       dB, magma_tally2_int_t lddb )
{
    magma_tally2blas_clacpy_q( uplo, m, n, dA, ldda, dB, lddb, magma_tally2_stream );
}


////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    CLACPY_BATCHED copies all or part of each two-dimensional matrix
    dAarray[i] to matrix dBarray[i], for 0 <= i < batchcount.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_tally2_uplo_t
            Specifies the part of each matrix dA to be copied to dB.
      -     = Magma_tally2Upper:      Upper triangular part
      -     = Magma_tally2Lower:      Lower triangular part
            Otherwise:  All of each matrix dA
    
    @param[in]
    m       INTEGER
            The number of rows of each matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of each matrix dA.  N >= 0.
    
    @param[in]
    dAarray COMPLEX* array, dimension (batchCount)
            Array of pointers to the matrices dA, where each dA is of dimension (LDDA,N).
            The M-by-N matrix dA.
            If UPLO = Magma_tally2Upper, only the upper triangle or trapezoid is accessed;
            if UPLO = Magma_tally2Lower, only the lower triangle or trapezoid is accessed.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of each array dA.  LDDA >= max(1,M).
    
    @param[out]
    dBarray COMPLEX* array, dimension (batchCount)
            Array of pointers to the matrices dB, where each dB is of dimension (LDDB,N).
            The M-by-N matrix dB.
            On exit, dB = dA in the locations specified by UPLO.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of each array dB.  LDDB >= max(1,M).
    
    @param[in]
    batchCount  Number of matrices in dAarray and dBarray.
    
    @param[in]
    queue   magma_tally2_queue_t
            Queue to execute in.

    @ingroup magma_tally2_caux2
    ********************************************************************/
extern "C" void
magma_tally2blas_clacpy_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr const dAarray[], magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr             dBarray[], magma_tally2_int_t lddb,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    if ( uplo != Magma_tally2Lower && uplo != Magma_tally2Upper && uplo != Magma_tally2Full )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < max(1,m))
        info = -5;
    else if ( lddb < max(1,m))
        info = -7;
    else if ( batchCount < 0 )
        info = -8;
    
    if ( info != 0 ) {
        magma_tally2_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 || n == 0 || batchCount == 0 ) {
        return;
    }
    
    dim3 threads( BLK_X, 1, 1 );
    dim3 grid( magma_tally2_ceildiv( m, BLK_X ), magma_tally2_ceildiv( n, BLK_Y ), batchCount );
    
    if ( uplo == Magma_tally2Lower ) {
        clacpy_lower_kernel_batched<<< grid, threads, 0, queue >>> ( m, n, dAarray, ldda, dBarray, lddb );
    }
    else if ( uplo == Magma_tally2Upper ) {
        clacpy_upper_kernel_batched<<< grid, threads, 0, queue >>> ( m, n, dAarray, ldda, dBarray, lddb );
    }
    else {
        clacpy_full_kernel_batched <<< grid, threads, 0, queue >>> ( m, n, dAarray, ldda, dBarray, lddb );
    }
}
