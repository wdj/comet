/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgeadd.cu normal z -> c, Fri Jan 30 19:00:08 2015
       @author Mark Gates
*/
#include "common_magma_tally4.h"

#define BLK_X 64
#define BLK_Y 32

/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.

    Code similar to claset.
*/
__global__
void cgeadd_full(
    int m, int n,
    magma_tally4FloatComplex alpha,
    const magma_tally4FloatComplex *dA, int ldda,
    magma_tally4FloatComplex       *dB, int lddb )
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
                dB[j*lddb] = alpha*dA[j*ldda] + dB[j*lddb];
            }
        }
        else {
            // partial block-column
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                dB[j*lddb] = alpha*dA[j*ldda] + dB[j*lddb];
            }
        }
    }
}


/**
    Purpose
    -------
    ZGEADD adds two matrices, dB = alpha*dA + dB.
    
    Arguments
    ---------
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    alpha   COMPLEX
            The scalar alpha.
            
    @param[in]
    dA      COMPLEX array, dimension (LDDA,N)
            The m by n matrix dA.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
            
    @param[in,out]
    dB      COMPLEX array, dimension (LDDB,N)
            The m by n matrix dB.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  LDDB >= max(1,M).
    
    @param[in]
    queue   magma_tally4_queue_t
            Queue to execute in.

    @ingroup magma_tally4_caux2
    ********************************************************************/
extern "C" void
magma_tally4blas_cgeadd_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < max(1,m))
        info = -5;
    else if ( lddb < max(1,m))
        info = -7;
    
    if ( info != 0 ) {
        magma_tally4_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 || n == 0 )
        return;
    
    dim3 threads( BLK_X );
    dim3 grid( (m + BLK_X - 1)/BLK_X, (n + BLK_Y - 1)/BLK_Y );
    
    cgeadd_full<<< grid, threads, 0, queue >>>
        ( m, n, alpha, dA, ldda, dB, lddb );
}


/**
    @see magma_tally4blas_cgeadd_q
    @ingroup magma_tally4_caux2
    ********************************************************************/
extern "C" void
magma_tally4blas_cgeadd(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dB, magma_tally4_int_t lddb )
{
    magma_tally4blas_cgeadd_q( m, n, alpha, dA, ldda, dB, lddb, magma_tally4_stream );
}
