/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Mark Gates
       
       @precisions normal z -> s d c

*/
#include "common_magma_tally3.h"

#define NB 64

/* ////////////////////////////////////////////////////////////////////////////
 -- GPU kernel for setting the k-1 super-diagonals to OFFDIAG
    and the main diagonal to DIAG.
    Divides matrix into min( ceil((m+k-1)/nb), ceil(n/nb) ) block-columns,
    with k threads in each block.
    Each thread iterates across one diagonal.
    Thread k-1 does the main diagonal, thread k-2 the first super-diagonal, etc.

      block 0           block 1
      0                           => skip above matrix
      1 0                         => skip above matrix
      2 1 0                       => skip above matrix
    [ 3 2 1 0         |         ]
    [   3 2 1 0       |         ]
    [     3 2 1 0     |         ]
    [       3 2 1 0   |         ]
    [         3 2 1 0 |         ]
    [           3 2 1 | 0       ]
    [             3 2 | 1 0     ]
    [               3 | 2 1 0   ]
    [                 | 3 2 1 0 ]
    [                 |   3 2 1 ]
                      |     3 2   => skip below matrix
                              3   => skip below matrix
    
    Thread assignment for m=10, n=12, k=4, nb=8. Each column is done in parallel.
*/
__global__
void zlaset_band_upper(
    int m, int n,
    magma_tally3DoubleComplex offdiag, magma_tally3DoubleComplex diag,
    magma_tally3DoubleComplex *A, int lda)
{
    int k   = blockDim.x;
    int ibx = blockIdx.x * NB;
    int ind = ibx + threadIdx.x - k + 1;
    
    A += ind + ibx*lda;
    
    magma_tally3DoubleComplex value = offdiag;
    if (threadIdx.x == k-1)
        value = diag;

    #pragma unroll
    for (int j=0; j < NB; j++) {
        if (ibx + j < n && ind + j >= 0 && ind + j < m) {
            A[j*(lda+1)] = value;
        }
    }
}

/* ////////////////////////////////////////////////////////////////////////////
 -- GPU kernel for setting the k-1 sub-diagonals to OFFDIAG
    and the main diagonal to DIAG.
    Divides matrix into min( ceil(m/nb), ceil(n/nb) ) block-columns,
    with k threads in each block.
    Each thread iterates across one diagonal.
    Thread 0 does the main diagonal, thread 1 the first sub-diagonal, etc.
    
      block 0           block 1
    [ 0               |         ]
    [ 1 0             |         ]
    [ 2 1 0           |         ]
    [ 3 2 1 0         |         ]
    [   3 2 1 0       |         ]
    [     3 2 1 0     |         ]
    [       3 2 1 0   |         ]
    [         3 2 1 0 |         ]
    [           3 2 1 | 0       ]
    [             3 2 | 1 0     ]
    [               3 | 2 1 0   ]
    [                   3 2 1 0 ]
    [                     3 2 1 ]
                            3 2   => skip below matrix
                              3   => skip below matrix
    
    Thread assignment for m=13, n=12, k=4, nb=8. Each column is done in parallel.
*/

__global__
void zlaset_band_lower(
    int m, int n,
    magma_tally3DoubleComplex offdiag, magma_tally3DoubleComplex diag,
    magma_tally3DoubleComplex *A, int lda)
{
    //int k   = blockDim.x;
    int ibx = blockIdx.x * NB;
    int ind = ibx + threadIdx.x;
    
    A += ind + ibx*lda;
    
    magma_tally3DoubleComplex value = offdiag;
    if (threadIdx.x == 0)
        value = diag;

    #pragma unroll
    for (int j=0; j < NB; j++) {
        if (ibx + j < n && ind + j < m) {
            A[j*(lda+1)] = value;
        }
    }
}


/**
    Purpose
    -------
    ZLASET_BAND_STREAM initializes the main diagonal of dA to DIAG,
    and the K-1 sub- or super-diagonals to OFFDIAG.
    
    This is the same as ZLASET_BAND, but adds queue argument.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_tally3_uplo_t
            Specifies the part of the matrix dA to be set.
      -     = Magma_tally3Upper:      Upper triangular part
      -     = Magma_tally3Lower:      Lower triangular part
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    k       INTEGER
            The number of diagonals to set, including the main diagonal.  K >= 0.
            Currently, K <= 1024 due to CUDA restrictions (max. number of threads per block).
    
    @param[in]
    offdiag COMPLEX_16
            Off-diagonal elements in the band are set to OFFDIAG.
    
    @param[in]
    diag    COMPLEX_16
            All the main diagonal elements are set to DIAG.
    
    @param[in]
    dA      COMPLEX_16 array, dimension (LDDA,N)
            The M-by-N matrix dA.
            If UPLO = Magma_tally3Upper, only the upper triangle or trapezoid is accessed;
            if UPLO = Magma_tally3Lower, only the lower triangle or trapezoid is accessed.
            On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n where i != j, abs(i-j) < k;
                     A(i,i) = BETA , 1 <= i <= min(m,n)
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[in]
    queue   magma_tally3_queue_t
            Stream to execute ZLASET in.
    
    @ingroup magma_tally3_zaux2
    ********************************************************************/
extern "C" void
magma_tally3blas_zlaset_band_q(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex offdiag, magma_tally3DoubleComplex diag,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_queue_t queue)
{
    magma_tally3_int_t info = 0;
    if ( uplo != Magma_tally3Lower && uplo != Magma_tally3Upper )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( k < 0 || k > 1024 )
        info = -4;
    else if ( ldda < max(1,m) )
        info = -6;
    
    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    if (uplo == Magma_tally3Upper) {
        dim3 threads( min(k,n) );
        dim3 grid( (min(m+k-1,n) - 1)/NB + 1 );
        zlaset_band_upper<<< grid, threads, 0, queue >>> (m, n, offdiag, diag, dA, ldda);
}
    else if (uplo == Magma_tally3Lower) {
        dim3 threads( min(k,m) );
        dim3 grid( (min(m,n) - 1)/NB + 1 );
        zlaset_band_lower<<< grid, threads, 0, queue >>> (m, n, offdiag, diag, dA, ldda);
    }
}


/**
    @see magma_tally3blas_zlaset_band_q
    @ingroup magma_tally3_zaux2
    ********************************************************************/
extern "C" void
magma_tally3blas_zlaset_band(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex offdiag, magma_tally3DoubleComplex diag,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda)
{
    magma_tally3blas_zlaset_band_q(uplo, m, n, k, offdiag, diag, dA, ldda, magma_tally3_stream);
}
