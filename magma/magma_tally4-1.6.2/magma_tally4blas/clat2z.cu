/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
       @author Mark Gates
*/
#include "common_magma_tally4.h"

#define PRECISION_z

#define BLK_X 64
#define BLK_Y 32


/*
    Divides matrix into ceil( n/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.
    Updates only the diagonal and below.
    Blocks that are fully above the diagonal exit immediately.
    
    Code similar to zlag2c and zlaset.
*/
__global__
void clat2z_lower(
    int n,
    const magma_tally4FloatComplex *SA, int ldsa,
    magma_tally4DoubleComplex      *A,  int lda )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (below diag) */
    bool full = (iby + BLK_Y <= n && (ind >= iby + BLK_Y));
    /* do only rows inside matrix, and blocks not above diag */
    if ( ind < n && ind + BLK_X > iby ) {
        A  += ind + iby*lda;
        SA += ind + iby*ldsa;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                A[j*lda] = cuComplexFloatToDouble( SA[j*ldsa] );
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n && ind >= iby+j; ++j ) {
                A[j*lda] = cuComplexFloatToDouble( SA[j*ldsa] );
            }
        }
    }
}


/*
    Similar to clat2z_full, but updates only the diagonal and above.
    Blocks that are fully below the diagonal exit immediately.
    
    Code similar to zlag2c and zlaset.
*/
__global__
void clat2z_upper(
    int n,
    const magma_tally4FloatComplex *SA, int ldsa,
    magma_tally4DoubleComplex      *A,  int lda )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column && (above diag) */
    bool full = (iby + BLK_Y <= n && (ind + BLK_X <= iby));
    /* do only rows inside matrix, and blocks not below diag */
    if ( ind < n && ind < iby + BLK_Y ) {
        A  += ind + iby*lda;
        SA += ind + iby*ldsa;
        if ( full ) {
            // full block-column, off-diagonal block
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                A[j*lda] = cuComplexFloatToDouble( SA[j*ldsa] );
            }
        }
        else {
            // either partial block-column or diagonal block
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                if ( ind <= iby+j ) {
                    A[j*lda] = cuComplexFloatToDouble( SA[j*ldsa] );
                }
            }
        }
    }
}


/**
    Purpose
    -------
    CLAT2Z_STREAM converts a single-complex matrix, SA,
                        to a double-complex matrix, A.

    Note that while it is possible to overflow while converting
    from double to single, it is not possible to overflow when
    converting from single to double.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
            Specifies the part of the matrix A to be converted.
      -     = Magma_tally4Upper:      Upper triangular part
      -     = Magma_tally4Lower:      Lower triangular part
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  n >= 0.
    
    @param[in]
    A       COMPLEX_16 array, dimension (LDA,n)
            On entry, the n-by-n coefficient matrix A.
    
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,n).
    
    @param[out]
    SA      COMPLEX array, dimension (LDSA,n)
            On exit, if INFO=0, the n-by-n coefficient matrix SA;
            if INFO > 0, the content of SA is unspecified.
    
    @param[in]
    ldsa    INTEGER
            The leading dimension of the array SA.  LDSA >= max(1,n).
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit.
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
    
    @param[in]
    queue   magma_tally4_queue_t
            Queue to execute in.
    
    @ingroup magma_tally4_zaux2
    ********************************************************************/
extern "C" void
magma_tally4blas_clat2z_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr SA, magma_tally4_int_t ldsa,
    magma_tally4DoubleComplex_ptr      A,  magma_tally4_int_t lda,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue )
{
    *info = 0;
    if ( uplo != Magma_tally4Lower && uplo != Magma_tally4Upper )
        *info = -1;
    else if ( n < 0 )
        *info = -2;
    else if ( lda < max(1,n) )
        *info = -4;
    else if ( ldsa < max(1,n) )
        *info = -6;
    
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return; //*info;
    }

    /* quick return */
    if ( n == 0 ) {
        return;
    }
    
    dim3 threads( BLK_X );
    dim3 grid( (n+BLK_X-1)/BLK_X, (n+BLK_Y-1)/BLK_Y );
    
    if (uplo == Magma_tally4Lower)
        clat2z_lower<<< grid, threads, 0, queue >>> (n, SA, ldsa, A, lda);
    else if (uplo == Magma_tally4Upper)                                         
        clat2z_upper<<< grid, threads, 0, queue >>> (n, SA, ldsa, A, lda);
}


/**
    @see magma_tally4blas_clat2z_q
    @ingroup magma_tally4_zaux2
    ********************************************************************/
extern "C" void
magma_tally4blas_clat2z(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr SA, magma_tally4_int_t ldsa,
    magma_tally4DoubleComplex_ptr      A,  magma_tally4_int_t lda,
    magma_tally4_int_t *info )
{
    magma_tally4blas_clat2z_q( uplo, n, SA, ldsa, A, lda, info, magma_tally4_stream );
}
