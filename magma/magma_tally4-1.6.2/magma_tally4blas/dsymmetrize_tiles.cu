/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zsymmetrize_tiles.cu normal z -> d, Fri Jan 30 19:00:09 2015
       @author Mark Gates
*/
#include "common_magma_tally4.h"

#define NB 64

/*
    Symmetrizes ntile tiles at a time, e.g., all diagonal tiles of a matrix.
    Grid is ntile x ceil(m/NB).
    Each tile is m x m, and is divided into block rows, each NB x m.
    Each block has NB threads.
    Each thread copies one row, iterating across all columns below diagonal.
    The bottom block of rows may be partially outside the matrix;
    if so, rows outside the matrix (i >= m) are disabled.
*/
__global__ void
dsymmetrize_tiles_lower( int m, double *dA, int ldda, int mstride, int nstride )
{
    // shift dA to tile's top-left corner
    dA += blockIdx.x*(mstride + nstride*ldda);
    
    // dA iterates across row i and dAT iterates down column i.
    int i = blockIdx.y*NB + threadIdx.x;
    double *dAT = dA;
    if ( i < m ) {
        dA  += i;
        dAT += i*ldda;
        double *dAend = dA + i*ldda;
        while( dA < dAend ) {
            *dAT = (*dA);  // upper := lower
            dA  += ldda;
            dAT += 1;
        }
    }
}


// only difference with _lower version is direction dA=dAT instead of dAT=dA.
__global__ void
dsymmetrize_tiles_upper( int m, double *dA, int ldda, int mstride, int nstride )
{
    // shift dA to tile's top-left corner
    dA += blockIdx.x*(mstride + nstride*ldda);
    
    // dA iterates across row i and dAT iterates down column i.
    int i = blockIdx.y*NB + threadIdx.x;
    double *dAT = dA;
    if ( i < m ) {
        dA  += i;
        dAT += i*ldda;
        double *dAend = dA + i*ldda;
        while( dA < dAend ) {
            *dA  = (*dAT);  // lower := upper
            dA  += ldda;
            dAT += 1;
        }
    }
}


/**
    Purpose
    -------
    
    DSYMMETRIZE_TILES copies lower triangle to upper triangle, or vice-versa,
    to make some blocks of dA into general representations of a symmetric block.
    This processes NTILE blocks, typically the diagonal blocks.
    Each block is offset by mstride rows and nstride columns from the previous block.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_tally4_uplo_t
            Specifies the part of the matrix dA that is valid on input.
      -     = Magma_tally4Upper:      Upper triangular part
      -     = Magma_tally4Lower:      Lower triangular part
    
    @param[in]
    m       INTEGER
            The number of rows & columns of each square block of dA.  M >= 0.
    
    @param[in,out]
    dA      DOUBLE_PRECISION array, dimension (LDDA,N)
            The matrix dA. N = m + nstride*(ntile-1).
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1, m + mstride*(ntile-1)).
    
    @param[in]
    ntile   INTEGER
            Number of blocks to symmetrize. ntile >= 0.
    
    @param[in]
    mstride INTEGER
            Row offset from start of one block to start of next block. mstride >= 0.
            Either (mstride >= m) or (nstride >= m), to prevent m-by-m tiles
            from overlapping.
    
    @param[in]
    nstride INTEGER
            Column offset from start of one block to start of next block. nstride >= 0.
    
    @param[in]
    queue   magma_tally4_queue_t
            Queue to execute in.

    @ingroup magma_tally4_daux2
    ********************************************************************/
extern "C" void
magma_tally4blas_dsymmetrize_tiles_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t ntile, magma_tally4_int_t mstride, magma_tally4_int_t nstride,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    if ( uplo != Magma_tally4Lower && uplo != Magma_tally4Upper )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( ldda < max(1,m + mstride*(ntile-1)) )
        info = -5;
    else if ( ntile < 0 )
        info = -6;
    else if ( mstride < 0 )
        info = -7;
    else if ( nstride < 0 )
        info = -8;
    else if ( mstride < m && nstride < m )  // only one must be >= m.
        info = -7;
    
    if ( info != 0 ) {
        magma_tally4_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 || ntile == 0 )
        return;
    
    dim3 threads( NB );
    dim3 grid( ntile, (m + NB - 1)/NB );
    
    //printf( "m %d, grid %d x %d, threads %d\n", m, grid.x, grid.y, threads.x );
    if ( uplo == Magma_tally4Upper ) {
        dsymmetrize_tiles_upper<<< grid, threads, 0, queue >>>( m, dA, ldda, mstride, nstride );
    }
    else {
        dsymmetrize_tiles_lower<<< grid, threads, 0, queue >>>( m, dA, ldda, mstride, nstride );
    }
}


/**
    @see magma_tally4blas_dsymmetrize_tiles_q
    @ingroup magma_tally4_daux2
    ********************************************************************/
extern "C" void
magma_tally4blas_dsymmetrize_tiles(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t ntile, magma_tally4_int_t mstride, magma_tally4_int_t nstride )
{
    magma_tally4blas_dsymmetrize_tiles_q( uplo, m, dA, ldda, ntile, mstride, nstride, magma_tally4_stream );
}
