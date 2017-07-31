/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zsymmetrize.cu normal z -> s, Fri Jan 30 19:00:09 2015
       @author Mark Gates
*/
#include "common_magma_tally2.h"

#define NB 64

/*
    Matrix is m x m, and is divided into block rows, each NB x m.
    Each block has NB threads.
    Each thread copies one row, iterating across all columns below diagonal.
    The bottom block of rows may be partially outside the matrix;
    if so, rows outside the matrix (i >= m) are disabled.
*/
__global__ void
ssymmetrize_lower( int m, float *dA, int ldda )
{
    // dA iterates across row i and dAT iterates down column i.
    int i = blockIdx.x*NB + threadIdx.x;
    float *dAT = dA;
    if ( i < m ) {
        dA  += i;
        dAT += i*ldda;
        float *dAend = dA + i*ldda;
        while( dA < dAend ) {
            *dAT = (*dA);  // upper := lower
            dA  += ldda;
            dAT += 1;
        }
    }
}


// only difference with _lower version is direction dA=dAT instead of dAT=dA.
__global__ void
ssymmetrize_upper( int m, float *dA, int ldda )
{
    // dA iterates across row i and dAT iterates down column i.
    int i = blockIdx.x*NB + threadIdx.x;
    float *dAT = dA;
    if ( i < m ) {
        dA  += i;
        dAT += i*ldda;
        float *dAend = dA + i*ldda;
        while( dA < dAend ) {
            *dA = (*dAT);  // lower := upper
            dA  += ldda;
            dAT += 1;
        }
    }
}


/**
    Purpose
    -------
    
    SSYMMETRIZE copies lower triangle to upper triangle, or vice-versa,
    to make dA a general representation of a symmetric matrix.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_tally2_uplo_t
            Specifies the part of the matrix dA that is valid on input.
      -     = Magma_tally2Upper:      Upper triangular part
      -     = Magma_tally2Lower:      Lower triangular part
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in,out]
    dA      REAL array, dimension (LDDA,N)
            The m by m matrix dA.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[in]
    queue   magma_tally2_queue_t
            Queue to execute in.

    @ingroup magma_tally2_saux2
    ********************************************************************/
extern "C" void
magma_tally2blas_ssymmetrize_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    if ( uplo != Magma_tally2Lower && uplo != Magma_tally2Upper )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( ldda < max(1,m) )
        info = -4;
    
    if ( info != 0 ) {
        magma_tally2_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 )
        return;
    
    
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    
    if ( uplo == Magma_tally2Upper ) {
        ssymmetrize_upper<<< grid, threads, 0, queue >>>( m, dA, ldda );
    }
    else {
        ssymmetrize_lower<<< grid, threads, 0, queue >>>( m, dA, ldda );
    }
}


/**
    @see magma_tally2blas_ssymmetrize_q
    @ingroup magma_tally2_saux2
    ********************************************************************/
extern "C" void
magma_tally2blas_ssymmetrize(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda )
{
    magma_tally2blas_ssymmetrize_q( uplo, m, dA, ldda, magma_tally2_stream );
}
