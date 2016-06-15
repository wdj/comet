/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zsymmetrize.cu normal z -> c, Fri Jan 30 19:00:09 2015
       @author Mark Gates
*/
#include "common_magma_tally3.h"

#define NB 64

/*
    Matrix is m x m, and is divided into block rows, each NB x m.
    Each block has NB threads.
    Each thread copies one row, iterating across all columns below diagonal.
    The bottom block of rows may be partially outside the matrix;
    if so, rows outside the matrix (i >= m) are disabled.
*/
__global__ void
csymmetrize_lower( int m, magma_tally3FloatComplex *dA, int ldda )
{
    // dA iterates across row i and dAT iterates down column i.
    int i = blockIdx.x*NB + threadIdx.x;
    magma_tally3FloatComplex *dAT = dA;
    if ( i < m ) {
        dA  += i;
        dAT += i*ldda;
        magma_tally3FloatComplex *dAend = dA + i*ldda;
        while( dA < dAend ) {
            *dAT = cuConjf(*dA);  // upper := lower
            dA  += ldda;
            dAT += 1;
        }
    }
}


// only difference with _lower version is direction dA=dAT instead of dAT=dA.
__global__ void
csymmetrize_upper( int m, magma_tally3FloatComplex *dA, int ldda )
{
    // dA iterates across row i and dAT iterates down column i.
    int i = blockIdx.x*NB + threadIdx.x;
    magma_tally3FloatComplex *dAT = dA;
    if ( i < m ) {
        dA  += i;
        dAT += i*ldda;
        magma_tally3FloatComplex *dAend = dA + i*ldda;
        while( dA < dAend ) {
            *dA = cuConjf(*dAT);  // lower := upper
            dA  += ldda;
            dAT += 1;
        }
    }
}


/**
    Purpose
    -------
    
    CSYMMETRIZE copies lower triangle to upper triangle, or vice-versa,
    to make dA a general representation of a symmetric matrix.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_tally3_uplo_t
            Specifies the part of the matrix dA that is valid on input.
      -     = Magma_tally3Upper:      Upper triangular part
      -     = Magma_tally3Lower:      Lower triangular part
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in,out]
    dA      COMPLEX array, dimension (LDDA,N)
            The m by m matrix dA.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[in]
    queue   magma_tally3_queue_t
            Queue to execute in.

    @ingroup magma_tally3_caux2
    ********************************************************************/
extern "C" void
magma_tally3blas_csymmetrize_q(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    if ( uplo != Magma_tally3Lower && uplo != Magma_tally3Upper )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( ldda < max(1,m) )
        info = -4;
    
    if ( info != 0 ) {
        magma_tally3_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 )
        return;
    
    
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    
    if ( uplo == Magma_tally3Upper ) {
        csymmetrize_upper<<< grid, threads, 0, queue >>>( m, dA, ldda );
    }
    else {
        csymmetrize_lower<<< grid, threads, 0, queue >>>( m, dA, ldda );
    }
}


/**
    @see magma_tally3blas_csymmetrize_q
    @ingroup magma_tally3_caux2
    ********************************************************************/
extern "C" void
magma_tally3blas_csymmetrize(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda )
{
    magma_tally3blas_csymmetrize_q( uplo, m, dA, ldda, magma_tally3_stream );
}