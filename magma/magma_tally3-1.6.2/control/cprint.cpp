/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from zprint.cpp normal z -> c, Fri Jan 30 19:00:20 2015

*/
#include "common_magma_tally3.h"

#define COMPLEX

/**
    Purpose
    -------

    magma_tally3_cprint prints a matrix that is located on the CPU host.
    The output is intended to be Matlab compatible, to be useful in debugging.

    Arguments
    ---------

    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    A       COMPLEX array, dimension (LDA,N), on the CPU host.
            The M-by-N matrix to be printed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @ingroup magma_tally3_caux2
    ********************************************************************/
extern "C"
void magma_tally3_cprint(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3FloatComplex *A, magma_tally3_int_t lda )
{
    #define A(i,j) (A + (i) + (j)*lda)
    
    magma_tally3_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( lda < max(1,m) )
        info = -4;
    
    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally3FloatComplex c_zero = MAGMA_tally3_C_ZERO;
    
    if ( m == 1 ) {
        printf( "[ " );
    }
    else {
        printf( "[\n" );
    }
    for( int i = 0; i < m; ++i ) {
        for( int j = 0; j < n; ++j ) {
            if ( MAGMA_tally3_C_EQUAL( *A(i,j), c_zero )) {
                #ifdef COMPLEX
                printf( "   0.              " );
                #else
                printf( "   0.    " );
                #endif
            }
            else {
                #ifdef COMPLEX
                printf( " %8.4f+%8.4fi", MAGMA_tally3_C_REAL( *A(i,j) ), MAGMA_tally3_C_IMAG( *A(i,j) ));
                #else
                printf( " %8.4f", MAGMA_tally3_C_REAL( *A(i,j) ));
                #endif
            }
        }
        if ( m > 1 ) {
            printf( "\n" );
        }
        else {
            printf( " " );
        }
    }
    printf( "];\n" );
}


/**
    Purpose
    -------
    magma_tally3_cprint_gpu prints a matrix that is located on the GPU device.
    Internally, it allocates CPU memory and copies the matrix to the CPU.
    The output is intended to be Matlab compatible, to be useful in debugging.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    dA      COMPLEX array, dimension (LDDA,N), on the GPU device.
            The M-by-N matrix to be printed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @ingroup magma_tally3_caux2
    ********************************************************************/
extern "C"
void magma_tally3_cprint_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3FloatComplex *dA, magma_tally3_int_t ldda )
{
    magma_tally3_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < max(1,m) )
        info = -4;
    
    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally3_int_t lda = m;
    magma_tally3FloatComplex* A;
    magma_tally3_cmalloc_cpu( &A, lda*n );
    magma_tally3_cgetmatrix( m, n, dA, ldda, A, lda );
    
    magma_tally3_cprint( m, n, A, lda );
    
    magma_tally3_free_cpu( A );
}
