/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @precisions normal z -> s d c

*/
#include "common_magma_minproduct.h"

#define COMPLEX

/**
    Purpose
    -------

    magma_minproduct_zprint prints a matrix that is located on the CPU host.
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
    A       COMPLEX_16 array, dimension (LDA,N), on the CPU host.
            The M-by-N matrix to be printed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @ingroup magma_minproduct_zaux2
    ********************************************************************/
extern "C"
void magma_minproduct_zprint(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductDoubleComplex *A, magma_minproduct_int_t lda )
{
    #define A(i,j) (A + (i) + (j)*lda)
    
    magma_minproduct_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( lda < max(1,m) )
        info = -4;
    
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_minproductDoubleComplex c_zero = MAGMA_minproduct_Z_ZERO;
    
    if ( m == 1 ) {
        printf( "[ " );
    }
    else {
        printf( "[\n" );
    }
    for( int i = 0; i < m; ++i ) {
        for( int j = 0; j < n; ++j ) {
            if ( MAGMA_minproduct_Z_EQUAL( *A(i,j), c_zero )) {
                #ifdef COMPLEX
                printf( "   0.              " );
                #else
                printf( "   0.    " );
                #endif
            }
            else {
                #ifdef COMPLEX
                printf( " %8.4f+%8.4fi", MAGMA_minproduct_Z_REAL( *A(i,j) ), MAGMA_minproduct_Z_IMAG( *A(i,j) ));
                #else
                printf( " %8.4f", MAGMA_minproduct_Z_REAL( *A(i,j) ));
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
    magma_minproduct_zprint_gpu prints a matrix that is located on the GPU device.
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
    dA      COMPLEX_16 array, dimension (LDDA,N), on the GPU device.
            The M-by-N matrix to be printed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @ingroup magma_minproduct_zaux2
    ********************************************************************/
extern "C"
void magma_minproduct_zprint_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductDoubleComplex *dA, magma_minproduct_int_t ldda )
{
    magma_minproduct_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < max(1,m) )
        info = -4;
    
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_minproduct_int_t lda = m;
    magma_minproductDoubleComplex* A;
    magma_minproduct_zmalloc_cpu( &A, lda*n );
    magma_minproduct_zgetmatrix( m, n, dA, ldda, A, lda );
    
    magma_minproduct_zprint( m, n, A, lda );
    
    magma_minproduct_free_cpu( A );
}
