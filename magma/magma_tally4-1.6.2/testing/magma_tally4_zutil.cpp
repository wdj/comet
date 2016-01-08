/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s

       @author Mark Gates

       Utilities for testing.
*/

#include "testings.h"

#define A(i,j)  A[i + j*lda]

// --------------------
// Make a matrix symmetric/Hermitian.
// Makes diagonal real.
// Sets Aji = conj( Aij ) for j < i, that is, copy lower triangle to upper triangle.
extern "C"
void magma_tally4_zmake_hermitian( magma_tally4_int_t N, magma_tally4DoubleComplex* A, magma_tally4_int_t lda )
{
    magma_tally4_int_t i, j;
    for( i=0; i<N; ++i ) {
        A(i,i) = MAGMA_tally4_Z_MAKE( MAGMA_tally4_Z_REAL( A(i,i) ), 0. );
        for( j=0; j<i; ++j ) {
            A(j,i) = MAGMA_tally4_Z_CNJG( A(i,j) );
        }
    }
}


// --------------------
// Make a matrix symmetric/Hermitian positive definite.
// Increases diagonal by N, and makes it real.
// Sets Aji = conj( Aij ) for j < i, that is, copy lower triangle to upper triangle.
extern "C"
void magma_tally4_zmake_hpd( magma_tally4_int_t N, magma_tally4DoubleComplex* A, magma_tally4_int_t lda )
{
    magma_tally4_int_t i, j;
    for( i=0; i<N; ++i ) {
        A(i,i) = MAGMA_tally4_Z_MAKE( MAGMA_tally4_Z_REAL( A(i,i) ) + N, 0. );
        for( j=0; j<i; ++j ) {
            A(j,i) = MAGMA_tally4_Z_CNJG( A(i,j) );
        }
    }
}
