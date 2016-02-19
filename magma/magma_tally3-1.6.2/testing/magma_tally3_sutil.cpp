/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3_zutil.cpp normal z -> s, Fri Jan 30 19:00:26 2015

       @author Mark Gates

       Utilities for testing.
*/

#include "testings.h"

#define A(i,j)  A[i + j*lda]

// --------------------
// Make a matrix symmetric/symmetric.
// Makes diagonal real.
// Sets Aji = conj( Aij ) for j < i, that is, copy lower triangle to upper triangle.
extern "C"
void magma_tally3_smake_symmetric( magma_tally3_int_t N, float* A, magma_tally3_int_t lda )
{
    magma_tally3_int_t i, j;
    for( i=0; i<N; ++i ) {
        A(i,i) = MAGMA_tally3_S_MAKE( MAGMA_tally3_S_REAL( A(i,i) ), 0. );
        for( j=0; j<i; ++j ) {
            A(j,i) = MAGMA_tally3_S_CNJG( A(i,j) );
        }
    }
}


// --------------------
// Make a matrix symmetric/symmetric positive definite.
// Increases diagonal by N, and makes it real.
// Sets Aji = conj( Aij ) for j < i, that is, copy lower triangle to upper triangle.
extern "C"
void magma_tally3_smake_hpd( magma_tally3_int_t N, float* A, magma_tally3_int_t lda )
{
    magma_tally3_int_t i, j;
    for( i=0; i<N; ++i ) {
        A(i,i) = MAGMA_tally3_S_MAKE( MAGMA_tally3_S_REAL( A(i,i) ) + N, 0. );
        for( j=0; j<i; ++j ) {
            A(j,i) = MAGMA_tally3_S_CNJG( A(i,j) );
        }
    }
}
