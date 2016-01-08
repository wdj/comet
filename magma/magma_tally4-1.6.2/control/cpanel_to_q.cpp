/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from zpanel_to_q.cpp normal z -> c, Fri Jan 30 19:00:20 2015
*/
#include "common_magma_tally4.h"

// -------------------------
// Put 0s in the upper triangular part of a panel and 1s on the diagonal.
// Stores previous values in work array, to be restored later with cq_to_panel.
extern "C"
void cpanel_to_q(magma_tally4_uplo_t uplo, magma_tally4_int_t ib, magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4FloatComplex *work)
{
    magma_tally4_int_t i, j, k = 0;
    magma_tally4FloatComplex *col;
    magma_tally4FloatComplex c_zero = MAGMA_tally4_C_ZERO;
    magma_tally4FloatComplex c_one  = MAGMA_tally4_C_ONE;
    
    if (uplo == Magma_tally4Upper) {
        for(i = 0; i < ib; ++i) {
            col = A + i*lda;
            for(j = 0; j < i; ++j) {
                work[k] = col[j];
                col [j] = c_zero;
                ++k;
            }
            
            work[k] = col[i];
            col [j] = c_one;
            ++k;
        }
    }
    else {
        for(i=0; i<ib; ++i) {
            col = A + i*lda;
            work[k] = col[i];
            col [i] = c_one;
            ++k;
            for(j=i+1; j<ib; ++j) {
                work[k] = col[j];
                col [j] = c_zero;
                ++k;
            }
        }
    }
}


// -------------------------
// Restores a panel, after call to cpanel_to_q.
extern "C"
void cq_to_panel(magma_tally4_uplo_t uplo, magma_tally4_int_t ib, magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4FloatComplex *work)
{
    magma_tally4_int_t i, j, k = 0;
    magma_tally4FloatComplex *col;
    
    if (uplo == Magma_tally4Upper) {
        for(i = 0; i < ib; ++i) {
            col = A + i*lda;
            for(j = 0; j <= i; ++j) {
                col[j] = work[k];
                ++k;
            }
        }
    }
    else {
        for(i = 0; i < ib; ++i) {
            col = A + i*lda;
            for(j = i; j < ib; ++j) {
                col[j] = work[k];
                ++k;
            }
        }
    }
}
