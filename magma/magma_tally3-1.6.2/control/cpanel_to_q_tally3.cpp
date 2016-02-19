/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from zpanel_to_q_tally3.cpp normal z -> c, Fri Jan 30 19:00:20 2015
*/
#include "common_magma_tally3.h"

// -------------------------
// Put 0s in the upper triangular part of a panel and 1s on the diagonal.
// Stores previous values in work array, to be restored later with cq_to_panel_tally3.
extern "C"
void cpanel_to_q_tally3(magma_tally3_uplo_t uplo, magma_tally3_int_t ib, magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3FloatComplex *work)
{
    magma_tally3_int_t i, j, k = 0;
    magma_tally3FloatComplex *col;
    magma_tally3FloatComplex c_zero = MAGMA_tally3_C_ZERO;
    magma_tally3FloatComplex c_one  = MAGMA_tally3_C_ONE;
    
    if (uplo == Magma_tally3Upper) {
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
// Restores a panel, after call to cpanel_to_q_tally3.
extern "C"
void cq_to_panel_tally3(magma_tally3_uplo_t uplo, magma_tally3_int_t ib, magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3FloatComplex *work)
{
    magma_tally3_int_t i, j, k = 0;
    magma_tally3FloatComplex *col;
    
    if (uplo == Magma_tally3Upper) {
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
