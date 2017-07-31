/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from zpanel_to_q_tally2.cpp normal z -> d, Fri Jan 30 19:00:20 2015
*/
#include "common_magma_tally2.h"

// -------------------------
// Put 0s in the upper triangular part of a panel and 1s on the diagonal.
// Stores previous values in work array, to be restored later with dq_to_panel_tally2.
extern "C"
void dpanel_to_q_tally2(magma_tally2_uplo_t uplo, magma_tally2_int_t ib, double *A, magma_tally2_int_t lda, double *work)
{
    magma_tally2_int_t i, j, k = 0;
    double *col;
    double c_zero = MAGMA_tally2_D_ZERO;
    double c_one  = MAGMA_tally2_D_ONE;
    
    if (uplo == Magma_tally2Upper) {
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
// Restores a panel, after call to dpanel_to_q_tally2.
extern "C"
void dq_to_panel_tally2(magma_tally2_uplo_t uplo, magma_tally2_int_t ib, double *A, magma_tally2_int_t lda, double *work)
{
    magma_tally2_int_t i, j, k = 0;
    double *col;
    
    if (uplo == Magma_tally2Upper) {
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
