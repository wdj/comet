/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Azzam Haidar 
*/
#include "common_magma_tally4.h"

extern "C"
void magma_tally4blas_zherk_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    double alpha,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t a_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC, magma_tally4_int_t lddc, magma_tally4_int_t c_offset)
{
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda + (a_offset) )
    #define dC(i_, j_) (dC + (i_) + (j_)*lddc)
    
    magma_tally4_trans_t transA;
    magma_tally4_trans_t transB;  
    magma_tally4DoubleComplex cbeta  = MAGMA_tally4_Z_MAKE( beta, 0. );
    magma_tally4DoubleComplex calpha = MAGMA_tally4_Z_MAKE( alpha, 0. );
    
    if (trans == Magma_tally4NoTrans) {
        transA = Magma_tally4NoTrans;
        transB = Magma_tally4_ConjTrans;
    } else {
        transA = Magma_tally4_ConjTrans;
        transB = Magma_tally4NoTrans;
    }

    if (uplo == Magma_tally4Upper) {
            printf("Error not supported\n");
            return;
    }

    magma_tally4_int_t ib, ioff;
    magma_tally4_int_t blockoffset = c_offset % nb;
    // loop over all blocks and does A * A**H
    // blockoffset is c_offset within first block; for subsequent blocks it is 0
    for( magma_tally4_int_t i = 0; i < n; i += ib ) {
        ib     = min( nb-blockoffset, n-i );  // block size
        ioff   = i + c_offset;                  // global index in parent matrix
        // C[i:n,i] += A[i:n,0] * A[i,0]'
        // printf( "zgemm  n=%4d, ib=%4d, k=%4d, i=%4d  ioff=%4d\n", n-i, ib, k, i, ioff );
        magma_tally4_zgemm( transA, transB, n-i, ib, k,
                     calpha, dA(i,0),       ldda,
                             dA(i,0),       ldda,
                     cbeta,  dC(ioff,ioff), lddc );
        blockoffset = 0;
    }
}
