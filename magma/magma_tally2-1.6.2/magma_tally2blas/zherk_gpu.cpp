/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Azzam Haidar 
*/
#include "common_magma_tally2.h"

extern "C"
void magma_tally2blas_zherk_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    double alpha,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t a_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC, magma_tally2_int_t lddc, magma_tally2_int_t c_offset)
{
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda + (a_offset) )
    #define dC(i_, j_) (dC + (i_) + (j_)*lddc)
    
    magma_tally2_trans_t transA;
    magma_tally2_trans_t transB;  
    magma_tally2DoubleComplex cbeta  = MAGMA_tally2_Z_MAKE( beta, 0. );
    magma_tally2DoubleComplex calpha = MAGMA_tally2_Z_MAKE( alpha, 0. );
    
    if (trans == Magma_tally2NoTrans) {
        transA = Magma_tally2NoTrans;
        transB = Magma_tally2_ConjTrans;
    } else {
        transA = Magma_tally2_ConjTrans;
        transB = Magma_tally2NoTrans;
    }

    if (uplo == Magma_tally2Upper) {
            printf("Error not supported\n");
            return;
    }

    magma_tally2_int_t ib, ioff;
    magma_tally2_int_t blockoffset = c_offset % nb;
    // loop over all blocks and does A * A**H
    // blockoffset is c_offset within first block; for subsequent blocks it is 0
    for( magma_tally2_int_t i = 0; i < n; i += ib ) {
        ib     = min( nb-blockoffset, n-i );  // block size
        ioff   = i + c_offset;                  // global index in parent matrix
        // C[i:n,i] += A[i:n,0] * A[i,0]'
        // printf( "zgemm  n=%4d, ib=%4d, k=%4d, i=%4d  ioff=%4d\n", n-i, ib, k, i, ioff );
        magma_tally2_zgemm( transA, transB, n-i, ib, k,
                     calpha, dA(i,0),       ldda,
                             dA(i,0),       ldda,
                     cbeta,  dC(ioff,ioff), lddc );
        blockoffset = 0;
    }
}
