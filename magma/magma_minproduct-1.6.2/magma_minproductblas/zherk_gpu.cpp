/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Azzam Haidar 
*/
#include "common_magma_minproduct.h"

extern "C"
void magma_minproductblas_zherk_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    double alpha,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    double beta,
    magma_minproductDoubleComplex_ptr dC, magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset)
{
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda + (a_offset) )
    #define dC(i_, j_) (dC + (i_) + (j_)*lddc)
    
    magma_minproduct_trans_t transA;
    magma_minproduct_trans_t transB;  
    magma_minproductDoubleComplex cbeta  = MAGMA_minproduct_Z_MAKE( beta, 0. );
    magma_minproductDoubleComplex calpha = MAGMA_minproduct_Z_MAKE( alpha, 0. );
    
    if (trans == Magma_minproductNoTrans) {
        transA = Magma_minproductNoTrans;
        transB = Magma_minproduct_ConjTrans;
    } else {
        transA = Magma_minproduct_ConjTrans;
        transB = Magma_minproductNoTrans;
    }

    if (uplo == Magma_minproductUpper) {
            printf("Error not supported\n");
            return;
    }

    magma_minproduct_int_t ib, ioff;
    magma_minproduct_int_t blockoffset = c_offset % nb;
    // loop over all blocks and does A * A**H
    // blockoffset is c_offset within first block; for subsequent blocks it is 0
    for( magma_minproduct_int_t i = 0; i < n; i += ib ) {
        ib     = min( nb-blockoffset, n-i );  // block size
        ioff   = i + c_offset;                  // global index in parent matrix
        // C[i:n,i] += A[i:n,0] * A[i,0]'
        // printf( "zgemm  n=%4d, ib=%4d, k=%4d, i=%4d  ioff=%4d\n", n-i, ib, k, i, ioff );
        magma_minproduct_zgemm( transA, transB, n-i, ib, k,
                     calpha, dA(i,0),       ldda,
                             dA(i,0),       ldda,
                     cbeta,  dC(ioff,ioff), lddc );
        blockoffset = 0;
    }
}
