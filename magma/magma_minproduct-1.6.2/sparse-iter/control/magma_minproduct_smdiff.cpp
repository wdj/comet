/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zmdiff.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_minproductsparse.h"

#define THRESHOLD 10e-99


/**
    Purpose
    -------

    Computes the Frobenius norm of the difference between the CSR matrices A
    and B. They do not need to share the same sparsity pattern!
        
            res = ||A-B||_F = sqrt( sum_ij (A_ij-B_ij)^2 )


    Arguments
    ---------

    @param[in]
    A           magma_minproduct_s_matrix
                sparse matrix in CSR

    @param[in]
    B           magma_minproduct_s_matrix
                sparse matrix in CSR
                
    @param[out]
    res         real_Double_t*
                residual
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_smdiff(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix B,
    real_Double_t *res,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    if( A.memory_location == Magma_minproduct_CPU && B.memory_location == Magma_minproduct_CPU
            && A.storage_type == Magma_minproduct_CSR && B.storage_type == Magma_minproduct_CSR ){
        real_Double_t tmp2;
        magma_minproduct_int_t i,j,k;
        *res = 0.0;
        
        for(i=0; i<A.num_rows; i++) {
            for(j=A.row[i]; j<A.row[i+1]; j++) {
                magma_minproduct_index_t localcol = A.col[j];
                for( k=B.row[i]; k<B.row[i+1]; k++) {
                    if (B.col[k] == localcol) {
                        tmp2 = (real_Double_t) fabs( MAGMA_minproduct_S_REAL(A.val[j] )
                                                        - MAGMA_minproduct_S_REAL(B.val[k]) );
    
                        (*res) = (*res) + tmp2* tmp2;
                    }
                }
            }
        }

        (*res) =  sqrt((*res));
    }
    else{
        printf("error: mdiff only supported for CSR matrices on the CPU.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
    return info;
}

