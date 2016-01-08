/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zmlumerge.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally4sparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Takes an strictly lower triangular matrix L and an upper triangular matrix U
    and merges them into a matrix A containing the upper and lower triangular
    parts.

    Arguments
    ---------

    @param[in]
    L           magma_tally4_c_matrix
                input strictly lower triangular matrix L

    @param[in]
    U           magma_tally4_c_matrix
                input upper triangular matrix U
    
    @param[out]
    A           magma_tally4_c_matrix*
                output matrix
                
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_caux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cmlumerge(
    magma_tally4_c_matrix L,
    magma_tally4_c_matrix U,
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;    

    if( L.storage_type == Magma_tally4_CSR && U.storage_type == Magma_tally4_CSR ){
        if( L.memory_location == Magma_tally4_CPU && U.memory_location == Magma_tally4_CPU ){
            
            CHECK( magma_tally4_cmtransfer( L, A, Magma_tally4_CPU, Magma_tally4_CPU, queue ));
            magma_tally4_free_cpu( A->col );
            magma_tally4_free_cpu( A->val );
            // make sure it is strictly lower triangular
            magma_tally4_int_t z = 0;
            for(magma_tally4_int_t i=0; i<A->num_rows; i++){
                for(magma_tally4_int_t j=L.row[i]; j<L.row[i+1]; j++){
                    if( L.col[j] < i ){// make sure it is strictly lower triangular
                        z++;
                    }
                }
                for(magma_tally4_int_t j=U.row[i]; j<U.row[i+1]; j++){
                    z++;
                }
            }
            A->nnz = z;
            // fill A with the new structure;
            CHECK( magma_tally4_index_malloc_cpu( &A->col, A->nnz ));
            CHECK( magma_tally4_cmalloc_cpu( &A->val, A->nnz ));
            z = 0;
            for(magma_tally4_int_t i=0; i<A->num_rows; i++){
                A->row[i] = z;
                for(magma_tally4_int_t j=L.row[i]; j<L.row[i+1]; j++){
                    if( L.col[j] < i ){// make sure it is strictly lower triangular
                        A->col[z] = L.col[j];
                        A->val[z] = L.val[j];
                        z++;
                    }
                }
                for(magma_tally4_int_t j=U.row[i]; j<U.row[i+1]; j++){
                    A->col[z] = U.col[j];
                    A->val[z] = U.val[j];
                    z++;
                }
            }
            A->row[A->num_rows] = z;
            A->nnz = z;
        }
        else{
            printf("error: matrix not on CPU.\n"); 
            info = MAGMA_tally4_ERR_NOT_SUPPORTED;
        }
    }
    else{
            printf("error: matrix in wrong format.\n"); 
            info = MAGMA_tally4_ERR_NOT_SUPPORTED;
    }
cleanup:
    if( info != 0 ){
        magma_tally4_cmfree( A, queue );
    }
    return info;
}





