/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zmlumerge.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

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
    L           magma_tally3_s_matrix
                input strictly lower triangular matrix L

    @param[in]
    U           magma_tally3_s_matrix
                input upper triangular matrix U
    
    @param[out]
    A           magma_tally3_s_matrix*
                output matrix
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_smlumerge(
    magma_tally3_s_matrix L,
    magma_tally3_s_matrix U,
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;    

    if( L.storage_type == Magma_tally3_CSR && U.storage_type == Magma_tally3_CSR ){
        if( L.memory_location == Magma_tally3_CPU && U.memory_location == Magma_tally3_CPU ){
            
            CHECK( magma_tally3_smtransfer( L, A, Magma_tally3_CPU, Magma_tally3_CPU, queue ));
            magma_tally3_free_cpu( A->col );
            magma_tally3_free_cpu( A->val );
            // make sure it is strictly lower triangular
            magma_tally3_int_t z = 0;
            for(magma_tally3_int_t i=0; i<A->num_rows; i++){
                for(magma_tally3_int_t j=L.row[i]; j<L.row[i+1]; j++){
                    if( L.col[j] < i ){// make sure it is strictly lower triangular
                        z++;
                    }
                }
                for(magma_tally3_int_t j=U.row[i]; j<U.row[i+1]; j++){
                    z++;
                }
            }
            A->nnz = z;
            // fill A with the new structure;
            CHECK( magma_tally3_index_malloc_cpu( &A->col, A->nnz ));
            CHECK( magma_tally3_smalloc_cpu( &A->val, A->nnz ));
            z = 0;
            for(magma_tally3_int_t i=0; i<A->num_rows; i++){
                A->row[i] = z;
                for(magma_tally3_int_t j=L.row[i]; j<L.row[i+1]; j++){
                    if( L.col[j] < i ){// make sure it is strictly lower triangular
                        A->col[z] = L.col[j];
                        A->val[z] = L.val[j];
                        z++;
                    }
                }
                for(magma_tally3_int_t j=U.row[i]; j<U.row[i+1]; j++){
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
            info = MAGMA_tally3_ERR_NOT_SUPPORTED;
        }
    }
    else{
            printf("error: matrix in wrong format.\n"); 
            info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
cleanup:
    if( info != 0 ){
        magma_tally3_smfree( A, queue );
    }
    return info;
}





