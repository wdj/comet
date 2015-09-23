/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Takes an strictly lower triangular matrix L and an upper triangular matrix U
    and merges them into a matrix A containing the upper and lower triangular
    parts.

    Arguments
    ---------

    @param[in]
    L           magma_minproduct_z_matrix
                input strictly lower triangular matrix L

    @param[in]
    U           magma_minproduct_z_matrix
                input upper triangular matrix U
    
    @param[out]
    A           magma_minproduct_z_matrix*
                output matrix
                
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zmlumerge(
    magma_minproduct_z_matrix L,
    magma_minproduct_z_matrix U,
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;    

    if( L.storage_type == Magma_minproduct_CSR && U.storage_type == Magma_minproduct_CSR ){
        if( L.memory_location == Magma_minproduct_CPU && U.memory_location == Magma_minproduct_CPU ){
            
            CHECK( magma_minproduct_zmtransfer( L, A, Magma_minproduct_CPU, Magma_minproduct_CPU, queue ));
            magma_minproduct_free_cpu( A->col );
            magma_minproduct_free_cpu( A->val );
            // make sure it is strictly lower triangular
            magma_minproduct_int_t z = 0;
            for(magma_minproduct_int_t i=0; i<A->num_rows; i++){
                for(magma_minproduct_int_t j=L.row[i]; j<L.row[i+1]; j++){
                    if( L.col[j] < i ){// make sure it is strictly lower triangular
                        z++;
                    }
                }
                for(magma_minproduct_int_t j=U.row[i]; j<U.row[i+1]; j++){
                    z++;
                }
            }
            A->nnz = z;
            // fill A with the new structure;
            CHECK( magma_minproduct_index_malloc_cpu( &A->col, A->nnz ));
            CHECK( magma_minproduct_zmalloc_cpu( &A->val, A->nnz ));
            z = 0;
            for(magma_minproduct_int_t i=0; i<A->num_rows; i++){
                A->row[i] = z;
                for(magma_minproduct_int_t j=L.row[i]; j<L.row[i+1]; j++){
                    if( L.col[j] < i ){// make sure it is strictly lower triangular
                        A->col[z] = L.col[j];
                        A->val[z] = L.val[j];
                        z++;
                    }
                }
                for(magma_minproduct_int_t j=U.row[i]; j<U.row[i+1]; j++){
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
            info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        }
    }
    else{
            printf("error: matrix in wrong format.\n"); 
            info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
cleanup:
    if( info != 0 ){
        magma_minproduct_zmfree( A, queue );
    }
    return info;
}





