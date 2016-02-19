/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zfree.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Free the memory of a magma_tally3_s_matrix.


    Arguments
    ---------

    @param[in,out]
    A           magma_tally3_s_matrix*
                matrix to free
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_smfree(
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue )
{
    if ( A->memory_location == Magma_tally3_CPU ) {
       if ( A->storage_type == Magma_tally3_ELL || A->storage_type == Magma_tally3_ELLPACKT ){
            magma_tally3_free_cpu( A->val );
            magma_tally3_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (A->storage_type == Magma_tally3_ELLD ) {
            magma_tally3_free_cpu( A->val );
            magma_tally3_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_ELLRT ) {
            magma_tally3_free_cpu( A->val );
            magma_tally3_free_cpu( A->row );
            magma_tally3_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_SELLP ) {
            magma_tally3_free_cpu( A->val );
            magma_tally3_free_cpu( A->row );
            magma_tally3_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_CSR || A->storage_type == Magma_tally3_CSC
                                        || A->storage_type == Magma_tally3_CSRD
                                        || A->storage_type == Magma_tally3_CSRL
                                        || A->storage_type == Magma_tally3_CSRU ) {
            magma_tally3_free_cpu( A->val );
            magma_tally3_free_cpu( A->col );
            magma_tally3_free_cpu( A->row );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (  A->storage_type == Magma_tally3_CSRCOO ) {
            magma_tally3_free_cpu( A->val );
            magma_tally3_free_cpu( A->col );
            magma_tally3_free_cpu( A->row );
            magma_tally3_free_cpu( A->rowidx );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_BCSR ) {
            magma_tally3_free_cpu( A->val );
            magma_tally3_free_cpu( A->col );
            magma_tally3_free_cpu( A->row );
            magma_tally3_free_cpu( A->blockinfo );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
            A->blockinfo = 0;
        }
        if ( A->storage_type == Magma_tally3_DENSE ) {
            magma_tally3_free_cpu( A->val );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        A->val = NULL;
        A->col = NULL;
        A->row = NULL;
        A->rowidx = NULL;
        A->blockinfo = NULL;
        A->diag = NULL;
        A->dval = NULL;
        A->dcol = NULL;
        A->drow = NULL;
        A->drowidx = NULL;
        A->ddiag = NULL;
    }

    if ( A->memory_location == Magma_tally3_DEV ) {
       if ( A->storage_type == Magma_tally3_ELL || A->storage_type == Magma_tally3_ELLPACKT ){
            if ( magma_tally3_free( A->dval ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->dcol ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_ELLD ) {
            if ( magma_tally3_free( A->dval ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->dcol ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_ELLRT ) {
            if ( magma_tally3_free( A->dval ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->drow ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->dcol ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_SELLP ) {
            if ( magma_tally3_free( A->dval ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->drow ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->dcol ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_CSR || A->storage_type == Magma_tally3_CSC
                                        || A->storage_type == Magma_tally3_CSRD
                                        || A->storage_type == Magma_tally3_CSRL
                                        || A->storage_type == Magma_tally3_CSRU ) {
            if ( magma_tally3_free( A->dval ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->drow ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->dcol ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (  A->storage_type == Magma_tally3_CSRCOO ) {
            if ( magma_tally3_free( A->dval ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->drow ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->dcol ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->drowidx ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_BCSR ) {
            if ( magma_tally3_free( A->dval ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->drow ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally3_free( A->dcol ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }
            magma_tally3_free_cpu( A->blockinfo );
            A->blockinfo = NULL;
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally3_DENSE ) {
            if ( magma_tally3_free( A->dval ) != MAGMA_tally3_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally3_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
                
        }
        A->val = NULL;
        A->col = NULL;
        A->row = NULL;
        A->rowidx = NULL;
        A->blockinfo = NULL;
        A->diag = NULL;
        A->dval = NULL;
        A->dcol = NULL;
        A->drow = NULL;
        A->drowidx = NULL;
        A->ddiag = NULL;
    }

    else {
        // printf("Memory Free Error.\n");
        return MAGMA_tally3_ERR_INVALID_PTR;
    }
    return MAGMA_tally3_SUCCESS;
}



   


