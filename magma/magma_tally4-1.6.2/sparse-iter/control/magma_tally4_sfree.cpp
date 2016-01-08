/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zfree.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    Free the memory of a magma_tally4_s_matrix.


    Arguments
    ---------

    @param[in,out]
    A           magma_tally4_s_matrix*
                matrix to free
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_saux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_smfree(
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue )
{
    if ( A->memory_location == Magma_tally4_CPU ) {
       if ( A->storage_type == Magma_tally4_ELL || A->storage_type == Magma_tally4_ELLPACKT ){
            magma_tally4_free_cpu( A->val );
            magma_tally4_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (A->storage_type == Magma_tally4_ELLD ) {
            magma_tally4_free_cpu( A->val );
            magma_tally4_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_ELLRT ) {
            magma_tally4_free_cpu( A->val );
            magma_tally4_free_cpu( A->row );
            magma_tally4_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_SELLP ) {
            magma_tally4_free_cpu( A->val );
            magma_tally4_free_cpu( A->row );
            magma_tally4_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_CSR || A->storage_type == Magma_tally4_CSC
                                        || A->storage_type == Magma_tally4_CSRD
                                        || A->storage_type == Magma_tally4_CSRL
                                        || A->storage_type == Magma_tally4_CSRU ) {
            magma_tally4_free_cpu( A->val );
            magma_tally4_free_cpu( A->col );
            magma_tally4_free_cpu( A->row );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (  A->storage_type == Magma_tally4_CSRCOO ) {
            magma_tally4_free_cpu( A->val );
            magma_tally4_free_cpu( A->col );
            magma_tally4_free_cpu( A->row );
            magma_tally4_free_cpu( A->rowidx );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_BCSR ) {
            magma_tally4_free_cpu( A->val );
            magma_tally4_free_cpu( A->col );
            magma_tally4_free_cpu( A->row );
            magma_tally4_free_cpu( A->blockinfo );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
            A->blockinfo = 0;
        }
        if ( A->storage_type == Magma_tally4_DENSE ) {
            magma_tally4_free_cpu( A->val );
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

    if ( A->memory_location == Magma_tally4_DEV ) {
       if ( A->storage_type == Magma_tally4_ELL || A->storage_type == Magma_tally4_ELLPACKT ){
            if ( magma_tally4_free( A->dval ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->dcol ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_ELLD ) {
            if ( magma_tally4_free( A->dval ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->dcol ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_ELLRT ) {
            if ( magma_tally4_free( A->dval ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->drow ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->dcol ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_SELLP ) {
            if ( magma_tally4_free( A->dval ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->drow ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->dcol ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_CSR || A->storage_type == Magma_tally4_CSC
                                        || A->storage_type == Magma_tally4_CSRD
                                        || A->storage_type == Magma_tally4_CSRL
                                        || A->storage_type == Magma_tally4_CSRU ) {
            if ( magma_tally4_free( A->dval ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->drow ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->dcol ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (  A->storage_type == Magma_tally4_CSRCOO ) {
            if ( magma_tally4_free( A->dval ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->drow ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->dcol ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->drowidx ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_BCSR ) {
            if ( magma_tally4_free( A->dval ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->drow ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally4_free( A->dcol ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
            }
            magma_tally4_free_cpu( A->blockinfo );
            A->blockinfo = NULL;
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally4_DENSE ) {
            if ( magma_tally4_free( A->dval ) != MAGMA_tally4_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally4_ERR_INVALID_PTR; 
                
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
        return MAGMA_tally4_ERR_INVALID_PTR;
    }
    return MAGMA_tally4_SUCCESS;
}



   


