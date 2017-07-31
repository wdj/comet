/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zfree.cpp normal z -> d, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Free the memory of a magma_tally2_d_matrix.


    Arguments
    ---------

    @param[in,out]
    A           magma_tally2_d_matrix*
                matrix to free
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_daux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_dmfree(
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue )
{
    if ( A->memory_location == Magma_tally2_CPU ) {
       if ( A->storage_type == Magma_tally2_ELL || A->storage_type == Magma_tally2_ELLPACKT ){
            magma_tally2_free_cpu( A->val );
            magma_tally2_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (A->storage_type == Magma_tally2_ELLD ) {
            magma_tally2_free_cpu( A->val );
            magma_tally2_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_ELLRT ) {
            magma_tally2_free_cpu( A->val );
            magma_tally2_free_cpu( A->row );
            magma_tally2_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_SELLP ) {
            magma_tally2_free_cpu( A->val );
            magma_tally2_free_cpu( A->row );
            magma_tally2_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_CSR || A->storage_type == Magma_tally2_CSC
                                        || A->storage_type == Magma_tally2_CSRD
                                        || A->storage_type == Magma_tally2_CSRL
                                        || A->storage_type == Magma_tally2_CSRU ) {
            magma_tally2_free_cpu( A->val );
            magma_tally2_free_cpu( A->col );
            magma_tally2_free_cpu( A->row );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (  A->storage_type == Magma_tally2_CSRCOO ) {
            magma_tally2_free_cpu( A->val );
            magma_tally2_free_cpu( A->col );
            magma_tally2_free_cpu( A->row );
            magma_tally2_free_cpu( A->rowidx );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_BCSR ) {
            magma_tally2_free_cpu( A->val );
            magma_tally2_free_cpu( A->col );
            magma_tally2_free_cpu( A->row );
            magma_tally2_free_cpu( A->blockinfo );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
            A->blockinfo = 0;
        }
        if ( A->storage_type == Magma_tally2_DENSE ) {
            magma_tally2_free_cpu( A->val );
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

    if ( A->memory_location == Magma_tally2_DEV ) {
       if ( A->storage_type == Magma_tally2_ELL || A->storage_type == Magma_tally2_ELLPACKT ){
            if ( magma_tally2_free( A->dval ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->dcol ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_ELLD ) {
            if ( magma_tally2_free( A->dval ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->dcol ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_ELLRT ) {
            if ( magma_tally2_free( A->dval ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->drow ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->dcol ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_SELLP ) {
            if ( magma_tally2_free( A->dval ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->drow ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->dcol ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_CSR || A->storage_type == Magma_tally2_CSC
                                        || A->storage_type == Magma_tally2_CSRD
                                        || A->storage_type == Magma_tally2_CSRL
                                        || A->storage_type == Magma_tally2_CSRU ) {
            if ( magma_tally2_free( A->dval ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->drow ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->dcol ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (  A->storage_type == Magma_tally2_CSRCOO ) {
            if ( magma_tally2_free( A->dval ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->drow ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->dcol ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->drowidx ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_BCSR ) {
            if ( magma_tally2_free( A->dval ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->drow ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            if ( magma_tally2_free( A->dcol ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
            }
            magma_tally2_free_cpu( A->blockinfo );
            A->blockinfo = NULL;
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_tally2_DENSE ) {
            if ( magma_tally2_free( A->dval ) != MAGMA_tally2_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_tally2_ERR_INVALID_PTR; 
                
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
        return MAGMA_tally2_ERR_INVALID_PTR;
    }
    return MAGMA_tally2_SUCCESS;
}



   


