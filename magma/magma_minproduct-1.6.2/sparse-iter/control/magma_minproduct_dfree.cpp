/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zfree.cpp normal z -> d, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    Free the memory of a magma_minproduct_d_matrix.


    Arguments
    ---------

    @param[in,out]
    A           magma_minproduct_d_matrix*
                matrix to free
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_dmfree(
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue )
{
    if ( A->memory_location == Magma_minproduct_CPU ) {
       if ( A->storage_type == Magma_minproduct_ELL || A->storage_type == Magma_minproduct_ELLPACKT ){
            magma_minproduct_free_cpu( A->val );
            magma_minproduct_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (A->storage_type == Magma_minproduct_ELLD ) {
            magma_minproduct_free_cpu( A->val );
            magma_minproduct_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_ELLRT ) {
            magma_minproduct_free_cpu( A->val );
            magma_minproduct_free_cpu( A->row );
            magma_minproduct_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_SELLP ) {
            magma_minproduct_free_cpu( A->val );
            magma_minproduct_free_cpu( A->row );
            magma_minproduct_free_cpu( A->col );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_CSR || A->storage_type == Magma_minproduct_CSC
                                        || A->storage_type == Magma_minproduct_CSRD
                                        || A->storage_type == Magma_minproduct_CSRL
                                        || A->storage_type == Magma_minproduct_CSRU ) {
            magma_minproduct_free_cpu( A->val );
            magma_minproduct_free_cpu( A->col );
            magma_minproduct_free_cpu( A->row );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (  A->storage_type == Magma_minproduct_CSRCOO ) {
            magma_minproduct_free_cpu( A->val );
            magma_minproduct_free_cpu( A->col );
            magma_minproduct_free_cpu( A->row );
            magma_minproduct_free_cpu( A->rowidx );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_BCSR ) {
            magma_minproduct_free_cpu( A->val );
            magma_minproduct_free_cpu( A->col );
            magma_minproduct_free_cpu( A->row );
            magma_minproduct_free_cpu( A->blockinfo );
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
            A->blockinfo = 0;
        }
        if ( A->storage_type == Magma_minproduct_DENSE ) {
            magma_minproduct_free_cpu( A->val );
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

    if ( A->memory_location == Magma_minproduct_DEV ) {
       if ( A->storage_type == Magma_minproduct_ELL || A->storage_type == Magma_minproduct_ELLPACKT ){
            if ( magma_minproduct_free( A->dval ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->dcol ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_ELLD ) {
            if ( magma_minproduct_free( A->dval ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->dcol ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_ELLRT ) {
            if ( magma_minproduct_free( A->dval ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->drow ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->dcol ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_SELLP ) {
            if ( magma_minproduct_free( A->dval ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->drow ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->dcol ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_CSR || A->storage_type == Magma_minproduct_CSC
                                        || A->storage_type == Magma_minproduct_CSRD
                                        || A->storage_type == Magma_minproduct_CSRL
                                        || A->storage_type == Magma_minproduct_CSRU ) {
            if ( magma_minproduct_free( A->dval ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->drow ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->dcol ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if (  A->storage_type == Magma_minproduct_CSRCOO ) {
            if ( magma_minproduct_free( A->dval ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->drow ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->dcol ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->drowidx ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }

            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_BCSR ) {
            if ( magma_minproduct_free( A->dval ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->drow ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            if ( magma_minproduct_free( A->dcol ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
            }
            magma_minproduct_free_cpu( A->blockinfo );
            A->blockinfo = NULL;
            A->num_rows = 0;
            A->num_cols = 0;
            A->nnz = 0;
        }
        if ( A->storage_type == Magma_minproduct_DENSE ) {
            if ( magma_minproduct_free( A->dval ) != MAGMA_minproduct_SUCCESS ) {
                printf("Memory Free Error.\n");
                return MAGMA_minproduct_ERR_INVALID_PTR; 
                
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
        return MAGMA_minproduct_ERR_INVALID_PTR;
    }
    return MAGMA_minproduct_SUCCESS;
}



   


