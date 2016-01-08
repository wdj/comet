/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions mixed zc -> ds
       @author Hartwig Anzt
*/
#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    convertes magma_tally4_c_vector from C to Z

    Arguments
    ---------

    @param[in]
    x           magma_tally4_c_vector
                input vector descriptor

    @param[out]
    y           magma_tally4_z_matrix*
                output vector descriptor
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_caux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_vector_clag2z(
    magma_tally4_c_vector x, magma_tally4_z_matrix *y,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info;
    if ( x.memory_location == Magma_tally4_DEV) {
        y->memory_location = x.memory_location;
        y->num_rows = x.num_rows;
        y->nnz = x.nnz;
        CHECK( magma_tally4_zmalloc( &y->val, x.num_rows ));
        magma_tally4blas_clag2z( x.num_rows, 1, x.dval, x.num_rows,
                                    y->val, x.num_rows, &info );
        return MAGMA_tally4_SUCCESS;
    }
    else if ( x.memory_location == Magma_tally4_CPU ) {
        y->memory_location = x.memory_location;
        y->num_rows = x.num_rows;
        y->nnz = x.nnz;
        CHECK( magma_tally4_zmalloc_cpu( &y->val, x.num_rows ));

        magma_tally4_int_t one= 1;
        magma_tally4_int_t info;
        lapackf77_clag2z( &x.num_rows, &one,
                       x.dval, &x.num_rows,
                       y->val, &x.num_rows, &info);
        return MAGMA_tally4_SUCCESS;

    }
    else
        return MAGMA_tally4_ERR_NOT_SUPPORTED;
}



/**
    Purpose
    -------

    convertes magma_tally4_c_sparse_matrix from C to Z

    Arguments
    ---------

    @param
    A           magma_tally4_c_sparse_matrix
                input matrix descriptor

    @param
    B           magma_tally4_z_matrix*
                output matrix descriptor
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_caux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sparse_matrix_clag2z(
    magma_tally4_c_sparse_matrix A, magma_tally4_z_matrix *B,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    if ( A.memory_location == Magma_tally4_DEV) {
        B->storage_type = A.storage_type;
        B->memory_location = A.memory_location;
        B->num_rows = A.num_rows;
        B->num_cols = A.num_cols;
        B->nnz = A.nnz;
        B->max_nnz_row = A.max_nnz_row;
        if ( A.storage_type == Magma_tally4_CSR ) {
            CHECK( magma_tally4_zmalloc( &B->val, A.nnz ));
            magma_tally4blas_clag2z_sparse( A.nnz, 1, A.dval, A.nnz,
                                            B->val, A.nnz, &info, queue );
            B->row = A.drow;
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_tally4_ELLPACK ) {
            CHECK( magma_tally4_zmalloc( &B->val, A.num_rows*A.max_nnz_row ));
            magma_tally4blas_clag2z_sparse( A.num_rows*A.max_nnz_row, 1, A.dval,
            A.num_rows*A.max_nnz_row, B->val, A.num_rows*A.max_nnz_row, &info,
             queue );
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_tally4_ELL ) {
            CHECK( magma_tally4_zmalloc( &B->val, A.num_rows*A.max_nnz_row ));
            magma_tally4blas_clag2z_sparse( A.num_rows*A.max_nnz_row, 1, A.dval,
            A.num_rows*A.max_nnz_row, B->val, A.num_rows*A.max_nnz_row, &info,
            queue );
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_tally4_DENSE ) {
            CHECK( magma_tally4_zmalloc( &B->val, A.num_rows*A.num_cols ));
            magma_tally4blas_clag2z_sparse( A.num_rows, A.num_cols, A.dval, A.num_rows,
                    B->val, A.num_rows, &info, queue );
        }
        else
            info = MAGMA_tally4_ERR_NOT_SUPPORTED;
    }
    else{
        info = MAGMA_tally4_ERR_NOT_SUPPORTED;
    }
cleanup:
        return info;
}

