/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions mixed zc -> ds
       @author Hartwig Anzt
*/
#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    convertes magma_minproduct_z_matrix from Z to C

    Arguments
    ---------

    @param
    x           magma_minproduct_z_matrix
                input vector descriptor

    @param
    y           magma_minproduct_c_vector*
                output vector descriptor
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_vector_zlag2c(
    magma_minproduct_z_matrix x, magma_minproduct_c_vector *y,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info;
    if ( x.memory_location == Magma_minproduct_DEV) {
        y->memory_location = x.memory_location;
        y->num_rows = x.num_rows;
        y->nnz = x.nnz;
        magma_minproduct_cmalloc( &y->val, x.num_rows );
        magma_minproductblas_zlag2c_sparse( x.num_rows, 1, x.dval, x.num_rows, y->val,
                    x.num_rows, &info, queue );
        return MAGMA_minproduct_SUCCESS;
    }
    else if ( x.memory_location == Magma_minproduct_CPU ) {
        y->memory_location = x.memory_location;
        y->num_rows = x.num_rows;
        y->nnz = x.nnz;
        magma_minproduct_cmalloc_cpu( &y->val, x.num_rows );

        magma_minproduct_int_t one= 1;
        magma_minproduct_int_t info;
        lapackf77_zlag2c( &x.num_rows, &one,
                       x.dval, &x.num_rows,
                       y->val, &x.num_rows, &info);
        return MAGMA_minproduct_SUCCESS;

    }
    else
        return MAGMA_minproduct_ERR_NOT_SUPPORTED;
}



/**
    Purpose
    -------

    convertes magma_minproduct_z_matrix from Z to C

    Arguments
    ---------

    @param
    A           magma_minproduct_z_matrix
                input matrix descriptor

    @param
    B           magma_minproduct_c_sparse_matrix*
                output matrix descriptor
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_z
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_sparse_matrix_zlag2c(
    magma_minproduct_z_matrix A, magma_minproduct_c_sparse_matrix *B,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info;
    if ( A.memory_location == Magma_minproduct_DEV) {
        B->storage_type = A.storage_type;
        B->memory_location = A.memory_location;
        B->num_rows = A.num_rows;
        B->num_cols = A.num_cols;
        B->nnz = A.nnz;
        B->max_nnz_row = A.max_nnz_row;
        if ( A.storage_type == Magma_minproduct_CSR ) {
            magma_minproduct_cmalloc( &B->val, A.nnz );
            magma_minproductblas_zlag2c_sparse( A.nnz, 1, A.dval, A.nnz, B->val,
                    A.nnz, &info, queue );
            B->row = A.drow;
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_minproduct_ELLPACK ) {
            magma_minproduct_cmalloc( &B->val, A.num_rows*A.max_nnz_row );
            magma_minproductblas_zlag2c_sparse( A.num_rows*A.max_nnz_row, 1, A.dval,
            A.num_rows*A.max_nnz_row, B->val, A.num_rows*A.max_nnz_row, &info, queue );
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_minproduct_ELL ) {
            magma_minproduct_cmalloc( &B->val, A.num_rows*A.max_nnz_row );
            magma_minproductblas_zlag2c_sparse(  A.num_rows*A.max_nnz_row, 1, A.dval,
            A.num_rows*A.max_nnz_row, B->val, A.num_rows*A.max_nnz_row, &info, queue );
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_minproduct_DENSE ) {
            magma_minproduct_cmalloc( &B->val, A.num_rows*A.num_cols );
            magma_minproductblas_zlag2c_sparse(  A.num_rows, A.num_cols, A.dval, A.num_rows,
                    B->val, A.num_rows, &info, queue );
        }
        else {
            info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        }
    }
    else {
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
    
cleanup:
    return info;
}

