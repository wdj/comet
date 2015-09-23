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

    convertes magma_minproduct_c_vector from C to Z

    Arguments
    ---------

    @param[in]
    x           magma_minproduct_c_vector
                input vector descriptor

    @param[out]
    y           magma_minproduct_z_matrix*
                output vector descriptor
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_vector_clag2z(
    magma_minproduct_c_vector x, magma_minproduct_z_matrix *y,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info;
    if ( x.memory_location == Magma_minproduct_DEV) {
        y->memory_location = x.memory_location;
        y->num_rows = x.num_rows;
        y->nnz = x.nnz;
        CHECK( magma_minproduct_zmalloc( &y->val, x.num_rows ));
        magma_minproductblas_clag2z( x.num_rows, 1, x.dval, x.num_rows,
                                    y->val, x.num_rows, &info );
        return MAGMA_minproduct_SUCCESS;
    }
    else if ( x.memory_location == Magma_minproduct_CPU ) {
        y->memory_location = x.memory_location;
        y->num_rows = x.num_rows;
        y->nnz = x.nnz;
        CHECK( magma_minproduct_zmalloc_cpu( &y->val, x.num_rows ));

        magma_minproduct_int_t one= 1;
        magma_minproduct_int_t info;
        lapackf77_clag2z( &x.num_rows, &one,
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

    convertes magma_minproduct_c_sparse_matrix from C to Z

    Arguments
    ---------

    @param
    A           magma_minproduct_c_sparse_matrix
                input matrix descriptor

    @param
    B           magma_minproduct_z_matrix*
                output matrix descriptor
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_sparse_matrix_clag2z(
    magma_minproduct_c_sparse_matrix A, magma_minproduct_z_matrix *B,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    if ( A.memory_location == Magma_minproduct_DEV) {
        B->storage_type = A.storage_type;
        B->memory_location = A.memory_location;
        B->num_rows = A.num_rows;
        B->num_cols = A.num_cols;
        B->nnz = A.nnz;
        B->max_nnz_row = A.max_nnz_row;
        if ( A.storage_type == Magma_minproduct_CSR ) {
            CHECK( magma_minproduct_zmalloc( &B->val, A.nnz ));
            magma_minproductblas_clag2z_sparse( A.nnz, 1, A.dval, A.nnz,
                                            B->val, A.nnz, &info, queue );
            B->row = A.drow;
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_minproduct_ELLPACK ) {
            CHECK( magma_minproduct_zmalloc( &B->val, A.num_rows*A.max_nnz_row ));
            magma_minproductblas_clag2z_sparse( A.num_rows*A.max_nnz_row, 1, A.dval,
            A.num_rows*A.max_nnz_row, B->val, A.num_rows*A.max_nnz_row, &info,
             queue );
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_minproduct_ELL ) {
            CHECK( magma_minproduct_zmalloc( &B->val, A.num_rows*A.max_nnz_row ));
            magma_minproductblas_clag2z_sparse( A.num_rows*A.max_nnz_row, 1, A.dval,
            A.num_rows*A.max_nnz_row, B->val, A.num_rows*A.max_nnz_row, &info,
            queue );
            B->col = A.dcol;
        }
        if ( A.storage_type == Magma_minproduct_DENSE ) {
            CHECK( magma_minproduct_zmalloc( &B->val, A.num_rows*A.num_cols ));
            magma_minproductblas_clag2z_sparse( A.num_rows, A.num_cols, A.dval, A.num_rows,
                    B->val, A.num_rows, &info, queue );
        }
        else
            info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
    else{
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
cleanup:
        return info;
}

