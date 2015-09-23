/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zmcsrpass_gpu.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    Passes a CSR matrix to MAGMA_minproduct (located on DEV).

    Arguments
    ---------

    @param[in]
    m           magma_minproduct_int_t
                number of rows

    @param[in]
    n           magma_minproduct_int_t
                number of columns

    @param[in]
    row         magma_minproductIndex_ptr
                row pointer

    @param[in]
    col         magma_minproductIndex_ptr
                column indices

    @param[in]
    val         magma_minproductFloat_ptr
                array containing matrix entries

    @param[out]
    A           magma_minproduct_s_matrix*
                matrix in magma_minproduct sparse matrix format
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_scsrset_gpu(
    magma_minproduct_int_t m,
    magma_minproduct_int_t n,
    magma_minproductIndex_ptr row,
    magma_minproductIndex_ptr col,
    magma_minproductFloat_ptr val,
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue )
{   
    A->num_rows = m;
    A->num_cols = n;
    magma_minproduct_index_t nnz;
    magma_minproduct_index_getvector( 1, row+m, 1, &nnz, 1 );
    A->nnz = (magma_minproduct_int_t) nnz;
    A->storage_type = Magma_minproduct_CSR;
    A->memory_location = Magma_minproduct_DEV;
    A->dval = val;
    A->dcol = col;
    A->drow = row;

    return MAGMA_minproduct_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_minproduct matrix to CSR structure (located on DEV).

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_s_matrix
                magma_minproduct sparse matrix in CSR format

    @param[out]
    m           magma_minproduct_int_t
                number of rows

    @param[out]
    n           magma_minproduct_int_t
                number of columns

    @param[out]
    row         magma_minproductIndex_ptr
                row pointer

    @param[out]
    col         magma_minproductIndex_ptr
                column indices

    @param[out]
    val         magma_minproductFloat_ptr
                array containing matrix entries

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_scsrget_gpu(
    magma_minproduct_s_matrix A,
    magma_minproduct_int_t *m,
    magma_minproduct_int_t *n,
    magma_minproductIndex_ptr *row,
    magma_minproductIndex_ptr *col,
    magma_minproductFloat_ptr *val,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_s_matrix A_DEV={Magma_minproduct_CSR}, A_CSR={Magma_minproduct_CSR};
    
    if ( A.memory_location == Magma_minproduct_DEV && A.storage_type == Magma_minproduct_CSR ) {
        *m = A.num_rows;
        *n = A.num_cols;
        *val = A.dval;
        *col = A.dcol;
        *row = A.drow;
    } else {
        CHECK( magma_minproduct_smconvert( A, &A_CSR, A.storage_type, Magma_minproduct_CSR, queue ));
        CHECK( magma_minproduct_smtransfer( A_CSR, &A_DEV, A.memory_location, Magma_minproduct_DEV, queue ));
        magma_minproduct_scsrget_gpu( A_DEV, m, n, row, col, val, queue );
    }
    
cleanup:
    magma_minproduct_smfree( &A_CSR, queue );
    magma_minproduct_smfree( &A_DEV, queue );
    return info;
}


