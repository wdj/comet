/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zmcsrpass.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    Passes a CSR matrix to MAGMA_minproduct.

    Arguments
    ---------

    @param[in]
    m           magma_minproduct_int_t
                number of rows

    @param[in]
    n           magma_minproduct_int_t
                number of columns

    @param[in]
    row         magma_minproduct_index_t*
                row pointer

    @param[in]
    col         magma_minproduct_index_t*
                column indices

    @param[in]
    val         magma_minproductFloatComplex*
                array containing matrix entries

    @param[out]
    A           magma_minproduct_c_matrix*
                matrix in magma_minproduct sparse matrix format
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_ccsrset(
    magma_minproduct_int_t m,
    magma_minproduct_int_t n,
    magma_minproduct_index_t *row,
    magma_minproduct_index_t *col,
    magma_minproductFloatComplex *val,
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue )
{
    A->num_rows = m;
    A->num_cols = n;
    A->nnz = row[m];
    A->storage_type = Magma_minproduct_CSR;
    A->memory_location = Magma_minproduct_CPU;
    A->val = val;
    A->col = col;
    A->row = row;
    A->fill_mode = Magma_minproduct_FULL;

    return MAGMA_minproduct_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_minproduct matrix to CSR structure.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_c_matrix
                magma_minproduct sparse matrix in CSR format

    @param[out]
    m           magma_minproduct_int_t
                number of rows

    @param[out]
    n           magma_minproduct_int_t
                number of columns

    @param[out]
    row         magma_minproduct_index_t*
                row pointer

    @param[out]
    col         magma_minproduct_index_t*
                column indices

    @param[out]
    val         magma_minproductFloatComplex*
                array containing matrix entries

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_ccsrget(
    magma_minproduct_c_matrix A,
    magma_minproduct_int_t *m,
    magma_minproduct_int_t *n,
    magma_minproduct_index_t **row,
    magma_minproduct_index_t **col,
    magma_minproductFloatComplex **val,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_c_matrix A_CPU={Magma_minproduct_CSR}, A_CSR={Magma_minproduct_CSR};
        
    if ( A.memory_location == Magma_minproduct_CPU && A.storage_type == Magma_minproduct_CSR ) {
        *m = A.num_rows;
        *n = A.num_cols;
        *val = A.val;
        *col = A.col;
        *row = A.row;
    } else {
        CHECK( magma_minproduct_cmtransfer( A, &A_CPU, A.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_cmconvert( A_CPU, &A_CSR, A_CPU.storage_type, Magma_minproduct_CSR, queue ));
        CHECK( magma_minproduct_ccsrget( A_CSR, m, n, row, col, val, queue ));
    }

cleanup:
    magma_minproduct_cmfree( &A_CSR, queue );
    magma_minproduct_cmfree( &A_CPU, queue );
    return info;
}


