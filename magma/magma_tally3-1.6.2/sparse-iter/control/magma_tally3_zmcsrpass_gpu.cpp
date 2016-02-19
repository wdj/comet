/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Passes a CSR matrix to MAGMA_tally3 (located on DEV).

    Arguments
    ---------

    @param[in]
    m           magma_tally3_int_t
                number of rows

    @param[in]
    n           magma_tally3_int_t
                number of columns

    @param[in]
    row         magma_tally3Index_ptr
                row pointer

    @param[in]
    col         magma_tally3Index_ptr
                column indices

    @param[in]
    val         magma_tally3DoubleComplex_ptr
                array containing matrix entries

    @param[out]
    A           magma_tally3_z_matrix*
                matrix in magma_tally3 sparse matrix format
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zcsrset_gpu(
    magma_tally3_int_t m,
    magma_tally3_int_t n,
    magma_tally3Index_ptr row,
    magma_tally3Index_ptr col,
    magma_tally3DoubleComplex_ptr val,
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue )
{   
    A->num_rows = m;
    A->num_cols = n;
    magma_tally3_index_t nnz;
    magma_tally3_index_getvector( 1, row+m, 1, &nnz, 1 );
    A->nnz = (magma_tally3_int_t) nnz;
    A->storage_type = Magma_tally3_CSR;
    A->memory_location = Magma_tally3_DEV;
    A->dval = val;
    A->dcol = col;
    A->drow = row;

    return MAGMA_tally3_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally3 matrix to CSR structure (located on DEV).

    Arguments
    ---------

    @param[in]
    A           magma_tally3_z_matrix
                magma_tally3 sparse matrix in CSR format

    @param[out]
    m           magma_tally3_int_t
                number of rows

    @param[out]
    n           magma_tally3_int_t
                number of columns

    @param[out]
    row         magma_tally3Index_ptr
                row pointer

    @param[out]
    col         magma_tally3Index_ptr
                column indices

    @param[out]
    val         magma_tally3DoubleComplex_ptr
                array containing matrix entries

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zcsrget_gpu(
    magma_tally3_z_matrix A,
    magma_tally3_int_t *m,
    magma_tally3_int_t *n,
    magma_tally3Index_ptr *row,
    magma_tally3Index_ptr *col,
    magma_tally3DoubleComplex_ptr *val,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_z_matrix A_DEV={Magma_tally3_CSR}, A_CSR={Magma_tally3_CSR};
    
    if ( A.memory_location == Magma_tally3_DEV && A.storage_type == Magma_tally3_CSR ) {
        *m = A.num_rows;
        *n = A.num_cols;
        *val = A.dval;
        *col = A.dcol;
        *row = A.drow;
    } else {
        CHECK( magma_tally3_zmconvert( A, &A_CSR, A.storage_type, Magma_tally3_CSR, queue ));
        CHECK( magma_tally3_zmtransfer( A_CSR, &A_DEV, A.memory_location, Magma_tally3_DEV, queue ));
        magma_tally3_zcsrget_gpu( A_DEV, m, n, row, col, val, queue );
    }
    
cleanup:
    magma_tally3_zmfree( &A_CSR, queue );
    magma_tally3_zmfree( &A_DEV, queue );
    return info;
}


