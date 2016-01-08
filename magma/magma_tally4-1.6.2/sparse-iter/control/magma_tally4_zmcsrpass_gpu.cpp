/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    Passes a CSR matrix to MAGMA_tally4 (located on DEV).

    Arguments
    ---------

    @param[in]
    m           magma_tally4_int_t
                number of rows

    @param[in]
    n           magma_tally4_int_t
                number of columns

    @param[in]
    row         magma_tally4Index_ptr
                row pointer

    @param[in]
    col         magma_tally4Index_ptr
                column indices

    @param[in]
    val         magma_tally4DoubleComplex_ptr
                array containing matrix entries

    @param[out]
    A           magma_tally4_z_matrix*
                matrix in magma_tally4 sparse matrix format
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C"
magma_tally4_int_t
magma_tally4_zcsrset_gpu(
    magma_tally4_int_t m,
    magma_tally4_int_t n,
    magma_tally4Index_ptr row,
    magma_tally4Index_ptr col,
    magma_tally4DoubleComplex_ptr val,
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue )
{   
    A->num_rows = m;
    A->num_cols = n;
    magma_tally4_index_t nnz;
    magma_tally4_index_getvector( 1, row+m, 1, &nnz, 1 );
    A->nnz = (magma_tally4_int_t) nnz;
    A->storage_type = Magma_tally4_CSR;
    A->memory_location = Magma_tally4_DEV;
    A->dval = val;
    A->dcol = col;
    A->drow = row;

    return MAGMA_tally4_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally4 matrix to CSR structure (located on DEV).

    Arguments
    ---------

    @param[in]
    A           magma_tally4_z_matrix
                magma_tally4 sparse matrix in CSR format

    @param[out]
    m           magma_tally4_int_t
                number of rows

    @param[out]
    n           magma_tally4_int_t
                number of columns

    @param[out]
    row         magma_tally4Index_ptr
                row pointer

    @param[out]
    col         magma_tally4Index_ptr
                column indices

    @param[out]
    val         magma_tally4DoubleComplex_ptr
                array containing matrix entries

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C"
magma_tally4_int_t
magma_tally4_zcsrget_gpu(
    magma_tally4_z_matrix A,
    magma_tally4_int_t *m,
    magma_tally4_int_t *n,
    magma_tally4Index_ptr *row,
    magma_tally4Index_ptr *col,
    magma_tally4DoubleComplex_ptr *val,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_z_matrix A_DEV={Magma_tally4_CSR}, A_CSR={Magma_tally4_CSR};
    
    if ( A.memory_location == Magma_tally4_DEV && A.storage_type == Magma_tally4_CSR ) {
        *m = A.num_rows;
        *n = A.num_cols;
        *val = A.dval;
        *col = A.dcol;
        *row = A.drow;
    } else {
        CHECK( magma_tally4_zmconvert( A, &A_CSR, A.storage_type, Magma_tally4_CSR, queue ));
        CHECK( magma_tally4_zmtransfer( A_CSR, &A_DEV, A.memory_location, Magma_tally4_DEV, queue ));
        magma_tally4_zcsrget_gpu( A_DEV, m, n, row, col, val, queue );
    }
    
cleanup:
    magma_tally4_zmfree( &A_CSR, queue );
    magma_tally4_zmfree( &A_DEV, queue );
    return info;
}


