/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zmcsrpass_gpu.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Passes a CSR matrix to MAGMA_tally2 (located on DEV).

    Arguments
    ---------

    @param[in]
    m           magma_tally2_int_t
                number of rows

    @param[in]
    n           magma_tally2_int_t
                number of columns

    @param[in]
    row         magma_tally2Index_ptr
                row pointer

    @param[in]
    col         magma_tally2Index_ptr
                column indices

    @param[in]
    val         magma_tally2Float_ptr
                array containing matrix entries

    @param[out]
    A           magma_tally2_s_matrix*
                matrix in magma_tally2 sparse matrix format
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_scsrset_gpu(
    magma_tally2_int_t m,
    magma_tally2_int_t n,
    magma_tally2Index_ptr row,
    magma_tally2Index_ptr col,
    magma_tally2Float_ptr val,
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue )
{   
    A->num_rows = m;
    A->num_cols = n;
    magma_tally2_index_t nnz;
    magma_tally2_index_getvector( 1, row+m, 1, &nnz, 1 );
    A->nnz = (magma_tally2_int_t) nnz;
    A->storage_type = Magma_tally2_CSR;
    A->memory_location = Magma_tally2_DEV;
    A->dval = val;
    A->dcol = col;
    A->drow = row;

    return MAGMA_tally2_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally2 matrix to CSR structure (located on DEV).

    Arguments
    ---------

    @param[in]
    A           magma_tally2_s_matrix
                magma_tally2 sparse matrix in CSR format

    @param[out]
    m           magma_tally2_int_t
                number of rows

    @param[out]
    n           magma_tally2_int_t
                number of columns

    @param[out]
    row         magma_tally2Index_ptr
                row pointer

    @param[out]
    col         magma_tally2Index_ptr
                column indices

    @param[out]
    val         magma_tally2Float_ptr
                array containing matrix entries

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_scsrget_gpu(
    magma_tally2_s_matrix A,
    magma_tally2_int_t *m,
    magma_tally2_int_t *n,
    magma_tally2Index_ptr *row,
    magma_tally2Index_ptr *col,
    magma_tally2Float_ptr *val,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_s_matrix A_DEV={Magma_tally2_CSR}, A_CSR={Magma_tally2_CSR};
    
    if ( A.memory_location == Magma_tally2_DEV && A.storage_type == Magma_tally2_CSR ) {
        *m = A.num_rows;
        *n = A.num_cols;
        *val = A.dval;
        *col = A.dcol;
        *row = A.drow;
    } else {
        CHECK( magma_tally2_smconvert( A, &A_CSR, A.storage_type, Magma_tally2_CSR, queue ));
        CHECK( magma_tally2_smtransfer( A_CSR, &A_DEV, A.memory_location, Magma_tally2_DEV, queue ));
        magma_tally2_scsrget_gpu( A_DEV, m, n, row, col, val, queue );
    }
    
cleanup:
    magma_tally2_smfree( &A_CSR, queue );
    magma_tally2_smfree( &A_DEV, queue );
    return info;
}


