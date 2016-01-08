/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zmcsrpass.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    Passes a CSR matrix to MAGMA_tally4.

    Arguments
    ---------

    @param[in]
    m           magma_tally4_int_t
                number of rows

    @param[in]
    n           magma_tally4_int_t
                number of columns

    @param[in]
    row         magma_tally4_index_t*
                row pointer

    @param[in]
    col         magma_tally4_index_t*
                column indices

    @param[in]
    val         float*
                array containing matrix entries

    @param[out]
    A           magma_tally4_s_matrix*
                matrix in magma_tally4 sparse matrix format
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_saux
    ********************************************************************/

extern "C"
magma_tally4_int_t
magma_tally4_scsrset(
    magma_tally4_int_t m,
    magma_tally4_int_t n,
    magma_tally4_index_t *row,
    magma_tally4_index_t *col,
    float *val,
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue )
{
    A->num_rows = m;
    A->num_cols = n;
    A->nnz = row[m];
    A->storage_type = Magma_tally4_CSR;
    A->memory_location = Magma_tally4_CPU;
    A->val = val;
    A->col = col;
    A->row = row;
    A->fill_mode = Magma_tally4_FULL;

    return MAGMA_tally4_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally4 matrix to CSR structure.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_s_matrix
                magma_tally4 sparse matrix in CSR format

    @param[out]
    m           magma_tally4_int_t
                number of rows

    @param[out]
    n           magma_tally4_int_t
                number of columns

    @param[out]
    row         magma_tally4_index_t*
                row pointer

    @param[out]
    col         magma_tally4_index_t*
                column indices

    @param[out]
    val         float*
                array containing matrix entries

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_saux
    ********************************************************************/

extern "C"
magma_tally4_int_t
magma_tally4_scsrget(
    magma_tally4_s_matrix A,
    magma_tally4_int_t *m,
    magma_tally4_int_t *n,
    magma_tally4_index_t **row,
    magma_tally4_index_t **col,
    float **val,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_s_matrix A_CPU={Magma_tally4_CSR}, A_CSR={Magma_tally4_CSR};
        
    if ( A.memory_location == Magma_tally4_CPU && A.storage_type == Magma_tally4_CSR ) {
        *m = A.num_rows;
        *n = A.num_cols;
        *val = A.val;
        *col = A.col;
        *row = A.row;
    } else {
        CHECK( magma_tally4_smtransfer( A, &A_CPU, A.memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_smconvert( A_CPU, &A_CSR, A_CPU.storage_type, Magma_tally4_CSR, queue ));
        CHECK( magma_tally4_scsrget( A_CSR, m, n, row, col, val, queue ));
    }

cleanup:
    magma_tally4_smfree( &A_CSR, queue );
    magma_tally4_smfree( &A_CPU, queue );
    return info;
}


