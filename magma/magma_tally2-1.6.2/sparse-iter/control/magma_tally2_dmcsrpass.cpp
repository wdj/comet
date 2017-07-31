/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zmcsrpass.cpp normal z -> d, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Passes a CSR matrix to MAGMA_tally2.

    Arguments
    ---------

    @param[in]
    m           magma_tally2_int_t
                number of rows

    @param[in]
    n           magma_tally2_int_t
                number of columns

    @param[in]
    row         magma_tally2_index_t*
                row pointer

    @param[in]
    col         magma_tally2_index_t*
                column indices

    @param[in]
    val         double*
                array containing matrix entries

    @param[out]
    A           magma_tally2_d_matrix*
                matrix in magma_tally2 sparse matrix format
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_daux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_dcsrset(
    magma_tally2_int_t m,
    magma_tally2_int_t n,
    magma_tally2_index_t *row,
    magma_tally2_index_t *col,
    double *val,
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue )
{
    A->num_rows = m;
    A->num_cols = n;
    A->nnz = row[m];
    A->storage_type = Magma_tally2_CSR;
    A->memory_location = Magma_tally2_CPU;
    A->val = val;
    A->col = col;
    A->row = row;
    A->fill_mode = Magma_tally2_FULL;

    return MAGMA_tally2_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally2 matrix to CSR structure.

    Arguments
    ---------

    @param[in]
    A           magma_tally2_d_matrix
                magma_tally2 sparse matrix in CSR format

    @param[out]
    m           magma_tally2_int_t
                number of rows

    @param[out]
    n           magma_tally2_int_t
                number of columns

    @param[out]
    row         magma_tally2_index_t*
                row pointer

    @param[out]
    col         magma_tally2_index_t*
                column indices

    @param[out]
    val         double*
                array containing matrix entries

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_daux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_dcsrget(
    magma_tally2_d_matrix A,
    magma_tally2_int_t *m,
    magma_tally2_int_t *n,
    magma_tally2_index_t **row,
    magma_tally2_index_t **col,
    double **val,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_d_matrix A_CPU={Magma_tally2_CSR}, A_CSR={Magma_tally2_CSR};
        
    if ( A.memory_location == Magma_tally2_CPU && A.storage_type == Magma_tally2_CSR ) {
        *m = A.num_rows;
        *n = A.num_cols;
        *val = A.val;
        *col = A.col;
        *row = A.row;
    } else {
        CHECK( magma_tally2_dmtransfer( A, &A_CPU, A.memory_location, Magma_tally2_CPU, queue ));
        CHECK( magma_tally2_dmconvert( A_CPU, &A_CSR, A_CPU.storage_type, Magma_tally2_CSR, queue ));
        CHECK( magma_tally2_dcsrget( A_CSR, m, n, row, col, val, queue ));
    }

cleanup:
    magma_tally2_dmfree( &A_CSR, queue );
    magma_tally2_dmfree( &A_CPU, queue );
    return info;
}


