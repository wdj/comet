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

    Passes a vector to MAGMA_tally4.

    Arguments
    ---------

    @param[in]
    m           magma_tally4_int_t
                number of rows

    @param[in]
    n           magma_tally4_int_t
                number of columns

    @param[in]
    val         magma_tally4DoubleComplex*
                array containing vector entries

    @param[out]
    v           magma_tally4_z_matrix*
                magma_tally4 vector
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C"
magma_tally4_int_t
magma_tally4_zvset(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *val,
    magma_tally4_z_matrix *v,
    magma_tally4_queue_t queue )
{
    v->num_rows = m;
    v->num_cols = n;
    v->nnz = m*n;
    v->memory_location = Magma_tally4_CPU;
    v->val = val;
    v->major = Magma_tally4ColMajor;
    v->storage_type = Magma_tally4_DENSE;

    return MAGMA_tally4_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally4 vector back.

    Arguments
    ---------

    @param[in]
    v           magma_tally4_z_matrix
                magma_tally4 vector

    @param[out]
    m           magma_tally4_int_t
                number of rows

    @param[out]
    n           magma_tally4_int_t
                number of columns

    @param[out]
    val         magma_tally4DoubleComplex*
                array containing vector entries

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C"
magma_tally4_int_t
magma_tally4_zvget(
    magma_tally4_z_matrix v,
    magma_tally4_int_t *m, magma_tally4_int_t *n,
    magma_tally4DoubleComplex **val,
    magma_tally4_queue_t queue )
{
    magma_tally4_z_matrix v_CPU={Magma_tally4_CSR};
    magma_tally4_int_t info =0;
    
    if ( v.memory_location == Magma_tally4_CPU ) {

        *m = v.num_rows;
        *n = v.num_cols;
        *val = v.val;
    } else {
        CHECK( magma_tally4_zmtransfer( v, &v_CPU, v.memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_zvget( v_CPU, m, n, val, queue ));
    }
    
cleanup:
    magma_tally4_zmfree( &v_CPU, queue );
    return info;
}


