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

    Passes a vector to MAGMA_tally4 (located on DEV).

    Arguments
    ---------

    @param[in]
    m           magma_tally4_int_t
                number of rows

    @param[in]
    n           magma_tally4_int_t
                number of columns

    @param[in]
    val         magma_tally4DoubleComplex_ptr
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
magma_tally4_zvset_dev(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr val,
    magma_tally4_z_matrix *v,
    magma_tally4_queue_t queue )
{
    v->num_rows = m;
    v->num_cols = n;
    v->nnz = m*n;
    v->memory_location = Magma_tally4_DEV;
    v->storage_type = Magma_tally4_DENSE;
    v->dval = val;
    v->major = Magma_tally4ColMajor;

    return MAGMA_tally4_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally4 vector back (located on DEV).

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
    val         magma_tally4DoubleComplex_ptr
                array containing vector entries

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C"
magma_tally4_int_t
magma_tally4_zvget_dev(
    magma_tally4_z_matrix v,
    magma_tally4_int_t *m, magma_tally4_int_t *n,
    magma_tally4DoubleComplex_ptr *val,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info =0;
    
    magma_tally4_z_matrix v_DEV={Magma_tally4_CSR};
    
    if ( v.memory_location == Magma_tally4_DEV ) {

        *m = v.num_rows;
        *n = v.num_cols;
        *val = v.dval;
    } else {
        CHECK( magma_tally4_zmtransfer( v, &v_DEV, v.memory_location, Magma_tally4_DEV, queue ));
        CHECK( magma_tally4_zvget_dev( v_DEV, m, n, val, queue ));
    }
    
cleanup:
    magma_tally4_zmfree( &v_DEV, queue );
    return info;
}


