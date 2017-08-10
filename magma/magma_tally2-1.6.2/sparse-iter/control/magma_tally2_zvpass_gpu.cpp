/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Passes a vector to MAGMA_tally2 (located on DEV).

    Arguments
    ---------

    @param[in]
    m           magma_tally2_int_t
                number of rows

    @param[in]
    n           magma_tally2_int_t
                number of columns

    @param[in]
    val         magma_tally2DoubleComplex_ptr
                array containing vector entries

    @param[out]
    v           magma_tally2_z_matrix*
                magma_tally2 vector
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zaux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_zvset_dev(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr val,
    magma_tally2_z_matrix *v,
    magma_tally2_queue_t queue )
{
    v->num_rows = m;
    v->num_cols = n;
    v->nnz = m*n;
    v->memory_location = Magma_tally2_DEV;
    v->storage_type = Magma_tally2_DENSE;
    v->dval = val;
    v->major = Magma_tally2ColMajor;

    return MAGMA_tally2_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally2 vector back (located on DEV).

    Arguments
    ---------

    @param[in]
    v           magma_tally2_z_matrix
                magma_tally2 vector

    @param[out]
    m           magma_tally2_int_t
                number of rows

    @param[out]
    n           magma_tally2_int_t
                number of columns

    @param[out]
    val         magma_tally2DoubleComplex_ptr
                array containing vector entries

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zaux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_zvget_dev(
    magma_tally2_z_matrix v,
    magma_tally2_int_t *m, magma_tally2_int_t *n,
    magma_tally2DoubleComplex_ptr *val,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info =0;
    
    magma_tally2_z_matrix v_DEV={Magma_tally2_CSR};
    
    if ( v.memory_location == Magma_tally2_DEV ) {

        *m = v.num_rows;
        *n = v.num_cols;
        *val = v.dval;
    } else {
        CHECK( magma_tally2_zmtransfer( v, &v_DEV, v.memory_location, Magma_tally2_DEV, queue ));
        CHECK( magma_tally2_zvget_dev( v_DEV, m, n, val, queue ));
    }
    
cleanup:
    magma_tally2_zmfree( &v_DEV, queue );
    return info;
}

