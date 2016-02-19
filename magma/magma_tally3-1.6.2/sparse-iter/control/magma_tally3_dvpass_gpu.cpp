/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zvpass_gpu.cpp normal z -> d, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Passes a vector to MAGMA_tally3 (located on DEV).

    Arguments
    ---------

    @param[in]
    m           magma_tally3_int_t
                number of rows

    @param[in]
    n           magma_tally3_int_t
                number of columns

    @param[in]
    val         magma_tally3Double_ptr
                array containing vector entries

    @param[out]
    v           magma_tally3_d_matrix*
                magma_tally3 vector
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_daux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_dvset_dev(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr val,
    magma_tally3_d_matrix *v,
    magma_tally3_queue_t queue )
{
    v->num_rows = m;
    v->num_cols = n;
    v->nnz = m*n;
    v->memory_location = Magma_tally3_DEV;
    v->storage_type = Magma_tally3_DENSE;
    v->dval = val;
    v->major = Magma_tally3ColMajor;

    return MAGMA_tally3_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally3 vector back (located on DEV).

    Arguments
    ---------

    @param[in]
    v           magma_tally3_d_matrix
                magma_tally3 vector

    @param[out]
    m           magma_tally3_int_t
                number of rows

    @param[out]
    n           magma_tally3_int_t
                number of columns

    @param[out]
    val         magma_tally3Double_ptr
                array containing vector entries

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_daux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_dvget_dev(
    magma_tally3_d_matrix v,
    magma_tally3_int_t *m, magma_tally3_int_t *n,
    magma_tally3Double_ptr *val,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info =0;
    
    magma_tally3_d_matrix v_DEV={Magma_tally3_CSR};
    
    if ( v.memory_location == Magma_tally3_DEV ) {

        *m = v.num_rows;
        *n = v.num_cols;
        *val = v.dval;
    } else {
        CHECK( magma_tally3_dmtransfer( v, &v_DEV, v.memory_location, Magma_tally3_DEV, queue ));
        CHECK( magma_tally3_dvget_dev( v_DEV, m, n, val, queue ));
    }
    
cleanup:
    magma_tally3_dmfree( &v_DEV, queue );
    return info;
}


