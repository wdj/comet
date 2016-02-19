/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zvpass.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Passes a vector to MAGMA_tally3.

    Arguments
    ---------

    @param[in]
    m           magma_tally3_int_t
                number of rows

    @param[in]
    n           magma_tally3_int_t
                number of columns

    @param[in]
    val         float*
                array containing vector entries

    @param[out]
    v           magma_tally3_s_matrix*
                magma_tally3 vector
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_svset(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *val,
    magma_tally3_s_matrix *v,
    magma_tally3_queue_t queue )
{
    v->num_rows = m;
    v->num_cols = n;
    v->nnz = m*n;
    v->memory_location = Magma_tally3_CPU;
    v->val = val;
    v->major = Magma_tally3ColMajor;
    v->storage_type = Magma_tally3_DENSE;

    return MAGMA_tally3_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally3 vector back.

    Arguments
    ---------

    @param[in]
    v           magma_tally3_s_matrix
                magma_tally3 vector

    @param[out]
    m           magma_tally3_int_t
                number of rows

    @param[out]
    n           magma_tally3_int_t
                number of columns

    @param[out]
    val         float*
                array containing vector entries

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_svget(
    magma_tally3_s_matrix v,
    magma_tally3_int_t *m, magma_tally3_int_t *n,
    float **val,
    magma_tally3_queue_t queue )
{
    magma_tally3_s_matrix v_CPU={Magma_tally3_CSR};
    magma_tally3_int_t info =0;
    
    if ( v.memory_location == Magma_tally3_CPU ) {

        *m = v.num_rows;
        *n = v.num_cols;
        *val = v.val;
    } else {
        CHECK( magma_tally3_smtransfer( v, &v_CPU, v.memory_location, Magma_tally3_CPU, queue ));
        CHECK( magma_tally3_svget( v_CPU, m, n, val, queue ));
    }
    
cleanup:
    magma_tally3_smfree( &v_CPU, queue );
    return info;
}


