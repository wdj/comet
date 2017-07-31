/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zvpass.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Passes a vector to MAGMA_tally2.

    Arguments
    ---------

    @param[in]
    m           magma_tally2_int_t
                number of rows

    @param[in]
    n           magma_tally2_int_t
                number of columns

    @param[in]
    val         float*
                array containing vector entries

    @param[out]
    v           magma_tally2_s_matrix*
                magma_tally2 vector
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_svset(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *val,
    magma_tally2_s_matrix *v,
    magma_tally2_queue_t queue )
{
    v->num_rows = m;
    v->num_cols = n;
    v->nnz = m*n;
    v->memory_location = Magma_tally2_CPU;
    v->val = val;
    v->major = Magma_tally2ColMajor;
    v->storage_type = Magma_tally2_DENSE;

    return MAGMA_tally2_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_tally2 vector back.

    Arguments
    ---------

    @param[in]
    v           magma_tally2_s_matrix
                magma_tally2 vector

    @param[out]
    m           magma_tally2_int_t
                number of rows

    @param[out]
    n           magma_tally2_int_t
                number of columns

    @param[out]
    val         float*
                array containing vector entries

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_svget(
    magma_tally2_s_matrix v,
    magma_tally2_int_t *m, magma_tally2_int_t *n,
    float **val,
    magma_tally2_queue_t queue )
{
    magma_tally2_s_matrix v_CPU={Magma_tally2_CSR};
    magma_tally2_int_t info =0;
    
    if ( v.memory_location == Magma_tally2_CPU ) {

        *m = v.num_rows;
        *n = v.num_cols;
        *val = v.val;
    } else {
        CHECK( magma_tally2_smtransfer( v, &v_CPU, v.memory_location, Magma_tally2_CPU, queue ));
        CHECK( magma_tally2_svget( v_CPU, m, n, val, queue ));
    }
    
cleanup:
    magma_tally2_smfree( &v_CPU, queue );
    return info;
}


