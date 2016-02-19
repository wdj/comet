/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zcustomprecond.cpp normal z -> s, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

#define PRECISION_s



/**
    Purpose
    -------

    This is an interface to the left solve for any custom preconditioner.
    It should compute x = FUNCTION(b)
    The vectors are located on the device.

    Arguments
    ---------

    @param[in]
    b           magma_tally3_s_matrix
                RHS

    @param[in,out]
    x           magma_tally3_s_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally3_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_sapplycustomprecond_l(
    magma_tally3_s_matrix b,
    magma_tally3_s_matrix *x,
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    // vector access via x.dval, y->dval
    
    return info;
}


/**
    Purpose
    -------

    This is an interface to the right solve for any custom preconditioner.
    It should compute x = FUNCTION(b)
    The vectors are located on the device.

    Arguments
    ---------

    @param[in]
    b           magma_tally3_s_matrix
                RHS

    @param[in,out]
    x           magma_tally3_s_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally3_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_sapplycustomprecond_r(
    magma_tally3_s_matrix b,
    magma_tally3_s_matrix *x,
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    // vector access via x.dval, y->dval
    // sizes are x.num_rows, x.num_cols
    
    
    return info;
}





