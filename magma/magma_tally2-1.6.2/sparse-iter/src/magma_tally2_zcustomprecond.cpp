/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"

#define PRECISION_z



/**
    Purpose
    -------

    This is an interface to the left solve for any custom preconditioner.
    It should compute x = FUNCTION(b)
    The vectors are located on the device.

    Arguments
    ---------

    @param[in]
    b           magma_tally2_z_matrix
                RHS

    @param[in,out]
    x           magma_tally2_z_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally2_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zaux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zapplycustomprecond_l(
    magma_tally2_z_matrix b,
    magma_tally2_z_matrix *x,
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
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
    b           magma_tally2_z_matrix
                RHS

    @param[in,out]
    x           magma_tally2_z_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally2_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zaux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zapplycustomprecond_r(
    magma_tally2_z_matrix b,
    magma_tally2_z_matrix *x,
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    // vector access via x.dval, y->dval
    // sizes are x.num_rows, x.num_cols
    
    
    return info;
}





