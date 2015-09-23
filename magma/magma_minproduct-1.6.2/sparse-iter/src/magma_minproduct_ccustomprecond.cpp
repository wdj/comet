/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zcustomprecond.cpp normal z -> c, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"

#define PRECISION_c



/**
    Purpose
    -------

    This is an interface to the left solve for any custom preconditioner.
    It should compute x = FUNCTION(b)
    The vectors are located on the device.

    Arguments
    ---------

    @param[in]
    b           magma_minproduct_c_matrix
                RHS

    @param[in,out]
    x           magma_minproduct_c_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_minproduct_c_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_capplycustomprecond_l(
    magma_minproduct_c_matrix b,
    magma_minproduct_c_matrix *x,
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
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
    b           magma_minproduct_c_matrix
                RHS

    @param[in,out]
    x           magma_minproduct_c_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_minproduct_c_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_capplycustomprecond_r(
    magma_minproduct_c_matrix b,
    magma_minproduct_c_matrix *x,
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    // vector access via x.dval, y->dval
    // sizes are x.num_rows, x.num_cols
    
    
    return info;
}





