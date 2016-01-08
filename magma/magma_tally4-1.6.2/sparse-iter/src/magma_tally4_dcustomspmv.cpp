/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zcustomspmv.cpp normal z -> d, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    This is an interface to any custom sparse matrix vector product.
    It should compute y = alpha*FUNCTION(x) + beta*y
    The vectors are located on the device, the scalars on the CPU.


    Arguments
    ---------

    @param[in]
    alpha       double
                scalar alpha

    @param[in]
    x           magma_tally4_d_matrix
                input vector x
                
    @param[in]
    beta        double
                scalar beta
    @param[out]
    y           magma_tally4_d_matrix
                output vector y
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dblas
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dcustomspmv(
    double alpha,
    magma_tally4_d_matrix x,
    double beta,
    magma_tally4_d_matrix y,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    // vector access via x.dval, y.dval
    // sizes are x.num_rows, x.num_cols
    

    return info;

}





