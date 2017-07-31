/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zcustomspmv.cpp normal z -> c, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    This is an interface to any custom sparse matrix vector product.
    It should compute y = alpha*FUNCTION(x) + beta*y
    The vectors are located on the device, the scalars on the CPU.


    Arguments
    ---------

    @param[in]
    alpha       magma_tally2FloatComplex
                scalar alpha

    @param[in]
    x           magma_tally2_c_matrix
                input vector x
                
    @param[in]
    beta        magma_tally2FloatComplex
                scalar beta
    @param[out]
    y           magma_tally2_c_matrix
                output vector y
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_cblas
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_ccustomspmv(
    magma_tally2FloatComplex alpha,
    magma_tally2_c_matrix x,
    magma_tally2FloatComplex beta,
    magma_tally2_c_matrix y,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    // vector access via x.dval, y.dval
    // sizes are x.num_rows, x.num_cols
    

    return info;

}





