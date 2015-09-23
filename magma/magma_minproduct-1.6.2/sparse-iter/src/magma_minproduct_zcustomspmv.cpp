/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    This is an interface to any custom sparse matrix vector product.
    It should compute y = alpha*FUNCTION(x) + beta*y
    The vectors are located on the device, the scalars on the CPU.


    Arguments
    ---------

    @param[in]
    alpha       magma_minproductDoubleComplex
                scalar alpha

    @param[in]
    x           magma_minproduct_z_matrix
                input vector x
                
    @param[in]
    beta        magma_minproductDoubleComplex
                scalar beta
    @param[out]
    y           magma_minproduct_z_matrix
                output vector y
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zblas
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zcustomspmv(
    magma_minproductDoubleComplex alpha,
    magma_minproduct_z_matrix x,
    magma_minproductDoubleComplex beta,
    magma_minproduct_z_matrix y,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    // vector access via x.dval, y.dval
    // sizes are x.num_rows, x.num_cols
    

    return info;

}





