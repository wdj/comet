/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zvtranspose.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    Transposes a vector from col to row major and vice versa.


    Arguments
    ---------

    @param[in]
    x           magma_minproduct_c_matrix
                input vector

    @param[out]
    y           magma_minproduct_c_matrix*
                output vector

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cvtranspose(
    magma_minproduct_c_matrix x,
    magma_minproduct_c_matrix *y,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_int_t    m = x.num_rows;
    magma_minproduct_int_t    n = x.num_cols;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    magma_minproduct_c_matrix x_d={Magma_minproduct_CSR}, y_d={Magma_minproduct_CSR};
            
    if ( x.memory_location == Magma_minproduct_DEV ) {
        CHECK( magma_minproduct_cvinit( y, Magma_minproduct_DEV, x.num_rows,x.num_cols, MAGMA_minproduct_C_ZERO, queue ));
        y->num_rows = x.num_rows;
        y->num_cols = x.num_cols;
        y->storage_type = x.storage_type;
        if ( x.major == Magma_minproductColMajor) {
            y->major = Magma_minproductRowMajor;
            magma_minproductblas_ctranspose( m, n, x.val, m, y->val, n );
        }
        else {
            y->major = Magma_minproductColMajor;
            magma_minproductblas_ctranspose( n, m, x.val, n, y->val, m );
        }
    } else {

        CHECK( magma_minproduct_cmtransfer( x, &x_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
        CHECK( magma_minproduct_cvtranspose( x_d, &y_d, queue ));
        CHECK( magma_minproduct_cmtransfer( y_d, y, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
    }
    
cleanup:
    if( info != 0 ){
        magma_minproduct_cmfree( y, queue );
    }
    magma_minproduct_cmfree( &x_d, queue );
    magma_minproduct_cmfree( &y_d, queue );
    magma_minproductblasSetKernelStream( orig_queue );
    return info;
}



   


