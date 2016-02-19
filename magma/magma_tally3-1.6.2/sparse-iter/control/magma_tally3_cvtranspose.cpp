/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zvtranspose.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Transposes a vector from col to row major and vice versa.


    Arguments
    ---------

    @param[in]
    x           magma_tally3_c_matrix
                input vector

    @param[out]
    y           magma_tally3_c_matrix*
                output vector

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_cvtranspose(
    magma_tally3_c_matrix x,
    magma_tally3_c_matrix *y,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_int_t    m = x.num_rows;
    magma_tally3_int_t    n = x.num_cols;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    magma_tally3_c_matrix x_d={Magma_tally3_CSR}, y_d={Magma_tally3_CSR};
            
    if ( x.memory_location == Magma_tally3_DEV ) {
        CHECK( magma_tally3_cvinit( y, Magma_tally3_DEV, x.num_rows,x.num_cols, MAGMA_tally3_C_ZERO, queue ));
        y->num_rows = x.num_rows;
        y->num_cols = x.num_cols;
        y->storage_type = x.storage_type;
        if ( x.major == Magma_tally3ColMajor) {
            y->major = Magma_tally3RowMajor;
            magma_tally3blas_ctranspose( m, n, x.val, m, y->val, n );
        }
        else {
            y->major = Magma_tally3ColMajor;
            magma_tally3blas_ctranspose( n, m, x.val, n, y->val, m );
        }
    } else {

        CHECK( magma_tally3_cmtransfer( x, &x_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
        CHECK( magma_tally3_cvtranspose( x_d, &y_d, queue ));
        CHECK( magma_tally3_cmtransfer( y_d, y, Magma_tally3_DEV, Magma_tally3_CPU, queue ));
    }
    
cleanup:
    if( info != 0 ){
        magma_tally3_cmfree( y, queue );
    }
    magma_tally3_cmfree( &x_d, queue );
    magma_tally3_cmfree( &y_d, queue );
    magma_tally3blasSetKernelStream( orig_queue );
    return info;
}



   


