/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zvtranspose.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    Transposes a vector from col to row major and vice versa.


    Arguments
    ---------

    @param[in]
    x           magma_tally4_c_matrix
                input vector

    @param[out]
    y           magma_tally4_c_matrix*
                output vector

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_caux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cvtranspose(
    magma_tally4_c_matrix x,
    magma_tally4_c_matrix *y,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_int_t    m = x.num_rows;
    magma_tally4_int_t    n = x.num_cols;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    magma_tally4_c_matrix x_d={Magma_tally4_CSR}, y_d={Magma_tally4_CSR};
            
    if ( x.memory_location == Magma_tally4_DEV ) {
        CHECK( magma_tally4_cvinit( y, Magma_tally4_DEV, x.num_rows,x.num_cols, MAGMA_tally4_C_ZERO, queue ));
        y->num_rows = x.num_rows;
        y->num_cols = x.num_cols;
        y->storage_type = x.storage_type;
        if ( x.major == Magma_tally4ColMajor) {
            y->major = Magma_tally4RowMajor;
            magma_tally4blas_ctranspose( m, n, x.val, m, y->val, n );
        }
        else {
            y->major = Magma_tally4ColMajor;
            magma_tally4blas_ctranspose( n, m, x.val, n, y->val, m );
        }
    } else {

        CHECK( magma_tally4_cmtransfer( x, &x_d, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
        CHECK( magma_tally4_cvtranspose( x_d, &y_d, queue ));
        CHECK( magma_tally4_cmtransfer( y_d, y, Magma_tally4_DEV, Magma_tally4_CPU, queue ));
    }
    
cleanup:
    if( info != 0 ){
        magma_tally4_cmfree( y, queue );
    }
    magma_tally4_cmfree( &x_d, queue );
    magma_tally4_cmfree( &y_d, queue );
    magma_tally4blasSetKernelStream( orig_queue );
    return info;
}



   


