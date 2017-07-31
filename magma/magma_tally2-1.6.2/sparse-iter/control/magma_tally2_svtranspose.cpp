/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zvtranspose.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Transposes a vector from col to row major and vice versa.


    Arguments
    ---------

    @param[in]
    x           magma_tally2_s_matrix
                input vector

    @param[out]
    y           magma_tally2_s_matrix*
                output vector

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_svtranspose(
    magma_tally2_s_matrix x,
    magma_tally2_s_matrix *y,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_int_t    m = x.num_rows;
    magma_tally2_int_t    n = x.num_cols;
    
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue=NULL;
    magma_tally2blasGetKernelStream( &orig_queue );

    magma_tally2_s_matrix x_d={Magma_tally2_CSR}, y_d={Magma_tally2_CSR};
            
    if ( x.memory_location == Magma_tally2_DEV ) {
        CHECK( magma_tally2_svinit( y, Magma_tally2_DEV, x.num_rows,x.num_cols, MAGMA_tally2_S_ZERO, queue ));
        y->num_rows = x.num_rows;
        y->num_cols = x.num_cols;
        y->storage_type = x.storage_type;
        if ( x.major == Magma_tally2ColMajor) {
            y->major = Magma_tally2RowMajor;
            magma_tally2blas_stranspose( m, n, x.val, m, y->val, n );
        }
        else {
            y->major = Magma_tally2ColMajor;
            magma_tally2blas_stranspose( n, m, x.val, n, y->val, m );
        }
    } else {

        CHECK( magma_tally2_smtransfer( x, &x_d, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
        CHECK( magma_tally2_svtranspose( x_d, &y_d, queue ));
        CHECK( magma_tally2_smtransfer( y_d, y, Magma_tally2_DEV, Magma_tally2_CPU, queue ));
    }
    
cleanup:
    if( info != 0 ){
        magma_tally2_smfree( y, queue );
    }
    magma_tally2_smfree( &x_d, queue );
    magma_tally2_smfree( &y_d, queue );
    magma_tally2blasSetKernelStream( orig_queue );
    return info;
}



   


