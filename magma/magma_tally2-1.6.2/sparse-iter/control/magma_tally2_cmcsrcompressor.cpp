/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zmcsrcompressor.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Removes zeros in a CSR matrix.

    Arguments
    ---------

    @param[in,out]
    A           magma_tally2_c_matrix*
                input/output matrix
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_caux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_cmcsrcompressor(
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;

    magma_tally2_c_matrix B={Magma_tally2_CSR};
    magma_tally2_c_matrix hA={Magma_tally2_CSR}, CSRA={Magma_tally2_CSR};
        
    if ( A->memory_location == Magma_tally2_CPU && A->storage_type == Magma_tally2_CSR ) {


        CHECK( magma_tally2_cmconvert( *A, &B, Magma_tally2_CSR, Magma_tally2_CSR, queue ));

        magma_tally2_free_cpu( A->row );
        magma_tally2_free_cpu( A->col );
        magma_tally2_free_cpu( A->val );
        CHECK( magma_tally2_c_csr_compressor(&B.val, &B.row, &B.col,
                       &A->val, &A->row, &A->col, &A->num_rows, queue ));
        A->nnz = A->row[A->num_rows];
    }
    else {

        magma_tally2_storage_t A_storage = A->storage_type;
        magma_tally2_location_t A_location = A->memory_location;
        CHECK( magma_tally2_cmtransfer( *A, &hA, A->memory_location, Magma_tally2_CPU, queue ));
        CHECK( magma_tally2_cmconvert( hA, &CSRA, hA.storage_type, Magma_tally2_CSR, queue ));

        CHECK( magma_tally2_cmcsrcompressor( &CSRA, queue ));

        magma_tally2_cmfree( &hA, queue );
        magma_tally2_cmfree( A, queue );
        CHECK( magma_tally2_cmconvert( CSRA, &hA, Magma_tally2_CSR, A_storage, queue ));
        CHECK( magma_tally2_cmtransfer( hA, A, Magma_tally2_CPU, A_location, queue ));
        magma_tally2_cmfree( &hA, queue );
        magma_tally2_cmfree( &CSRA, queue );
    }
    
cleanup:
    magma_tally2_cmfree( &hA, queue );
    magma_tally2_cmfree( &CSRA, queue );
    magma_tally2_cmfree( &B, queue );
    return info;
}


