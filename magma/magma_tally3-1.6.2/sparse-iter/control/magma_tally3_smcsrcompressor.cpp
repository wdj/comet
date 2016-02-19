/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zmcsrcompressor.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Removes zeros in a CSR matrix.

    Arguments
    ---------

    @param[in,out]
    A           magma_tally3_s_matrix*
                input/output matrix
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_smcsrcompressor(
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    magma_tally3_s_matrix B={Magma_tally3_CSR};
    magma_tally3_s_matrix hA={Magma_tally3_CSR}, CSRA={Magma_tally3_CSR};
        
    if ( A->memory_location == Magma_tally3_CPU && A->storage_type == Magma_tally3_CSR ) {


        CHECK( magma_tally3_smconvert( *A, &B, Magma_tally3_CSR, Magma_tally3_CSR, queue ));

        magma_tally3_free_cpu( A->row );
        magma_tally3_free_cpu( A->col );
        magma_tally3_free_cpu( A->val );
        CHECK( magma_tally3_s_csr_compressor(&B.val, &B.row, &B.col,
                       &A->val, &A->row, &A->col, &A->num_rows, queue ));
        A->nnz = A->row[A->num_rows];
    }
    else {

        magma_tally3_storage_t A_storage = A->storage_type;
        magma_tally3_location_t A_location = A->memory_location;
        CHECK( magma_tally3_smtransfer( *A, &hA, A->memory_location, Magma_tally3_CPU, queue ));
        CHECK( magma_tally3_smconvert( hA, &CSRA, hA.storage_type, Magma_tally3_CSR, queue ));

        CHECK( magma_tally3_smcsrcompressor( &CSRA, queue ));

        magma_tally3_smfree( &hA, queue );
        magma_tally3_smfree( A, queue );
        CHECK( magma_tally3_smconvert( CSRA, &hA, Magma_tally3_CSR, A_storage, queue ));
        CHECK( magma_tally3_smtransfer( hA, A, Magma_tally3_CPU, A_location, queue ));
        magma_tally3_smfree( &hA, queue );
        magma_tally3_smfree( &CSRA, queue );
    }
    
cleanup:
    magma_tally3_smfree( &hA, queue );
    magma_tally3_smfree( &CSRA, queue );
    magma_tally3_smfree( &B, queue );
    return info;
}


