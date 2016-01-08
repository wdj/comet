/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zmcsrcompressor.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    Removes zeros in a CSR matrix.

    Arguments
    ---------

    @param[in,out]
    A           magma_tally4_s_matrix*
                input/output matrix
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_saux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_smcsrcompressor(
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;

    magma_tally4_s_matrix B={Magma_tally4_CSR};
    magma_tally4_s_matrix hA={Magma_tally4_CSR}, CSRA={Magma_tally4_CSR};
        
    if ( A->memory_location == Magma_tally4_CPU && A->storage_type == Magma_tally4_CSR ) {


        CHECK( magma_tally4_smconvert( *A, &B, Magma_tally4_CSR, Magma_tally4_CSR, queue ));

        magma_tally4_free_cpu( A->row );
        magma_tally4_free_cpu( A->col );
        magma_tally4_free_cpu( A->val );
        CHECK( magma_tally4_s_csr_compressor(&B.val, &B.row, &B.col,
                       &A->val, &A->row, &A->col, &A->num_rows, queue ));
        A->nnz = A->row[A->num_rows];
    }
    else {

        magma_tally4_storage_t A_storage = A->storage_type;
        magma_tally4_location_t A_location = A->memory_location;
        CHECK( magma_tally4_smtransfer( *A, &hA, A->memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_smconvert( hA, &CSRA, hA.storage_type, Magma_tally4_CSR, queue ));

        CHECK( magma_tally4_smcsrcompressor( &CSRA, queue ));

        magma_tally4_smfree( &hA, queue );
        magma_tally4_smfree( A, queue );
        CHECK( magma_tally4_smconvert( CSRA, &hA, Magma_tally4_CSR, A_storage, queue ));
        CHECK( magma_tally4_smtransfer( hA, A, Magma_tally4_CPU, A_location, queue ));
        magma_tally4_smfree( &hA, queue );
        magma_tally4_smfree( &CSRA, queue );
    }
    
cleanup:
    magma_tally4_smfree( &hA, queue );
    magma_tally4_smfree( &CSRA, queue );
    magma_tally4_smfree( &B, queue );
    return info;
}


