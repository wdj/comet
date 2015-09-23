/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zmcsrcompressor.cpp normal z -> d, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    Removes zeros in a CSR matrix.

    Arguments
    ---------

    @param[in,out]
    A           magma_minproduct_d_matrix*
                input/output matrix
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_dmcsrcompressor(
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;

    magma_minproduct_d_matrix B={Magma_minproduct_CSR};
    magma_minproduct_d_matrix hA={Magma_minproduct_CSR}, CSRA={Magma_minproduct_CSR};
        
    if ( A->memory_location == Magma_minproduct_CPU && A->storage_type == Magma_minproduct_CSR ) {


        CHECK( magma_minproduct_dmconvert( *A, &B, Magma_minproduct_CSR, Magma_minproduct_CSR, queue ));

        magma_minproduct_free_cpu( A->row );
        magma_minproduct_free_cpu( A->col );
        magma_minproduct_free_cpu( A->val );
        CHECK( magma_minproduct_d_csr_compressor(&B.val, &B.row, &B.col,
                       &A->val, &A->row, &A->col, &A->num_rows, queue ));
        A->nnz = A->row[A->num_rows];
    }
    else {

        magma_minproduct_storage_t A_storage = A->storage_type;
        magma_minproduct_location_t A_location = A->memory_location;
        CHECK( magma_minproduct_dmtransfer( *A, &hA, A->memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_dmconvert( hA, &CSRA, hA.storage_type, Magma_minproduct_CSR, queue ));

        CHECK( magma_minproduct_dmcsrcompressor( &CSRA, queue ));

        magma_minproduct_dmfree( &hA, queue );
        magma_minproduct_dmfree( A, queue );
        CHECK( magma_minproduct_dmconvert( CSRA, &hA, Magma_minproduct_CSR, A_storage, queue ));
        CHECK( magma_minproduct_dmtransfer( hA, A, Magma_minproduct_CPU, A_location, queue ));
        magma_minproduct_dmfree( &hA, queue );
        magma_minproduct_dmfree( &CSRA, queue );
    }
    
cleanup:
    magma_minproduct_dmfree( &hA, queue );
    magma_minproduct_dmfree( &CSRA, queue );
    magma_minproduct_dmfree( &B, queue );
    return info;
}


