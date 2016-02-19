/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zmscale.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Scales a matrix.

    Arguments
    ---------

    @param[in,out]
    A           magma_tally3_s_matrix*
                input/output matrix

    @param[in]
    scaling     magma_tally3_scale_t
                scaling type (unit rownorm / unit diagonal)

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_smscale(
    magma_tally3_s_matrix *A,
    magma_tally3_scale_t scaling,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    float *tmp=NULL;
    
    magma_tally3_s_matrix hA={Magma_tally3_CSR}, CSRA={Magma_tally3_CSR};

    if ( A->memory_location == Magma_tally3_CPU && A->storage_type == Magma_tally3_CSRCOO ) {
        if ( scaling == Magma_tally3_NOSCALE ) {
            // no scale
            ;
        }
        else if ( scaling == Magma_tally3_UNITROW ) {
            // scale to unit rownorm
            CHECK( magma_tally3_smalloc_cpu( &tmp, A->num_rows ));
            for( magma_tally3_int_t z=0; z<A->num_rows; z++ ) {
                float s = MAGMA_tally3_S_MAKE( 0.0, 0.0 );
                for( magma_tally3_int_t f=A->row[z]; f<A->row[z+1]; f++ )
                    s+= MAGMA_tally3_S_REAL(A->val[f])*MAGMA_tally3_S_REAL(A->val[f]);
                tmp[z] = MAGMA_tally3_S_MAKE( 1.0/sqrt(  MAGMA_tally3_S_REAL( s )  ), 0.0 );
            }        printf("inhere1\n");
            for( magma_tally3_int_t z=0; z<A->nnz; z++ ) {
                A->val[z] = A->val[z] * tmp[A->col[z]] * tmp[A->rowidx[z]];
            }
        }
        else if (scaling == Magma_tally3_UNITDIAG ) {
            // scale to unit diagonal
            CHECK( magma_tally3_smalloc_cpu( &tmp, A->num_rows ));
            for( magma_tally3_int_t z=0; z<A->num_rows; z++ ) {
                float s = MAGMA_tally3_S_MAKE( 0.0, 0.0 );
                for( magma_tally3_int_t f=A->row[z]; f<A->row[z+1]; f++ ) {
                    if ( A->col[f]== z ) {
                        // add some identity matrix
                        //A->val[f] = A->val[f] +  MAGMA_tally3_S_MAKE( 100000.0, 0.0 );
                        s = A->val[f];
                    }
                }
                if ( s == MAGMA_tally3_S_MAKE( 0.0, 0.0 ) ){
                    printf("error: zero diagonal element.\n");
                    info = MAGMA_tally3_ERR;
                }
                tmp[z] = MAGMA_tally3_S_MAKE( 1.0/sqrt(  MAGMA_tally3_S_REAL( s )  ), 0.0 );
                   
            }
            for( magma_tally3_int_t z=0; z<A->nnz; z++ ) {
                A->val[z] = A->val[z] * tmp[A->col[z]] * tmp[A->rowidx[z]];
            }
        }
        else{
            printf( "error: scaling not supported.\n" );
            info = MAGMA_tally3_ERR_NOT_SUPPORTED;
        }
    }
    else {
        magma_tally3_storage_t A_storage = A->storage_type;
        magma_tally3_location_t A_location = A->memory_location;
        CHECK( magma_tally3_smtransfer( *A, &hA, A->memory_location, Magma_tally3_CPU, queue ));
        CHECK( magma_tally3_smconvert( hA, &CSRA, hA.storage_type, Magma_tally3_CSRCOO, queue ));

        CHECK( magma_tally3_smscale( &CSRA, scaling, queue ));

        magma_tally3_smfree( &hA, queue );
        magma_tally3_smfree( A, queue );
        CHECK( magma_tally3_smconvert( CSRA, &hA, Magma_tally3_CSRCOO, A_storage, queue ));
        CHECK( magma_tally3_smtransfer( hA, A, Magma_tally3_CPU, A_location, queue ));
    }
    
cleanup:
    magma_tally3_free_cpu( tmp );
    magma_tally3_smfree( &hA, queue );
    magma_tally3_smfree( &CSRA, queue );
    return info;
}


/**
    Purpose
    -------

    Adds a multiple of the Identity matrix to a matrix: A = A+add * I

    Arguments
    ---------

    @param[in,out]
    A           magma_tally3_s_matrix*
                input/output matrix

    @param[in]
    add         float
                scaling for the identity matrix
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_smdiagadd(
    magma_tally3_s_matrix *A,
    float add,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_s_matrix hA={Magma_tally3_CSR}, CSRA={Magma_tally3_CSR};
    
    if ( A->memory_location == Magma_tally3_CPU && A->storage_type == Magma_tally3_CSRCOO ) {
        for( magma_tally3_int_t z=0; z<A->nnz; z++ ) {
            if ( A->col[z]== A->rowidx[z] ) {
                // add some identity matrix
                A->val[z] = A->val[z] +  add;
            }
        }
    }
    else {
        magma_tally3_storage_t A_storage = A->storage_type;
        magma_tally3_location_t A_location = A->memory_location;
        CHECK( magma_tally3_smtransfer( *A, &hA, A->memory_location, Magma_tally3_CPU, queue ));
        CHECK( magma_tally3_smconvert( hA, &CSRA, hA.storage_type, Magma_tally3_CSRCOO, queue ));

        CHECK( magma_tally3_smdiagadd( &CSRA, add, queue ));

        magma_tally3_smfree( &hA, queue );
        magma_tally3_smfree( A, queue );
        CHECK( magma_tally3_smconvert( CSRA, &hA, Magma_tally3_CSRCOO, A_storage, queue ));
        CHECK( magma_tally3_smtransfer( hA, A, Magma_tally3_CPU, A_location, queue ));
    }
    
cleanup:
    magma_tally3_smfree( &hA, queue );
    magma_tally3_smfree( &CSRA, queue );
    return info;
}



