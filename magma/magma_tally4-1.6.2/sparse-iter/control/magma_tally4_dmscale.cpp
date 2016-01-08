/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zmscale.cpp normal z -> d, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally4sparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Scales a matrix.

    Arguments
    ---------

    @param[in,out]
    A           magma_tally4_d_matrix*
                input/output matrix

    @param[in]
    scaling     magma_tally4_scale_t
                scaling type (unit rownorm / unit diagonal)

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_daux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dmscale(
    magma_tally4_d_matrix *A,
    magma_tally4_scale_t scaling,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    double *tmp=NULL;
    
    magma_tally4_d_matrix hA={Magma_tally4_CSR}, CSRA={Magma_tally4_CSR};

    if ( A->memory_location == Magma_tally4_CPU && A->storage_type == Magma_tally4_CSRCOO ) {
        if ( scaling == Magma_tally4_NOSCALE ) {
            // no scale
            ;
        }
        else if ( scaling == Magma_tally4_UNITROW ) {
            // scale to unit rownorm
            CHECK( magma_tally4_dmalloc_cpu( &tmp, A->num_rows ));
            for( magma_tally4_int_t z=0; z<A->num_rows; z++ ) {
                double s = MAGMA_tally4_D_MAKE( 0.0, 0.0 );
                for( magma_tally4_int_t f=A->row[z]; f<A->row[z+1]; f++ )
                    s+= MAGMA_tally4_D_REAL(A->val[f])*MAGMA_tally4_D_REAL(A->val[f]);
                tmp[z] = MAGMA_tally4_D_MAKE( 1.0/sqrt(  MAGMA_tally4_D_REAL( s )  ), 0.0 );
            }        printf("inhere1\n");
            for( magma_tally4_int_t z=0; z<A->nnz; z++ ) {
                A->val[z] = A->val[z] * tmp[A->col[z]] * tmp[A->rowidx[z]];
            }
        }
        else if (scaling == Magma_tally4_UNITDIAG ) {
            // scale to unit diagonal
            CHECK( magma_tally4_dmalloc_cpu( &tmp, A->num_rows ));
            for( magma_tally4_int_t z=0; z<A->num_rows; z++ ) {
                double s = MAGMA_tally4_D_MAKE( 0.0, 0.0 );
                for( magma_tally4_int_t f=A->row[z]; f<A->row[z+1]; f++ ) {
                    if ( A->col[f]== z ) {
                        // add some identity matrix
                        //A->val[f] = A->val[f] +  MAGMA_tally4_D_MAKE( 100000.0, 0.0 );
                        s = A->val[f];
                    }
                }
                if ( s == MAGMA_tally4_D_MAKE( 0.0, 0.0 ) ){
                    printf("error: zero diagonal element.\n");
                    info = MAGMA_tally4_ERR;
                }
                tmp[z] = MAGMA_tally4_D_MAKE( 1.0/sqrt(  MAGMA_tally4_D_REAL( s )  ), 0.0 );
                   
            }
            for( magma_tally4_int_t z=0; z<A->nnz; z++ ) {
                A->val[z] = A->val[z] * tmp[A->col[z]] * tmp[A->rowidx[z]];
            }
        }
        else{
            printf( "error: scaling not supported.\n" );
            info = MAGMA_tally4_ERR_NOT_SUPPORTED;
        }
    }
    else {
        magma_tally4_storage_t A_storage = A->storage_type;
        magma_tally4_location_t A_location = A->memory_location;
        CHECK( magma_tally4_dmtransfer( *A, &hA, A->memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_dmconvert( hA, &CSRA, hA.storage_type, Magma_tally4_CSRCOO, queue ));

        CHECK( magma_tally4_dmscale( &CSRA, scaling, queue ));

        magma_tally4_dmfree( &hA, queue );
        magma_tally4_dmfree( A, queue );
        CHECK( magma_tally4_dmconvert( CSRA, &hA, Magma_tally4_CSRCOO, A_storage, queue ));
        CHECK( magma_tally4_dmtransfer( hA, A, Magma_tally4_CPU, A_location, queue ));
    }
    
cleanup:
    magma_tally4_free_cpu( tmp );
    magma_tally4_dmfree( &hA, queue );
    magma_tally4_dmfree( &CSRA, queue );
    return info;
}


/**
    Purpose
    -------

    Adds a multiple of the Identity matrix to a matrix: A = A+add * I

    Arguments
    ---------

    @param[in,out]
    A           magma_tally4_d_matrix*
                input/output matrix

    @param[in]
    add         double
                scaling for the identity matrix
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_daux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dmdiagadd(
    magma_tally4_d_matrix *A,
    double add,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_d_matrix hA={Magma_tally4_CSR}, CSRA={Magma_tally4_CSR};
    
    if ( A->memory_location == Magma_tally4_CPU && A->storage_type == Magma_tally4_CSRCOO ) {
        for( magma_tally4_int_t z=0; z<A->nnz; z++ ) {
            if ( A->col[z]== A->rowidx[z] ) {
                // add some identity matrix
                A->val[z] = A->val[z] +  add;
            }
        }
    }
    else {
        magma_tally4_storage_t A_storage = A->storage_type;
        magma_tally4_location_t A_location = A->memory_location;
        CHECK( magma_tally4_dmtransfer( *A, &hA, A->memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_dmconvert( hA, &CSRA, hA.storage_type, Magma_tally4_CSRCOO, queue ));

        CHECK( magma_tally4_dmdiagadd( &CSRA, add, queue ));

        magma_tally4_dmfree( &hA, queue );
        magma_tally4_dmfree( A, queue );
        CHECK( magma_tally4_dmconvert( CSRA, &hA, Magma_tally4_CSRCOO, A_storage, queue ));
        CHECK( magma_tally4_dmtransfer( hA, A, Magma_tally4_CPU, A_location, queue ));
    }
    
cleanup:
    magma_tally4_dmfree( &hA, queue );
    magma_tally4_dmfree( &CSRA, queue );
    return info;
}



