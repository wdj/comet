/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zmscale.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Scales a matrix.

    Arguments
    ---------

    @param[in,out]
    A           magma_tally2_c_matrix*
                input/output matrix

    @param[in]
    scaling     magma_tally2_scale_t
                scaling type (unit rownorm / unit diagonal)

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_caux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_cmscale(
    magma_tally2_c_matrix *A,
    magma_tally2_scale_t scaling,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2FloatComplex *tmp=NULL;
    
    magma_tally2_c_matrix hA={Magma_tally2_CSR}, CSRA={Magma_tally2_CSR};

    if ( A->memory_location == Magma_tally2_CPU && A->storage_type == Magma_tally2_CSRCOO ) {
        if ( scaling == Magma_tally2_NOSCALE ) {
            // no scale
            ;
        }
        else if ( scaling == Magma_tally2_UNITROW ) {
            // scale to unit rownorm
            CHECK( magma_tally2_cmalloc_cpu( &tmp, A->num_rows ));
            for( magma_tally2_int_t z=0; z<A->num_rows; z++ ) {
                magma_tally2FloatComplex s = MAGMA_tally2_C_MAKE( 0.0, 0.0 );
                for( magma_tally2_int_t f=A->row[z]; f<A->row[z+1]; f++ )
                    s+= MAGMA_tally2_C_REAL(A->val[f])*MAGMA_tally2_C_REAL(A->val[f]);
                tmp[z] = MAGMA_tally2_C_MAKE( 1.0/sqrt(  MAGMA_tally2_C_REAL( s )  ), 0.0 );
            }        printf("inhere1\n");
            for( magma_tally2_int_t z=0; z<A->nnz; z++ ) {
                A->val[z] = A->val[z] * tmp[A->col[z]] * tmp[A->rowidx[z]];
            }
        }
        else if (scaling == Magma_tally2_UNITDIAG ) {
            // scale to unit diagonal
            CHECK( magma_tally2_cmalloc_cpu( &tmp, A->num_rows ));
            for( magma_tally2_int_t z=0; z<A->num_rows; z++ ) {
                magma_tally2FloatComplex s = MAGMA_tally2_C_MAKE( 0.0, 0.0 );
                for( magma_tally2_int_t f=A->row[z]; f<A->row[z+1]; f++ ) {
                    if ( A->col[f]== z ) {
                        // add some identity matrix
                        //A->val[f] = A->val[f] +  MAGMA_tally2_C_MAKE( 100000.0, 0.0 );
                        s = A->val[f];
                    }
                }
                if ( s == MAGMA_tally2_C_MAKE( 0.0, 0.0 ) ){
                    printf("error: zero diagonal element.\n");
                    info = MAGMA_tally2_ERR;
                }
                tmp[z] = MAGMA_tally2_C_MAKE( 1.0/sqrt(  MAGMA_tally2_C_REAL( s )  ), 0.0 );
                   
            }
            for( magma_tally2_int_t z=0; z<A->nnz; z++ ) {
                A->val[z] = A->val[z] * tmp[A->col[z]] * tmp[A->rowidx[z]];
            }
        }
        else{
            printf( "error: scaling not supported.\n" );
            info = MAGMA_tally2_ERR_NOT_SUPPORTED;
        }
    }
    else {
        magma_tally2_storage_t A_storage = A->storage_type;
        magma_tally2_location_t A_location = A->memory_location;
        CHECK( magma_tally2_cmtransfer( *A, &hA, A->memory_location, Magma_tally2_CPU, queue ));
        CHECK( magma_tally2_cmconvert( hA, &CSRA, hA.storage_type, Magma_tally2_CSRCOO, queue ));

        CHECK( magma_tally2_cmscale( &CSRA, scaling, queue ));

        magma_tally2_cmfree( &hA, queue );
        magma_tally2_cmfree( A, queue );
        CHECK( magma_tally2_cmconvert( CSRA, &hA, Magma_tally2_CSRCOO, A_storage, queue ));
        CHECK( magma_tally2_cmtransfer( hA, A, Magma_tally2_CPU, A_location, queue ));
    }
    
cleanup:
    magma_tally2_free_cpu( tmp );
    magma_tally2_cmfree( &hA, queue );
    magma_tally2_cmfree( &CSRA, queue );
    return info;
}


/**
    Purpose
    -------

    Adds a multiple of the Identity matrix to a matrix: A = A+add * I

    Arguments
    ---------

    @param[in,out]
    A           magma_tally2_c_matrix*
                input/output matrix

    @param[in]
    add         magma_tally2FloatComplex
                scaling for the identity matrix
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_caux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_cmdiagadd(
    magma_tally2_c_matrix *A,
    magma_tally2FloatComplex add,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_c_matrix hA={Magma_tally2_CSR}, CSRA={Magma_tally2_CSR};
    
    if ( A->memory_location == Magma_tally2_CPU && A->storage_type == Magma_tally2_CSRCOO ) {
        for( magma_tally2_int_t z=0; z<A->nnz; z++ ) {
            if ( A->col[z]== A->rowidx[z] ) {
                // add some identity matrix
                A->val[z] = A->val[z] +  add;
            }
        }
    }
    else {
        magma_tally2_storage_t A_storage = A->storage_type;
        magma_tally2_location_t A_location = A->memory_location;
        CHECK( magma_tally2_cmtransfer( *A, &hA, A->memory_location, Magma_tally2_CPU, queue ));
        CHECK( magma_tally2_cmconvert( hA, &CSRA, hA.storage_type, Magma_tally2_CSRCOO, queue ));

        CHECK( magma_tally2_cmdiagadd( &CSRA, add, queue ));

        magma_tally2_cmfree( &hA, queue );
        magma_tally2_cmfree( A, queue );
        CHECK( magma_tally2_cmconvert( CSRA, &hA, Magma_tally2_CSRCOO, A_storage, queue ));
        CHECK( magma_tally2_cmtransfer( hA, A, Magma_tally2_CPU, A_location, queue ));
    }
    
cleanup:
    magma_tally2_cmfree( &hA, queue );
    magma_tally2_cmfree( &CSRA, queue );
    return info;
}



