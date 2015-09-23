/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Scales a matrix.

    Arguments
    ---------

    @param[in,out]
    A           magma_minproduct_z_matrix*
                input/output matrix

    @param[in]
    scaling     magma_minproduct_scale_t
                scaling type (unit rownorm / unit diagonal)

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zmscale(
    magma_minproduct_z_matrix *A,
    magma_minproduct_scale_t scaling,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproductDoubleComplex *tmp=NULL;
    
    magma_minproduct_z_matrix hA={Magma_minproduct_CSR}, CSRA={Magma_minproduct_CSR};

    if ( A->memory_location == Magma_minproduct_CPU && A->storage_type == Magma_minproduct_CSRCOO ) {
        if ( scaling == Magma_minproduct_NOSCALE ) {
            // no scale
            ;
        }
        else if ( scaling == Magma_minproduct_UNITROW ) {
            // scale to unit rownorm
            CHECK( magma_minproduct_zmalloc_cpu( &tmp, A->num_rows ));
            for( magma_minproduct_int_t z=0; z<A->num_rows; z++ ) {
                magma_minproductDoubleComplex s = MAGMA_minproduct_Z_MAKE( 0.0, 0.0 );
                for( magma_minproduct_int_t f=A->row[z]; f<A->row[z+1]; f++ )
                    s+= MAGMA_minproduct_Z_REAL(A->val[f])*MAGMA_minproduct_Z_REAL(A->val[f]);
                tmp[z] = MAGMA_minproduct_Z_MAKE( 1.0/sqrt(  MAGMA_minproduct_Z_REAL( s )  ), 0.0 );
            }        printf("inhere1\n");
            for( magma_minproduct_int_t z=0; z<A->nnz; z++ ) {
                A->val[z] = A->val[z] * tmp[A->col[z]] * tmp[A->rowidx[z]];
            }
        }
        else if (scaling == Magma_minproduct_UNITDIAG ) {
            // scale to unit diagonal
            CHECK( magma_minproduct_zmalloc_cpu( &tmp, A->num_rows ));
            for( magma_minproduct_int_t z=0; z<A->num_rows; z++ ) {
                magma_minproductDoubleComplex s = MAGMA_minproduct_Z_MAKE( 0.0, 0.0 );
                for( magma_minproduct_int_t f=A->row[z]; f<A->row[z+1]; f++ ) {
                    if ( A->col[f]== z ) {
                        // add some identity matrix
                        //A->val[f] = A->val[f] +  MAGMA_minproduct_Z_MAKE( 100000.0, 0.0 );
                        s = A->val[f];
                    }
                }
                if ( s == MAGMA_minproduct_Z_MAKE( 0.0, 0.0 ) ){
                    printf("error: zero diagonal element.\n");
                    info = MAGMA_minproduct_ERR;
                }
                tmp[z] = MAGMA_minproduct_Z_MAKE( 1.0/sqrt(  MAGMA_minproduct_Z_REAL( s )  ), 0.0 );
                   
            }
            for( magma_minproduct_int_t z=0; z<A->nnz; z++ ) {
                A->val[z] = A->val[z] * tmp[A->col[z]] * tmp[A->rowidx[z]];
            }
        }
        else{
            printf( "error: scaling not supported.\n" );
            info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        }
    }
    else {
        magma_minproduct_storage_t A_storage = A->storage_type;
        magma_minproduct_location_t A_location = A->memory_location;
        CHECK( magma_minproduct_zmtransfer( *A, &hA, A->memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_zmconvert( hA, &CSRA, hA.storage_type, Magma_minproduct_CSRCOO, queue ));

        CHECK( magma_minproduct_zmscale( &CSRA, scaling, queue ));

        magma_minproduct_zmfree( &hA, queue );
        magma_minproduct_zmfree( A, queue );
        CHECK( magma_minproduct_zmconvert( CSRA, &hA, Magma_minproduct_CSRCOO, A_storage, queue ));
        CHECK( magma_minproduct_zmtransfer( hA, A, Magma_minproduct_CPU, A_location, queue ));
    }
    
cleanup:
    magma_minproduct_free_cpu( tmp );
    magma_minproduct_zmfree( &hA, queue );
    magma_minproduct_zmfree( &CSRA, queue );
    return info;
}


/**
    Purpose
    -------

    Adds a multiple of the Identity matrix to a matrix: A = A+add * I

    Arguments
    ---------

    @param[in,out]
    A           magma_minproduct_z_matrix*
                input/output matrix

    @param[in]
    add         magma_minproductDoubleComplex
                scaling for the identity matrix
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zmdiagadd(
    magma_minproduct_z_matrix *A,
    magma_minproductDoubleComplex add,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_z_matrix hA={Magma_minproduct_CSR}, CSRA={Magma_minproduct_CSR};
    
    if ( A->memory_location == Magma_minproduct_CPU && A->storage_type == Magma_minproduct_CSRCOO ) {
        for( magma_minproduct_int_t z=0; z<A->nnz; z++ ) {
            if ( A->col[z]== A->rowidx[z] ) {
                // add some identity matrix
                A->val[z] = A->val[z] +  add;
            }
        }
    }
    else {
        magma_minproduct_storage_t A_storage = A->storage_type;
        magma_minproduct_location_t A_location = A->memory_location;
        CHECK( magma_minproduct_zmtransfer( *A, &hA, A->memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_zmconvert( hA, &CSRA, hA.storage_type, Magma_minproduct_CSRCOO, queue ));

        CHECK( magma_minproduct_zmdiagadd( &CSRA, add, queue ));

        magma_minproduct_zmfree( &hA, queue );
        magma_minproduct_zmfree( A, queue );
        CHECK( magma_minproduct_zmconvert( CSRA, &hA, Magma_minproduct_CSRCOO, A_storage, queue ));
        CHECK( magma_minproduct_zmtransfer( hA, A, Magma_minproduct_CPU, A_location, queue ));
    }
    
cleanup:
    magma_minproduct_zmfree( &hA, queue );
    magma_minproduct_zmfree( &CSRA, queue );
    return info;
}



