/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @precisions normal z -> s d c
*/

#include "common_magma_minproductsparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * x = b
    via the block asynchronous iteration method on GPU.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_z_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_z_matrix
                RHS b

    @param[in,out]
    x           magma_minproduct_z_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_minproduct_z_solver_par*
                solver parameters

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgesv
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zbaiter(
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix b,
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
        
    // prepare solver feedback
    solver_par->solver = Magma_minproduct_BAITER;
    solver_par->info = MAGMA_minproduct_SUCCESS;

    // initial residual
    real_Double_t tempo1, tempo2;
    double residual;
    magma_minproduct_int_t localiter = 1;
    
    magma_minproduct_z_matrix Ah={Magma_minproduct_CSR}, ACSR={Magma_minproduct_CSR}, A_d={Magma_minproduct_CSR}, D={Magma_minproduct_CSR}, 
                    R={Magma_minproduct_CSR}, D_d={Magma_minproduct_CSR}, R_d={Magma_minproduct_CSR};
    
    CHECK( magma_minproduct_zresidual( A, b, *x, &residual, queue ));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;



    CHECK( magma_minproduct_zmtransfer( A, &Ah, A.memory_location, Magma_minproduct_CPU, queue ));
    CHECK( magma_minproduct_zmconvert( Ah, &ACSR, Ah.storage_type, Magma_minproduct_CSR, queue ));

    CHECK( magma_minproduct_zmtransfer( ACSR, &A_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

    // setup
    CHECK( magma_minproduct_zcsrsplit( 256, ACSR, &D, &R, queue ));
    CHECK( magma_minproduct_zmtransfer( D, &D_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_zmtransfer( R, &R_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));


    tempo1 = magma_minproduct_sync_wtime( queue );

    // block-asynchronous iteration iterator
    for( int iter=0; iter<solver_par->maxiter; iter++)
        CHECK( magma_minproduct_zbajac_csr( localiter, D_d, R_d, b, x, queue ));

    tempo2 = magma_minproduct_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_minproduct_zresidual( A, b, *x, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res ){
        info = MAGMA_minproduct_SUCCESS;
    }
    else{
        info = MAGMA_minproduct_DIVERGENCE;
    }
    
cleanup:
    magma_minproduct_zmfree(&D, queue );
    magma_minproduct_zmfree(&R, queue );
    magma_minproduct_zmfree(&D_d, queue );
    magma_minproduct_zmfree(&R_d, queue );
    magma_minproduct_zmfree(&A_d, queue );
    magma_minproduct_zmfree(&ACSR, queue );
    magma_minproduct_zmfree(&Ah, queue );

    solver_par->info = info;
    return info;
}   /* magma_minproduct_zbaiter */

