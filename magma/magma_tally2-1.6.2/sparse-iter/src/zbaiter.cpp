/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @precisions normal z -> s d c
*/

#include "common_magma_tally2sparse.h"

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
    A           magma_tally2_z_matrix
                input matrix A

    @param[in]
    b           magma_tally2_z_matrix
                RHS b

    @param[in,out]
    x           magma_tally2_z_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally2_z_solver_par*
                solver parameters

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zgesv
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zbaiter(
    magma_tally2_z_matrix A,
    magma_tally2_z_matrix b,
    magma_tally2_z_matrix *x,
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
        
    // prepare solver feedback
    solver_par->solver = Magma_tally2_BAITER;
    solver_par->info = MAGMA_tally2_SUCCESS;

    // initial residual
    real_Double_t tempo1, tempo2;
    double residual;
    magma_tally2_int_t localiter = 1;
    
    magma_tally2_z_matrix Ah={Magma_tally2_CSR}, ACSR={Magma_tally2_CSR}, A_d={Magma_tally2_CSR}, D={Magma_tally2_CSR}, 
                    R={Magma_tally2_CSR}, D_d={Magma_tally2_CSR}, R_d={Magma_tally2_CSR};
    
    CHECK( magma_tally2_zresidual( A, b, *x, &residual, queue ));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;



    CHECK( magma_tally2_zmtransfer( A, &Ah, A.memory_location, Magma_tally2_CPU, queue ));
    CHECK( magma_tally2_zmconvert( Ah, &ACSR, Ah.storage_type, Magma_tally2_CSR, queue ));

    CHECK( magma_tally2_zmtransfer( ACSR, &A_d, Magma_tally2_CPU, Magma_tally2_DEV, queue ));

    // setup
    CHECK( magma_tally2_zcsrsplit( 256, ACSR, &D, &R, queue ));
    CHECK( magma_tally2_zmtransfer( D, &D_d, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    CHECK( magma_tally2_zmtransfer( R, &R_d, Magma_tally2_CPU, Magma_tally2_DEV, queue ));


    tempo1 = magma_tally2_sync_wtime( queue );

    // block-asynchronous iteration iterator
    for( int iter=0; iter<solver_par->maxiter; iter++)
        CHECK( magma_tally2_zbajac_csr( localiter, D_d, R_d, b, x, queue ));

    tempo2 = magma_tally2_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_tally2_zresidual( A, b, *x, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res ){
        info = MAGMA_tally2_SUCCESS;
    }
    else{
        info = MAGMA_tally2_DIVERGENCE;
    }
    
cleanup:
    magma_tally2_zmfree(&D, queue );
    magma_tally2_zmfree(&R, queue );
    magma_tally2_zmfree(&D_d, queue );
    magma_tally2_zmfree(&R_d, queue );
    magma_tally2_zmfree(&A_d, queue );
    magma_tally2_zmfree(&ACSR, queue );
    magma_tally2_zmfree(&Ah, queue );

    solver_par->info = info;
    return info;
}   /* magma_tally2_zbaiter */

