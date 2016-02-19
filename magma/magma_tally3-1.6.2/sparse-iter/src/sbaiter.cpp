/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zbaiter.cpp normal z -> s, Sun May  3 11:22:59 2015
*/

#include "common_magma_tally3sparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * x = b
    via the block asynchronous iteration method on GPU.

    Arguments
    ---------

    @param[in]
    A           magma_tally3_s_matrix
                input matrix A

    @param[in]
    b           magma_tally3_s_matrix
                RHS b

    @param[in,out]
    x           magma_tally3_s_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally3_s_solver_par*
                solver parameters

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_sgesv
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_sbaiter(
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix b,
    magma_tally3_s_matrix *x,
    magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
        
    // prepare solver feedback
    solver_par->solver = Magma_tally3_BAITER;
    solver_par->info = MAGMA_tally3_SUCCESS;

    // initial residual
    real_Double_t tempo1, tempo2;
    float residual;
    magma_tally3_int_t localiter = 1;
    
    magma_tally3_s_matrix Ah={Magma_tally3_CSR}, ACSR={Magma_tally3_CSR}, A_d={Magma_tally3_CSR}, D={Magma_tally3_CSR}, 
                    R={Magma_tally3_CSR}, D_d={Magma_tally3_CSR}, R_d={Magma_tally3_CSR};
    
    CHECK( magma_tally3_sresidual( A, b, *x, &residual, queue ));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;



    CHECK( magma_tally3_smtransfer( A, &Ah, A.memory_location, Magma_tally3_CPU, queue ));
    CHECK( magma_tally3_smconvert( Ah, &ACSR, Ah.storage_type, Magma_tally3_CSR, queue ));

    CHECK( magma_tally3_smtransfer( ACSR, &A_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));

    // setup
    CHECK( magma_tally3_scsrsplit( 256, ACSR, &D, &R, queue ));
    CHECK( magma_tally3_smtransfer( D, &D_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
    CHECK( magma_tally3_smtransfer( R, &R_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));


    tempo1 = magma_tally3_sync_wtime( queue );

    // block-asynchronous iteration iterator
    for( int iter=0; iter<solver_par->maxiter; iter++)
        CHECK( magma_tally3_sbajac_csr( localiter, D_d, R_d, b, x, queue ));

    tempo2 = magma_tally3_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_tally3_sresidual( A, b, *x, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res ){
        info = MAGMA_tally3_SUCCESS;
    }
    else{
        info = MAGMA_tally3_DIVERGENCE;
    }
    
cleanup:
    magma_tally3_smfree(&D, queue );
    magma_tally3_smfree(&R, queue );
    magma_tally3_smfree(&D_d, queue );
    magma_tally3_smfree(&R_d, queue );
    magma_tally3_smfree(&A_d, queue );
    magma_tally3_smfree(&ACSR, queue );
    magma_tally3_smfree(&Ah, queue );

    solver_par->info = info;
    return info;
}   /* magma_tally3_sbaiter */

