/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zbaiter.cpp normal z -> d, Sun May  3 11:22:59 2015
*/

#include "common_magma_tally4sparse.h"

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
    A           magma_tally4_d_matrix
                input matrix A

    @param[in]
    b           magma_tally4_d_matrix
                RHS b

    @param[in,out]
    x           magma_tally4_d_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally4_d_solver_par*
                solver parameters

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dgesv
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dbaiter(
    magma_tally4_d_matrix A,
    magma_tally4_d_matrix b,
    magma_tally4_d_matrix *x,
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
        
    // prepare solver feedback
    solver_par->solver = Magma_tally4_BAITER;
    solver_par->info = MAGMA_tally4_SUCCESS;

    // initial residual
    real_Double_t tempo1, tempo2;
    double residual;
    magma_tally4_int_t localiter = 1;
    
    magma_tally4_d_matrix Ah={Magma_tally4_CSR}, ACSR={Magma_tally4_CSR}, A_d={Magma_tally4_CSR}, D={Magma_tally4_CSR}, 
                    R={Magma_tally4_CSR}, D_d={Magma_tally4_CSR}, R_d={Magma_tally4_CSR};
    
    CHECK( magma_tally4_dresidual( A, b, *x, &residual, queue ));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;



    CHECK( magma_tally4_dmtransfer( A, &Ah, A.memory_location, Magma_tally4_CPU, queue ));
    CHECK( magma_tally4_dmconvert( Ah, &ACSR, Ah.storage_type, Magma_tally4_CSR, queue ));

    CHECK( magma_tally4_dmtransfer( ACSR, &A_d, Magma_tally4_CPU, Magma_tally4_DEV, queue ));

    // setup
    CHECK( magma_tally4_dcsrsplit( 256, ACSR, &D, &R, queue ));
    CHECK( magma_tally4_dmtransfer( D, &D_d, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
    CHECK( magma_tally4_dmtransfer( R, &R_d, Magma_tally4_CPU, Magma_tally4_DEV, queue ));


    tempo1 = magma_tally4_sync_wtime( queue );

    // block-asynchronous iteration iterator
    for( int iter=0; iter<solver_par->maxiter; iter++)
        CHECK( magma_tally4_dbajac_csr( localiter, D_d, R_d, b, x, queue ));

    tempo2 = magma_tally4_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_tally4_dresidual( A, b, *x, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res ){
        info = MAGMA_tally4_SUCCESS;
    }
    else{
        info = MAGMA_tally4_DIVERGENCE;
    }
    
cleanup:
    magma_tally4_dmfree(&D, queue );
    magma_tally4_dmfree(&R, queue );
    magma_tally4_dmfree(&D_d, queue );
    magma_tally4_dmfree(&R_d, queue );
    magma_tally4_dmfree(&A_d, queue );
    magma_tally4_dmfree(&ACSR, queue );
    magma_tally4_dmfree(&Ah, queue );

    solver_par->info = info;
    return info;
}   /* magma_tally4_dbaiter */

