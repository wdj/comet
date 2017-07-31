/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zbaiter.cpp normal z -> c, Sun May  3 11:22:59 2015
*/

#include "common_magma_tally2sparse.h"

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
    A           magma_tally2_c_matrix
                input matrix A

    @param[in]
    b           magma_tally2_c_matrix
                RHS b

    @param[in,out]
    x           magma_tally2_c_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally2_c_solver_par*
                solver parameters

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_cgesv
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_cbaiter(
    magma_tally2_c_matrix A,
    magma_tally2_c_matrix b,
    magma_tally2_c_matrix *x,
    magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
        
    // prepare solver feedback
    solver_par->solver = Magma_tally2_BAITER;
    solver_par->info = MAGMA_tally2_SUCCESS;

    // initial residual
    real_Double_t tempo1, tempo2;
    float residual;
    magma_tally2_int_t localiter = 1;
    
    magma_tally2_c_matrix Ah={Magma_tally2_CSR}, ACSR={Magma_tally2_CSR}, A_d={Magma_tally2_CSR}, D={Magma_tally2_CSR}, 
                    R={Magma_tally2_CSR}, D_d={Magma_tally2_CSR}, R_d={Magma_tally2_CSR};
    
    CHECK( magma_tally2_cresidual( A, b, *x, &residual, queue ));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;



    CHECK( magma_tally2_cmtransfer( A, &Ah, A.memory_location, Magma_tally2_CPU, queue ));
    CHECK( magma_tally2_cmconvert( Ah, &ACSR, Ah.storage_type, Magma_tally2_CSR, queue ));

    CHECK( magma_tally2_cmtransfer( ACSR, &A_d, Magma_tally2_CPU, Magma_tally2_DEV, queue ));

    // setup
    CHECK( magma_tally2_ccsrsplit( 256, ACSR, &D, &R, queue ));
    CHECK( magma_tally2_cmtransfer( D, &D_d, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    CHECK( magma_tally2_cmtransfer( R, &R_d, Magma_tally2_CPU, Magma_tally2_DEV, queue ));


    tempo1 = magma_tally2_sync_wtime( queue );

    // block-asynchronous iteration iterator
    for( int iter=0; iter<solver_par->maxiter; iter++)
        CHECK( magma_tally2_cbajac_csr( localiter, D_d, R_d, b, x, queue ));

    tempo2 = magma_tally2_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_tally2_cresidual( A, b, *x, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res ){
        info = MAGMA_tally2_SUCCESS;
    }
    else{
        info = MAGMA_tally2_DIVERGENCE;
    }
    
cleanup:
    magma_tally2_cmfree(&D, queue );
    magma_tally2_cmfree(&R, queue );
    magma_tally2_cmfree(&D_d, queue );
    magma_tally2_cmfree(&R_d, queue );
    magma_tally2_cmfree(&A_d, queue );
    magma_tally2_cmfree(&ACSR, queue );
    magma_tally2_cmfree(&Ah, queue );

    solver_par->info = info;
    return info;
}   /* magma_tally2_cbaiter */

