/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zbaiter.cpp normal z -> s, Sun May  3 11:22:59 2015
*/

#include "common_magma_minproductsparse.h"

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
    A           magma_minproduct_s_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_s_matrix
                RHS b

    @param[in,out]
    x           magma_minproduct_s_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_minproduct_s_solver_par*
                solver parameters

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_sgesv
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_sbaiter(
    magma_minproduct_s_matrix A,
    magma_minproduct_s_matrix b,
    magma_minproduct_s_matrix *x,
    magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
        
    // prepare solver feedback
    solver_par->solver = Magma_minproduct_BAITER;
    solver_par->info = MAGMA_minproduct_SUCCESS;

    // initial residual
    real_Double_t tempo1, tempo2;
    float residual;
    magma_minproduct_int_t localiter = 1;
    
    magma_minproduct_s_matrix Ah={Magma_minproduct_CSR}, ACSR={Magma_minproduct_CSR}, A_d={Magma_minproduct_CSR}, D={Magma_minproduct_CSR}, 
                    R={Magma_minproduct_CSR}, D_d={Magma_minproduct_CSR}, R_d={Magma_minproduct_CSR};
    
    CHECK( magma_minproduct_sresidual( A, b, *x, &residual, queue ));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;



    CHECK( magma_minproduct_smtransfer( A, &Ah, A.memory_location, Magma_minproduct_CPU, queue ));
    CHECK( magma_minproduct_smconvert( Ah, &ACSR, Ah.storage_type, Magma_minproduct_CSR, queue ));

    CHECK( magma_minproduct_smtransfer( ACSR, &A_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

    // setup
    CHECK( magma_minproduct_scsrsplit( 256, ACSR, &D, &R, queue ));
    CHECK( magma_minproduct_smtransfer( D, &D_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_smtransfer( R, &R_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));


    tempo1 = magma_minproduct_sync_wtime( queue );

    // block-asynchronous iteration iterator
    for( int iter=0; iter<solver_par->maxiter; iter++)
        CHECK( magma_minproduct_sbajac_csr( localiter, D_d, R_d, b, x, queue ));

    tempo2 = magma_minproduct_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_minproduct_sresidual( A, b, *x, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res ){
        info = MAGMA_minproduct_SUCCESS;
    }
    else{
        info = MAGMA_minproduct_DIVERGENCE;
    }
    
cleanup:
    magma_minproduct_smfree(&D, queue );
    magma_minproduct_smfree(&R, queue );
    magma_minproduct_smfree(&D_d, queue );
    magma_minproduct_smfree(&R_d, queue );
    magma_minproduct_smfree(&A_d, queue );
    magma_minproduct_smfree(&ACSR, queue );
    magma_minproduct_smfree(&Ah, queue );

    solver_par->info = info;
    return info;
}   /* magma_minproduct_sbaiter */

