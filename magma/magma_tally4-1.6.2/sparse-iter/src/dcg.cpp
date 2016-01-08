/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zcg.cpp normal z -> d, Sun May  3 11:22:59 2015
*/

#include "common_magma_tally4sparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a real symmetric N-by-N positive definite matrix A.
    This is a GPU implementation of the Conjugate Gradient method.

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

    @ingroup magma_tally4sparse_dposv
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dcg(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, magma_tally4_d_matrix *x,
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_tally4_CG;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_tally4_SUCCESS;

    // local variables
    double c_zero = MAGMA_tally4_D_ZERO, c_one = MAGMA_tally4_D_ONE;
    
    magma_tally4_int_t dofs = A.num_rows * b.num_cols;

    // GPU workspace
    magma_tally4_d_matrix r={Magma_tally4_CSR}, p={Magma_tally4_CSR}, q={Magma_tally4_CSR};
    CHECK( magma_tally4_dvinit( &r, Magma_tally4_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally4_dvinit( &p, Magma_tally4_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally4_dvinit( &q, Magma_tally4_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    
    // solver variables
    double alpha, beta;
    double nom, nom0, r0, betanom, betanomsq, den;

    // solver setup
    CHECK(  magma_tally4_dresidualvec( A, b, *x, &r, &nom0, queue));
    magma_tally4_dcopy( dofs, r.dval, 1, p.dval, 1 );                    // p = r
    betanom = nom0;
    nom  = nom0 * nom0;                                // nom = r' * r
    CHECK( magma_tally4_d_spmv( c_one, A, p, c_zero, q, queue ));             // q = A p
    den = MAGMA_tally4_D_REAL( magma_tally4_ddot(dofs, p.dval, 1, q.dval, 1) );// den = p dot q
    solver_par->init_res = nom0;
    
    if ( (r0 = nom * solver_par->epsilon) < ATOLERANCE )
        r0 = ATOLERANCE;
    if ( nom < r0 ) {
        solver_par->final_res = solver_par->init_res;
        solver_par->iter_res = solver_par->init_res;
        goto cleanup;
    }
    // check positive definite
    if (den <= 0.0) {
        printf("Operator A is not postive definite. (Ar,r) = %f\n", den);
        magma_tally4blasSetKernelStream( orig_queue );
        info = MAGMA_tally4_NONSPD; 
        goto cleanup;
    }

    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_tally4_sync_wtime( queue );
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = (real_Double_t)nom0;
        solver_par->timing[0] = 0.0;
    }
    
    solver_par->numiter = 0;
    // start iteration
    do
    {
        solver_par->numiter++;
        alpha = MAGMA_tally4_D_MAKE(nom/den, 0.);
        magma_tally4_daxpy(dofs,  alpha, p.dval, 1, x->dval, 1);     // x = x + alpha p
        magma_tally4_daxpy(dofs, -alpha, q.dval, 1, r.dval, 1);      // r = r - alpha q
        betanom = magma_tally4_dnrm2(dofs, r.dval, 1);             // betanom = || r ||
        betanomsq = betanom * betanom;                      // betanoms = r' * r

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_tally4_sync_wtime( queue );
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }

        if (  betanom  < r0 ) {
            break;
        }

        beta = MAGMA_tally4_D_MAKE(betanomsq/nom, 0.);           // beta = betanoms/nom
        magma_tally4_dscal(dofs, beta, p.dval, 1);                // p = beta*p
        magma_tally4_daxpy(dofs, c_one, r.dval, 1, p.dval, 1);     // p = p + r
        CHECK( magma_tally4_d_spmv( c_one, A, p, c_zero, q, queue ));   // q = A p
        den = MAGMA_tally4_D_REAL(magma_tally4_ddot(dofs, p.dval, 1, q.dval, 1));
                // den = p dot q
        nom = betanomsq;
    }
    while ( solver_par->numiter+1 <= solver_par->maxiter );
    
    tempo2 = magma_tally4_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    double residual;
    CHECK(  magma_tally4_dresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter ) {
        solver_par->info = MAGMA_tally4_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_tally4_SLOW_CONVERGENCE;
        if( solver_par->iter_res < solver_par->epsilon*solver_par->init_res ){
            info = MAGMA_tally4_SUCCESS;
        }
    }
    else {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_tally4_DIVERGENCE;
    }
    
cleanup:
    magma_tally4_dmfree(&r, queue );
    magma_tally4_dmfree(&p, queue );
    magma_tally4_dmfree(&q, queue );

    magma_tally4blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally4_dcg */


