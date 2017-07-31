/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zcg.cpp normal z -> s, Sun May  3 11:22:59 2015
*/

#include "common_magma_tally2sparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


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
    A           magma_tally2_s_matrix
                input matrix A

    @param[in]
    b           magma_tally2_s_matrix
                RHS b

    @param[in,out]
    x           magma_tally2_s_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally2_s_solver_par*
                solver parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_sposv
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_scg(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, magma_tally2_s_matrix *x,
    magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue=NULL;
    magma_tally2blasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_tally2_CG;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_tally2_SUCCESS;

    // local variables
    float c_zero = MAGMA_tally2_S_ZERO, c_one = MAGMA_tally2_S_ONE;
    
    magma_tally2_int_t dofs = A.num_rows * b.num_cols;

    // GPU workspace
    magma_tally2_s_matrix r={Magma_tally2_CSR}, p={Magma_tally2_CSR}, q={Magma_tally2_CSR};
    CHECK( magma_tally2_svinit( &r, Magma_tally2_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally2_svinit( &p, Magma_tally2_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally2_svinit( &q, Magma_tally2_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    
    // solver variables
    float alpha, beta;
    float nom, nom0, r0, betanom, betanomsq, den;

    // solver setup
    CHECK(  magma_tally2_sresidualvec( A, b, *x, &r, &nom0, queue));
    magma_tally2_scopy( dofs, r.dval, 1, p.dval, 1 );                    // p = r
    betanom = nom0;
    nom  = nom0 * nom0;                                // nom = r' * r
    CHECK( magma_tally2_s_spmv( c_one, A, p, c_zero, q, queue ));             // q = A p
    den = MAGMA_tally2_S_REAL( magma_tally2_sdot(dofs, p.dval, 1, q.dval, 1) );// den = p dot q
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
        magma_tally2blasSetKernelStream( orig_queue );
        info = MAGMA_tally2_NONSPD; 
        goto cleanup;
    }

    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_tally2_sync_wtime( queue );
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = (real_Double_t)nom0;
        solver_par->timing[0] = 0.0;
    }
    
    solver_par->numiter = 0;
    // start iteration
    do
    {
        solver_par->numiter++;
        alpha = MAGMA_tally2_S_MAKE(nom/den, 0.);
        magma_tally2_saxpy(dofs,  alpha, p.dval, 1, x->dval, 1);     // x = x + alpha p
        magma_tally2_saxpy(dofs, -alpha, q.dval, 1, r.dval, 1);      // r = r - alpha q
        betanom = magma_tally2_snrm2(dofs, r.dval, 1);             // betanom = || r ||
        betanomsq = betanom * betanom;                      // betanoms = r' * r

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_tally2_sync_wtime( queue );
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

        beta = MAGMA_tally2_S_MAKE(betanomsq/nom, 0.);           // beta = betanoms/nom
        magma_tally2_sscal(dofs, beta, p.dval, 1);                // p = beta*p
        magma_tally2_saxpy(dofs, c_one, r.dval, 1, p.dval, 1);     // p = p + r
        CHECK( magma_tally2_s_spmv( c_one, A, p, c_zero, q, queue ));   // q = A p
        den = MAGMA_tally2_S_REAL(magma_tally2_sdot(dofs, p.dval, 1, q.dval, 1));
                // den = p dot q
        nom = betanomsq;
    }
    while ( solver_par->numiter+1 <= solver_par->maxiter );
    
    tempo2 = magma_tally2_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    float residual;
    CHECK(  magma_tally2_sresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter ) {
        solver_par->info = MAGMA_tally2_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_tally2_SLOW_CONVERGENCE;
        if( solver_par->iter_res < solver_par->epsilon*solver_par->init_res ){
            info = MAGMA_tally2_SUCCESS;
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
        info = MAGMA_tally2_DIVERGENCE;
    }
    
cleanup:
    magma_tally2_smfree(&r, queue );
    magma_tally2_smfree(&p, queue );
    magma_tally2_smfree(&q, queue );

    magma_tally2blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally2_scg */


