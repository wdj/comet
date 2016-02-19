/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zbicgstab.cpp normal z -> s, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a general matrix.
    This is a GPU implementation of the Biconjugate Gradient Stabilized method.

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
magma_tally3_sbicgstab(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, magma_tally3_s_matrix *x,
    magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_tally3_BICGSTAB;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_tally3_SUCCESS;

    // some useful variables
    float c_zero = MAGMA_tally3_S_ZERO, c_one = MAGMA_tally3_S_ONE,
                                            c_mone = MAGMA_tally3_S_NEG_ONE;
    
    magma_tally3_int_t dofs = A.num_rows * b.num_cols;

    // workspace
    magma_tally3_s_matrix r={Magma_tally3_CSR}, rr={Magma_tally3_CSR}, p={Magma_tally3_CSR}, v={Magma_tally3_CSR}, s={Magma_tally3_CSR}, t={Magma_tally3_CSR};
    CHECK( magma_tally3_svinit( &r, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_svinit( &rr,Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_svinit( &p, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_svinit( &v, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_svinit( &s, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_svinit( &t, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));

    
    // solver variables
    float alpha, beta, omega, rho_old, rho_new;
    float nom, betanom, nom0, r0, den, res;

    // solver setup
    CHECK(  magma_tally3_sresidualvec( A, b, *x, &r, &nom0, queue));
    magma_tally3_scopy( dofs, r.dval, 1, rr.dval, 1 );                  // rr = r
    betanom = nom0;
    nom = nom0*nom0;
    rho_new = magma_tally3_sdot( dofs, r.dval, 1, r.dval, 1 );             // rho=<rr,r>
    rho_old = omega = alpha = MAGMA_tally3_S_MAKE( 1.0, 0. );
    solver_par->init_res = nom0;

    CHECK( magma_tally3_s_spmv( c_one, A, r, c_zero, v, queue ));              // z = A r
    den = MAGMA_tally3_S_REAL( magma_tally3_sdot(dofs, v.dval, 1, r.dval, 1) ); // den = z' * r

    if ( (r0 = nom * solver_par->epsilon) < ATOLERANCE )
        r0 = ATOLERANCE;
    if ( nom < r0 ) {
        solver_par->final_res = solver_par->init_res;
        solver_par->iter_res = solver_par->init_res;
        goto cleanup;
    }

    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_tally3_sync_wtime( queue );
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = nom0;
        solver_par->timing[0] = 0.0;
    }

    solver_par->numiter = 0;
    // start iteration
    do
    {
        solver_par->numiter++;

        rho_new = magma_tally3_sdot( dofs, rr.dval, 1, r.dval, 1 );  // rho=<rr,r>
        beta = rho_new/rho_old * alpha/omega;   // beta=rho/rho_old *alpha/omega
        magma_tally3_sscal( dofs, beta, p.dval, 1 );                 // p = beta*p
        magma_tally3_saxpy( dofs, c_mone * omega * beta, v.dval, 1 , p.dval, 1 );
                                                        // p = p-omega*beta*v
        magma_tally3_saxpy( dofs, c_one, r.dval, 1, p.dval, 1 );      // p = p+r
        CHECK( magma_tally3_s_spmv( c_one, A, p, c_zero, v, queue ));      // v = Ap

        alpha = rho_new / magma_tally3_sdot( dofs, rr.dval, 1, v.dval, 1 );
        magma_tally3_scopy( dofs, r.dval, 1 , s.dval, 1 );            // s=r
        magma_tally3_saxpy( dofs, c_mone * alpha, v.dval, 1 , s.dval, 1 ); // s=s-alpha*v

        CHECK( magma_tally3_s_spmv( c_one, A, s, c_zero, t, queue ));       // t=As
        omega = magma_tally3_sdot( dofs, t.dval, 1, s.dval, 1 )   // omega = <s,t>/<t,t>
                   / magma_tally3_sdot( dofs, t.dval, 1, t.dval, 1 );

        magma_tally3_saxpy( dofs, alpha, p.dval, 1 , x->dval, 1 );     // x=x+alpha*p
        magma_tally3_saxpy( dofs, omega, s.dval, 1 , x->dval, 1 );     // x=x+omega*s

        magma_tally3_scopy( dofs, s.dval, 1 , r.dval, 1 );             // r=s
        magma_tally3_saxpy( dofs, c_mone * omega, t.dval, 1 , r.dval, 1 ); // r=r-omega*t
        res = betanom = magma_tally3_snrm2( dofs, r.dval, 1 );

        nom = betanom*betanom;
        rho_old = rho_new;                                    // rho_old=rho

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_tally3_sync_wtime( queue );
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }

        if ( res/nom0  < solver_par->epsilon ) {
            break;
        }
    }
    while ( solver_par->numiter+1 <= solver_par->maxiter );
    
    tempo2 = magma_tally3_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    float residual;
    CHECK(  magma_tally3_sresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->iter_res = res;
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter ) {
        solver_par->info = MAGMA_tally3_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_tally3_SLOW_CONVERGENCE;
        if( solver_par->iter_res < solver_par->epsilon*solver_par->init_res ){
            info = MAGMA_tally3_SUCCESS;
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
        info = MAGMA_tally3_DIVERGENCE;
    }
    
cleanup:
    magma_tally3_smfree(&r, queue );
    magma_tally3_smfree(&rr, queue );
    magma_tally3_smfree(&p, queue );
    magma_tally3_smfree(&v, queue );
    magma_tally3_smfree(&s, queue );
    magma_tally3_smfree(&t, queue );

    magma_tally3blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally3_sbicgstab */


