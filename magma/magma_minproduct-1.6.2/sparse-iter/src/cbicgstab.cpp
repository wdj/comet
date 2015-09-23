/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zbicgstab.cpp normal z -> c, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"

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
    A           magma_minproduct_c_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_c_matrix
                RHS b

    @param[in,out]
    x           magma_minproduct_c_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_minproduct_c_solver_par*
                solver parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_cgesv
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cbicgstab(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, magma_minproduct_c_matrix *x,
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_minproduct_BICGSTAB;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_minproduct_SUCCESS;

    // some useful variables
    magma_minproductFloatComplex c_zero = MAGMA_minproduct_C_ZERO, c_one = MAGMA_minproduct_C_ONE,
                                            c_mone = MAGMA_minproduct_C_NEG_ONE;
    
    magma_minproduct_int_t dofs = A.num_rows * b.num_cols;

    // workspace
    magma_minproduct_c_matrix r={Magma_minproduct_CSR}, rr={Magma_minproduct_CSR}, p={Magma_minproduct_CSR}, v={Magma_minproduct_CSR}, s={Magma_minproduct_CSR}, t={Magma_minproduct_CSR};
    CHECK( magma_minproduct_cvinit( &r, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_cvinit( &rr,Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_cvinit( &p, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_cvinit( &v, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_cvinit( &s, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_cvinit( &t, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));

    
    // solver variables
    magma_minproductFloatComplex alpha, beta, omega, rho_old, rho_new;
    float nom, betanom, nom0, r0, den, res;

    // solver setup
    CHECK(  magma_minproduct_cresidualvec( A, b, *x, &r, &nom0, queue));
    magma_minproduct_ccopy( dofs, r.dval, 1, rr.dval, 1 );                  // rr = r
    betanom = nom0;
    nom = nom0*nom0;
    rho_new = magma_minproduct_cdotc( dofs, r.dval, 1, r.dval, 1 );             // rho=<rr,r>
    rho_old = omega = alpha = MAGMA_minproduct_C_MAKE( 1.0, 0. );
    solver_par->init_res = nom0;

    CHECK( magma_minproduct_c_spmv( c_one, A, r, c_zero, v, queue ));              // z = A r
    den = MAGMA_minproduct_C_REAL( magma_minproduct_cdotc(dofs, v.dval, 1, r.dval, 1) ); // den = z' * r

    if ( (r0 = nom * solver_par->epsilon) < ATOLERANCE )
        r0 = ATOLERANCE;
    if ( nom < r0 ) {
        solver_par->final_res = solver_par->init_res;
        solver_par->iter_res = solver_par->init_res;
        goto cleanup;
    }

    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_minproduct_sync_wtime( queue );
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = nom0;
        solver_par->timing[0] = 0.0;
    }

    solver_par->numiter = 0;
    // start iteration
    do
    {
        solver_par->numiter++;

        rho_new = magma_minproduct_cdotc( dofs, rr.dval, 1, r.dval, 1 );  // rho=<rr,r>
        beta = rho_new/rho_old * alpha/omega;   // beta=rho/rho_old *alpha/omega
        magma_minproduct_cscal( dofs, beta, p.dval, 1 );                 // p = beta*p
        magma_minproduct_caxpy( dofs, c_mone * omega * beta, v.dval, 1 , p.dval, 1 );
                                                        // p = p-omega*beta*v
        magma_minproduct_caxpy( dofs, c_one, r.dval, 1, p.dval, 1 );      // p = p+r
        CHECK( magma_minproduct_c_spmv( c_one, A, p, c_zero, v, queue ));      // v = Ap

        alpha = rho_new / magma_minproduct_cdotc( dofs, rr.dval, 1, v.dval, 1 );
        magma_minproduct_ccopy( dofs, r.dval, 1 , s.dval, 1 );            // s=r
        magma_minproduct_caxpy( dofs, c_mone * alpha, v.dval, 1 , s.dval, 1 ); // s=s-alpha*v

        CHECK( magma_minproduct_c_spmv( c_one, A, s, c_zero, t, queue ));       // t=As
        omega = magma_minproduct_cdotc( dofs, t.dval, 1, s.dval, 1 )   // omega = <s,t>/<t,t>
                   / magma_minproduct_cdotc( dofs, t.dval, 1, t.dval, 1 );

        magma_minproduct_caxpy( dofs, alpha, p.dval, 1 , x->dval, 1 );     // x=x+alpha*p
        magma_minproduct_caxpy( dofs, omega, s.dval, 1 , x->dval, 1 );     // x=x+omega*s

        magma_minproduct_ccopy( dofs, s.dval, 1 , r.dval, 1 );             // r=s
        magma_minproduct_caxpy( dofs, c_mone * omega, t.dval, 1 , r.dval, 1 ); // r=r-omega*t
        res = betanom = magma_minproduct_scnrm2( dofs, r.dval, 1 );

        nom = betanom*betanom;
        rho_old = rho_new;                                    // rho_old=rho

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_minproduct_sync_wtime( queue );
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
    
    tempo2 = magma_minproduct_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    float residual;
    CHECK(  magma_minproduct_cresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->iter_res = res;
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter ) {
        solver_par->info = MAGMA_minproduct_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_minproduct_SLOW_CONVERGENCE;
        if( solver_par->iter_res < solver_par->epsilon*solver_par->init_res ){
            info = MAGMA_minproduct_SUCCESS;
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
        info = MAGMA_minproduct_DIVERGENCE;
    }
    
cleanup:
    magma_minproduct_cmfree(&r, queue );
    magma_minproduct_cmfree(&rr, queue );
    magma_minproduct_cmfree(&p, queue );
    magma_minproduct_cmfree(&v, queue );
    magma_minproduct_cmfree(&s, queue );
    magma_minproduct_cmfree(&t, queue );

    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_cbicgstab */


