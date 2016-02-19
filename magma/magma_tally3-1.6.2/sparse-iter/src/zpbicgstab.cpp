/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"


#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a complex N-by-N general matrix.
    This is a GPU implementation of the preconditioned
    Biconjugate Gradient Stabelized method.

    Arguments
    ---------

    @param[in]
    A           magma_tally3_z_matrix
                input matrix A

    @param[in]
    b           magma_tally3_z_matrix
                RHS b

    @param[in,out]
    x           magma_tally3_z_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally3_z_solver_par*
                solver parameters

    @param[in]
    precond_par magma_tally3_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zgesv
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_zpbicgstab(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, magma_tally3_z_matrix *x,
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_tally3_PBICGSTAB;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_tally3_SUCCESS;

    // some useful variables
    magma_tally3DoubleComplex c_zero = MAGMA_tally3_Z_ZERO, c_one = MAGMA_tally3_Z_ONE,
                                            c_mone = MAGMA_tally3_Z_NEG_ONE;
    
    magma_tally3_int_t dofs = A.num_rows*b.num_cols;

    // workspace
    magma_tally3_z_matrix r={Magma_tally3_CSR}, rr={Magma_tally3_CSR}, p={Magma_tally3_CSR}, v={Magma_tally3_CSR}, s={Magma_tally3_CSR}, t={Magma_tally3_CSR}, ms={Magma_tally3_CSR}, mt={Magma_tally3_CSR}, y={Magma_tally3_CSR}, z={Magma_tally3_CSR};
    CHECK( magma_tally3_zvinit( &r, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &rr,Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &p, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &v, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &s, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &t, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &ms,Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &mt,Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &y, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_zvinit( &z, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));

    
    // solver variables
    magma_tally3DoubleComplex alpha, beta, omega, rho_old, rho_new;
    double nom, betanom, nom0, r0, den, res;

    // solver setup
    CHECK(  magma_tally3_zresidualvec( A, b, *x, &r, &nom0, queue));
    magma_tally3_zcopy( dofs, r.dval, 1, rr.dval, 1 );                  // rr = r
    betanom = nom0;
    nom = nom0*nom0;
    rho_new = omega = alpha = MAGMA_tally3_Z_MAKE( 1.0, 0. );
    solver_par->init_res = nom0;

    CHECK( magma_tally3_z_spmv( c_one, A, r, c_zero, v, queue ));              // z = A r
    den = MAGMA_tally3_Z_REAL( magma_tally3_zdotc(dofs, v.dval, 1, r.dval, 1) ); // den = z' * r

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
        rho_old = rho_new;                                   // rho_old=rho
        rho_new = magma_tally3_zdotc( dofs, rr.dval, 1, r.dval, 1 );  // rho=<rr,r>
        beta = rho_new/rho_old * alpha/omega;   // beta=rho/rho_old *alpha/omega
        magma_tally3_zscal( dofs, beta, p.dval, 1 );                 // p = beta*p
        magma_tally3_zaxpy( dofs, c_mone * omega * beta, v.dval, 1 , p.dval, 1 );
                                                        // p = p-omega*beta*v
        magma_tally3_zaxpy( dofs, c_one, r.dval, 1, p.dval, 1 );      // p = p+r

        // preconditioner
        CHECK( magma_tally3_z_applyprecond_left( A, p, &mt, precond_par, queue ));
        CHECK( magma_tally3_z_applyprecond_right( A, mt, &y, precond_par, queue ));

        CHECK( magma_tally3_z_spmv( c_one, A, y, c_zero, v, queue ));      // v = Ap

        alpha = rho_new / magma_tally3_zdotc( dofs, rr.dval, 1, v.dval, 1 );
        magma_tally3_zcopy( dofs, r.dval, 1 , s.dval, 1 );            // s=r
        magma_tally3_zaxpy( dofs, c_mone * alpha, v.dval, 1 , s.dval, 1 ); // s=s-alpha*v

        // preconditioner
        CHECK( magma_tally3_z_applyprecond_left( A, s, &ms, precond_par, queue ));
        CHECK( magma_tally3_z_applyprecond_right( A, ms, &z, precond_par, queue ));

        CHECK( magma_tally3_z_spmv( c_one, A, z, c_zero, t, queue ));       // t=As

        // preconditioner
        CHECK( magma_tally3_z_applyprecond_left( A, s, &ms, precond_par, queue ));
        CHECK( magma_tally3_z_applyprecond_left( A, t, &mt, precond_par, queue ));

        // omega = <ms,mt>/<mt,mt>
        omega = magma_tally3_zdotc( dofs, mt.dval, 1, ms.dval, 1 )
                   / magma_tally3_zdotc( dofs, mt.dval, 1, mt.dval, 1 );

        magma_tally3_zaxpy( dofs, alpha, y.dval, 1 , x->dval, 1 );     // x=x+alpha*p
        magma_tally3_zaxpy( dofs, omega, z.dval, 1 , x->dval, 1 );     // x=x+omega*s

        magma_tally3_zcopy( dofs, s.dval, 1 , r.dval, 1 );             // r=s
        magma_tally3_zaxpy( dofs, c_mone * omega, t.dval, 1 , r.dval, 1 ); // r=r-omega*t
        res = betanom = magma_tally3_dznrm2( dofs, r.dval, 1 );

        nom = betanom*betanom;


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
    double residual;
    CHECK(  magma_tally3_zresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->final_res = residual;
    solver_par->iter_res = res;

    if ( solver_par->numiter < solver_par->maxiter ) {
        info = MAGMA_tally3_SUCCESS;
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
    magma_tally3_zmfree(&r, queue );
    magma_tally3_zmfree(&rr, queue );
    magma_tally3_zmfree(&p, queue );
    magma_tally3_zmfree(&v, queue );
    magma_tally3_zmfree(&s, queue );
    magma_tally3_zmfree(&t, queue );
    magma_tally3_zmfree(&ms, queue );
    magma_tally3_zmfree(&mt, queue );
    magma_tally3_zmfree(&y, queue );
    magma_tally3_zmfree(&z, queue );

    magma_tally3blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally3_zbicgstab */


