/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zpcg.cpp normal z -> d, Sun May  3 11:22:59 2015
*/

#include "common_magma_minproductsparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a real symmetric N-by-N positive definite matrix A.
    This is a GPU implementation of the preconditioned Conjugate
    Gradient method.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                RHS b

    @param[in,out]
    x           magma_minproduct_d_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_minproduct_d_solver_par*
                solver parameters

    @param[in]
    precond_par magma_minproduct_d_preconditioner*
                preconditioner
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_dposv
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_dpcg(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_minproduct_PCG;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_minproduct_SUCCESS;
    
    // solver variables
    double alpha, beta;
    double nom, nom0, r0, gammaold, gammanew, den, res;

    // local variables
    double c_zero = MAGMA_minproduct_D_ZERO, c_one = MAGMA_minproduct_D_ONE;
    
    magma_minproduct_int_t dofs = A.num_rows* b.num_cols;

    // GPU workspace
    magma_minproduct_d_matrix r={Magma_minproduct_CSR}, rt={Magma_minproduct_CSR}, p={Magma_minproduct_CSR}, q={Magma_minproduct_CSR}, h={Magma_minproduct_CSR};
    CHECK( magma_minproduct_dvinit( &r, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_dvinit( &rt,Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_dvinit( &p, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_dvinit( &q, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_dvinit( &h, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    

    // solver setup
    CHECK(  magma_minproduct_dresidualvec( A, b, *x, &r, &nom0, queue));

    // preconditioner
    CHECK( magma_minproduct_d_applyprecond_left( A, r, &rt, precond_par, queue ));
    CHECK( magma_minproduct_d_applyprecond_right( A, rt, &h, precond_par, queue ));

    magma_minproduct_dcopy( dofs, h.dval, 1, p.dval, 1 );                    // p = h
    nom = MAGMA_minproduct_D_REAL( magma_minproduct_ddot(dofs, r.dval, 1, h.dval, 1) );
    CHECK( magma_minproduct_d_spmv( c_one, A, p, c_zero, q, queue ));             // q = A p
    den = MAGMA_minproduct_D_REAL( magma_minproduct_ddot(dofs, p.dval, 1, q.dval, 1) );// den = p dot q
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
        info = MAGMA_minproduct_NONSPD;
        goto cleanup;
    }

    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_minproduct_sync_wtime( queue );
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = (real_Double_t)nom0;
        solver_par->timing[0] = 0.0;
    }
    
    solver_par->numiter = 0;
    // start iteration
    do
    {
        solver_par->numiter++;

        // preconditioner
        CHECK( magma_minproduct_d_applyprecond_left( A, r, &rt, precond_par, queue ));
        CHECK( magma_minproduct_d_applyprecond_right( A, rt, &h, precond_par, queue ));

        gammanew = MAGMA_minproduct_D_REAL( magma_minproduct_ddot(dofs, r.dval, 1, h.dval, 1) );
                                                            // gn = < r,h>

        if ( solver_par->numiter==1 ) {
            magma_minproduct_dcopy( dofs, h.dval, 1, p.dval, 1 );                    // p = h
        } else {
            beta = MAGMA_minproduct_D_MAKE(gammanew/gammaold, 0.);       // beta = gn/go
            magma_minproduct_dscal(dofs, beta, p.dval, 1);            // p = beta*p
            magma_minproduct_daxpy(dofs, c_one, h.dval, 1, p.dval, 1); // p = p + h
        }

        CHECK( magma_minproduct_d_spmv( c_one, A, p, c_zero, q, queue ));   // q = A p
        den = MAGMA_minproduct_D_REAL(magma_minproduct_ddot(dofs, p.dval, 1, q.dval, 1));
                // den = p dot q

        alpha = MAGMA_minproduct_D_MAKE(gammanew/den, 0.);
        magma_minproduct_daxpy(dofs,  alpha, p.dval, 1, x->dval, 1);     // x = x + alpha p
        magma_minproduct_daxpy(dofs, -alpha, q.dval, 1, r.dval, 1);      // r = r - alpha q
        gammaold = gammanew;

        res = magma_minproduct_dnrm2( dofs, r.dval, 1 );
        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_minproduct_sync_wtime( queue );
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }


        if (  res/nom0  < solver_par->epsilon ) {
            break;
        }
    }
    while ( solver_par->numiter+1 <= solver_par->maxiter );
    
    tempo2 = magma_minproduct_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    double residual;
    CHECK(  magma_minproduct_dresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->iter_res = res;
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter ) {
        info = MAGMA_minproduct_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
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
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_minproduct_DIVERGENCE;
    }
    
cleanup:
    magma_minproduct_dmfree(&r, queue );
    magma_minproduct_dmfree(&rt, queue );
    magma_minproduct_dmfree(&p, queue );
    magma_minproduct_dmfree(&q, queue );
    magma_minproduct_dmfree(&h, queue );

    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_dcg */


