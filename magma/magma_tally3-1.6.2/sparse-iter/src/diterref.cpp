/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from ziterref.cpp normal z -> d, Sun May  3 11:22:59 2015
*/

#include "common_magma_tally3sparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a real symmetric N-by-N positive definite matrix A.
    This is a GPU implementation of the Iterative Refinement method.
    The inner solver is passed via the preconditioner argument.

    Arguments
    ---------

    @param[in]
    A           magma_tally3_d_matrix
                input matrix A

    @param[in]
    b           magma_tally3_d_matrix
                RHS b

    @param[in,out]
    x           magma_tally3_d_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally3_d_solver_par*
                solver parameters

    @param[in,out]
    precond_par magma_tally3_d_preconditioner*
                inner solver
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_dgesv
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_diterref(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, magma_tally3_d_matrix *x,
    magma_tally3_d_solver_par *solver_par, magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );
    
    // some useful variables
    double c_zero = MAGMA_tally3_D_ZERO, c_one = MAGMA_tally3_D_ONE,
                                                c_mone = MAGMA_tally3_D_NEG_ONE;

    // prepare solver feedback
    solver_par->solver = Magma_tally3_ITERREF;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_tally3_SUCCESS;
    
    magma_tally3_int_t dofs = A.num_rows*b.num_cols;

    // solver variables
    double nom, nom0, r0;
    
    // workspace
    magma_tally3_d_matrix r={Magma_tally3_CSR}, z={Magma_tally3_CSR};
    CHECK( magma_tally3_dvinit( &r, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_tally3_dvinit( &z, Magma_tally3_DEV, A.num_rows, b.num_cols, c_zero, queue ));

    double residual;
    CHECK( magma_tally3_dresidual( A, b, *x, &residual, queue ));
    solver_par->init_res = residual;
   

    // solver setup
    magma_tally3_dscal( dofs, c_zero, x->dval, 1) ;                    // x = 0
    //CHECK(  magma_tally3_dresidualvec( A, b, *x, &r, nom, queue));
    magma_tally3_dcopy( dofs, b.dval, 1, r.dval, 1 );                    // r = b
    nom0 = magma_tally3_dnrm2(dofs, r.dval, 1);                       // nom0 = || r ||
    nom = nom0 * nom0;
    solver_par->init_res = nom0;

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
    
    // start iteration
    for( solver_par->numiter= 1; solver_par->numiter<solver_par->maxiter;
                                                    solver_par->numiter++ ) {

        magma_tally3_dscal( dofs, MAGMA_tally3_D_MAKE(1./nom, 0.), r.dval, 1) ;  // scale it
        CHECK( magma_tally3_d_precond( A, r, &z, precond_par, queue )); // inner solver:  A * z = r
        magma_tally3_dscal( dofs, MAGMA_tally3_D_MAKE(nom, 0.), z.dval, 1) ;  // scale it
        magma_tally3_daxpy(dofs,  c_one, z.dval, 1, x->dval, 1);        // x = x + z
        CHECK( magma_tally3_d_spmv( c_mone, A, *x, c_zero, r, queue ));      // r = - A x
        magma_tally3_daxpy(dofs,  c_one, b.dval, 1, r.dval, 1);         // r = r + b
        nom = magma_tally3_dnrm2(dofs, r.dval, 1);                    // nom = || r ||

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_tally3_sync_wtime( queue );
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) nom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }

        if (  nom  < r0 ) {
            break;
        }
    }
    tempo2 = magma_tally3_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_tally3_dresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->final_res = residual;
    solver_par->iter_res = nom;

    if ( solver_par->numiter < solver_par->maxiter ) {
        solver_par->info = MAGMA_tally3_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) nom;
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
                        = (real_Double_t) nom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_tally3_DIVERGENCE;
    }
    
cleanup:
    magma_tally3_dmfree(&r, queue );
    magma_tally3_dmfree(&z, queue );


    magma_tally3blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally3_diterref */


