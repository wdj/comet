/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from ziterref.cpp normal z -> s, Sun May  3 11:22:59 2015
*/

#include "common_magma_minproductsparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


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

    @param[in,out]
    precond_par magma_minproduct_s_preconditioner*
                inner solver
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_sgesv
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_siterref(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, magma_minproduct_s_matrix *x,
    magma_minproduct_s_solver_par *solver_par, magma_minproduct_s_preconditioner *precond_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );
    
    // some useful variables
    float c_zero = MAGMA_minproduct_S_ZERO, c_one = MAGMA_minproduct_S_ONE,
                                                c_mone = MAGMA_minproduct_S_NEG_ONE;

    // prepare solver feedback
    solver_par->solver = Magma_minproduct_ITERREF;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_minproduct_SUCCESS;
    
    magma_minproduct_int_t dofs = A.num_rows*b.num_cols;

    // solver variables
    float nom, nom0, r0;
    
    // workspace
    magma_minproduct_s_matrix r={Magma_minproduct_CSR}, z={Magma_minproduct_CSR};
    CHECK( magma_minproduct_svinit( &r, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_svinit( &z, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));

    float residual;
    CHECK( magma_minproduct_sresidual( A, b, *x, &residual, queue ));
    solver_par->init_res = residual;
   

    // solver setup
    magma_minproduct_sscal( dofs, c_zero, x->dval, 1) ;                    // x = 0
    //CHECK(  magma_minproduct_sresidualvec( A, b, *x, &r, nom, queue));
    magma_minproduct_scopy( dofs, b.dval, 1, r.dval, 1 );                    // r = b
    nom0 = magma_minproduct_snrm2(dofs, r.dval, 1);                       // nom0 = || r ||
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
    tempo1 = magma_minproduct_sync_wtime( queue );
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = nom0;
        solver_par->timing[0] = 0.0;
    }
    
    // start iteration
    for( solver_par->numiter= 1; solver_par->numiter<solver_par->maxiter;
                                                    solver_par->numiter++ ) {

        magma_minproduct_sscal( dofs, MAGMA_minproduct_S_MAKE(1./nom, 0.), r.dval, 1) ;  // scale it
        CHECK( magma_minproduct_s_precond( A, r, &z, precond_par, queue )); // inner solver:  A * z = r
        magma_minproduct_sscal( dofs, MAGMA_minproduct_S_MAKE(nom, 0.), z.dval, 1) ;  // scale it
        magma_minproduct_saxpy(dofs,  c_one, z.dval, 1, x->dval, 1);        // x = x + z
        CHECK( magma_minproduct_s_spmv( c_mone, A, *x, c_zero, r, queue ));      // r = - A x
        magma_minproduct_saxpy(dofs,  c_one, b.dval, 1, r.dval, 1);         // r = r + b
        nom = magma_minproduct_snrm2(dofs, r.dval, 1);                    // nom = || r ||

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_minproduct_sync_wtime( queue );
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
    tempo2 = magma_minproduct_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_minproduct_sresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->final_res = residual;
    solver_par->iter_res = nom;

    if ( solver_par->numiter < solver_par->maxiter ) {
        solver_par->info = MAGMA_minproduct_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) nom;
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
                        = (real_Double_t) nom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_minproduct_DIVERGENCE;
    }
    
cleanup:
    magma_minproduct_smfree(&r, queue );
    magma_minproduct_smfree(&z, queue );


    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_siterref */


