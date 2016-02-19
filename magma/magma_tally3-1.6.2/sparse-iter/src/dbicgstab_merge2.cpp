/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zbicgstab_merge2.cpp normal z -> d, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )

#define  q(i)     (q.dval + (i)*dofs)

/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a general matrix.
    This is a GPU implementation of the Biconjugate Gradient Stabilized method.
    The difference to magma_tally3_dbicgstab is that we use specifically designed kernels
    merging multiple operations into one kernel.

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
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_dgesv
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_dbicgstab_merge2(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b,
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_tally3_BICGSTABMERGE2;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_tally3_SUCCESS;
    
    // solver variables
    double alpha, beta, omega, rho_old, rho_new, *skp_h={0};
    double nom, nom0, betanom, r0, den;

    // some useful variables
    double c_zero = MAGMA_tally3_D_ZERO, c_one = MAGMA_tally3_D_ONE;
    
    magma_tally3_int_t dofs = A.num_rows;

    // GPU stream
    magma_tally3_queue_t stream[2]={0};
    magma_tally3_event_t event[1]={0};
    magma_tally3_queue_create( &stream[0] );
    magma_tally3_queue_create( &stream[1] );
    magma_tally3_event_create( &event[0] );

    // workspace
    magma_tally3_d_matrix q={Magma_tally3_CSR}, r={Magma_tally3_CSR}, rr={Magma_tally3_CSR}, p={Magma_tally3_CSR}, v={Magma_tally3_CSR}, s={Magma_tally3_CSR}, t={Magma_tally3_CSR};
    double *d1=NULL, *d2=NULL, *skp=NULL;
    d1 = NULL;
    d2 = NULL;
    skp = NULL;
    
    CHECK( magma_tally3_dmalloc( &d1, dofs*(2) ));
    CHECK( magma_tally3_dmalloc( &d2, dofs*(2) ));

    // array for the parameters
    CHECK( magma_tally3_dmalloc( &skp, 8 ));
    // skp = [alpha|beta|omega|rho_old|rho|nom|tmp1|tmp2]
    CHECK( magma_tally3_dvinit( &q, Magma_tally3_DEV, dofs*6, 1, c_zero, queue ));

    // q = rr|r|p|v|s|t
    rr.memory_location = Magma_tally3_DEV; rr.dval = NULL; rr.num_rows = rr.nnz = dofs;
    r.memory_location = Magma_tally3_DEV; r.dval = NULL; r.num_rows = r.nnz = dofs;
    p.memory_location = Magma_tally3_DEV; p.dval = NULL; p.num_rows = p.nnz = dofs;
    v.memory_location = Magma_tally3_DEV; v.dval = NULL; v.num_rows = v.nnz = dofs;
    s.memory_location = Magma_tally3_DEV; s.dval = NULL; s.num_rows = s.nnz = dofs;
    t.memory_location = Magma_tally3_DEV; t.dval = NULL; t.num_rows = t.nnz = dofs;

    rr.dval = q(0);
    r.dval = q(1);
    p.dval = q(2);
    v.dval = q(3);
    s.dval = q(4);
    t.dval = q(5);

    // solver setup
    magma_tally3_dscal( dofs, c_zero, x->dval, 1) ;                            // x = 0
    CHECK(  magma_tally3_dresidualvec( A, b, *x, &r, &nom0, queue));
    magma_tally3_dcopy( dofs, r.dval, 1, q(0), 1 );                            // rr = r
    magma_tally3_dcopy( dofs, r.dval, 1, q(1), 1 );                            // q = r
    betanom = nom0;
    nom = nom0*nom0;
    rho_new = magma_tally3_ddot( dofs, r.dval, 1, r.dval, 1 );             // rho=<rr,r>
    rho_old = omega = alpha = MAGMA_tally3_D_MAKE( 1.0, 0. );
    beta = rho_new;
    solver_par->init_res = nom0;
    // array on host for the parameters
    CHECK( magma_tally3_dmalloc_cpu( &skp_h, 8 ));
    
    skp_h[0]=alpha;
    skp_h[1]=beta;
    skp_h[2]=omega;
    skp_h[3]=rho_old;
    skp_h[4]=rho_new;
    skp_h[5]=MAGMA_tally3_D_MAKE(nom, 0.0);
    magma_tally3_dsetvector( 8, skp_h, 1, skp, 1 );

    CHECK( magma_tally3_d_spmv( c_one, A, r, c_zero, v, queue ));             // z = A r
    den = MAGMA_tally3_D_REAL( magma_tally3_ddot(dofs, v.dval, 1, r.dval, 1) );// den = z dot r

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

        magma_tally3blasSetKernelStream(stream[0]);

        // computes p=r+beta*(p-omega*v)
        CHECK( magma_tally3_dbicgmerge1( dofs, skp, v.dval, r.dval, p.dval, queue ));
        CHECK( magma_tally3_dbicgmerge_spmv1(  A, d1, d2, q(2), q(0), q(3), skp, queue ));
        CHECK( magma_tally3_dbicgmerge2( dofs, skp, r.dval, v.dval, s.dval, queue )); // s=r-alpha*v
        CHECK( magma_tally3_dbicgmerge_spmv2( A, d1, d2, q(4), q(5), skp, queue ));
        CHECK( magma_tally3_dbicgmerge_xrbeta( dofs, d1, d2, q(0), q(1), q(2),
                                                    q(4), q(5), x->dval, skp, queue ));

        // check stopping criterion (asynchronous copy)
        magma_tally3_dgetvector_async( 1 , skp+5, 1,
                                                        skp_h+5, 1, stream[1] );

        betanom = sqrt(MAGMA_tally3_D_REAL(skp_h[5]));

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_tally3_sync_wtime( queue );
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
    }
    while ( solver_par->numiter+1 <= solver_par->maxiter );
    
    tempo2 = magma_tally3_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    double residual;
    CHECK( magma_tally3_dresidual( A, b, *x, &residual, queue ));
    solver_par->iter_res = betanom;
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
    magma_tally3_dmfree(&q, queue );  // frees all vectors
    magma_tally3_free(d1);
    magma_tally3_free(d2);
    magma_tally3_free( skp );
    magma_tally3_free_cpu( skp_h );

    magma_tally3blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* dbicgstab_merge2 */


