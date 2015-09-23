/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @precisions normal z -> s d c
*/

#include "common_magma_minproductsparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a complex Hermitian N-by-N positive definite matrix A.
    This is a GPU implementation of the Conjugate Gradient method in variant,
    where multiple operations are merged into one compute kernel.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_z_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_z_matrix
                RHS b

    @param[in,out]
    x           magma_minproduct_z_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_minproduct_z_solver_par*
                solver parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zposv
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zcg_merge(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, magma_minproduct_z_matrix *x,
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // prepare solver feedback
    solver_par->solver = Magma_minproduct_CGMERGE;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_minproduct_SUCCESS;
    
    // solver variables
    magma_minproductDoubleComplex alpha, beta, gamma, rho, tmp1, *skp_h={0};
    double nom, nom0, r0, betanom, den;

    // some useful variables
    magma_minproductDoubleComplex c_zero = MAGMA_minproduct_Z_ZERO, c_one = MAGMA_minproduct_Z_ONE;
    magma_minproduct_int_t dofs = A.num_rows*b.num_cols;

    magma_minproduct_z_matrix r={Magma_minproduct_CSR}, d={Magma_minproduct_CSR}, z={Magma_minproduct_CSR};
    magma_minproductDoubleComplex *d1=NULL, *d2=NULL, *skp=NULL;

    // GPU stream
    magma_minproduct_queue_t stream[2]={0};
    magma_minproduct_event_t event[1]={0};
    magma_minproduct_queue_create( &stream[0] );
    magma_minproduct_queue_create( &stream[1] );
    magma_minproduct_event_create( &event[0] );

    // GPU workspace
    CHECK( magma_minproduct_zvinit( &r, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_zvinit( &d, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_minproduct_zvinit( &z, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    
    CHECK( magma_minproduct_zmalloc( &d1, dofs*(1) ));
    CHECK( magma_minproduct_zmalloc( &d2, dofs*(1) ));
    // array for the parameters
    CHECK( magma_minproduct_zmalloc( &skp, 6 ));
    // skp = [alpha|beta|gamma|rho|tmp1|tmp2]

    // solver setup
    magma_minproduct_zscal( dofs, c_zero, x->dval, 1) ;                     // x = 0
    //CHECK(  magma_minproduct_zresidualvec( A, b, *x, &r, nom0, queue));
    magma_minproduct_zcopy( dofs, b.dval, 1, r.dval, 1 );                    // r = b
    magma_minproduct_zcopy( dofs, r.dval, 1, d.dval, 1 );                    // d = r
    nom0 = betanom = magma_minproduct_dznrm2( dofs, r.dval, 1 );
    nom = nom0 * nom0;                                           // nom = r' * r
    CHECK( magma_minproduct_z_spmv( c_one, A, d, c_zero, z, queue ));              // z = A d
    den = MAGMA_minproduct_Z_REAL( magma_minproduct_zdotc(dofs, d.dval, 1, z.dval, 1) ); // den = d'* z
    solver_par->init_res = nom0;
    
    // array on host for the parameters
    CHECK( magma_minproduct_zmalloc_cpu( &skp_h, 6 ));
    
    alpha = rho = gamma = tmp1 = c_one;
    beta =  magma_minproduct_zdotc(dofs, r.dval, 1, r.dval, 1);
    skp_h[0]=alpha;
    skp_h[1]=beta;
    skp_h[2]=gamma;
    skp_h[3]=rho;
    skp_h[4]=tmp1;
    skp_h[5]=MAGMA_minproduct_Z_MAKE(nom, 0.0);

    magma_minproduct_zsetvector( 6, skp_h, 1, skp, 1 );
    
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
        solver_par->res_vec[0] = (real_Double_t) nom0;
        solver_par->timing[0] = 0.0;
    }
    
    solver_par->numiter = 0;
    // start iteration
    do
    {
        solver_par->numiter++;

        magma_minproductblasSetKernelStream(stream[0]);
        
        // computes SpMV and dot product
        CHECK( magma_minproduct_zcgmerge_spmv1(  A, d1, d2, d.dval, z.dval, skp, queue ));
            
        // updates x, r, computes scalars and updates d
        CHECK( magma_minproduct_zcgmerge_xrbeta( dofs, d1, d2, x->dval, r.dval, d.dval, z.dval, skp, queue ));

        // check stopping criterion (asynchronous copy)
        magma_minproduct_zgetvector_async( 1 , skp+1, 1,
                                                    skp_h+1, 1, stream[1] );
        betanom = sqrt(MAGMA_minproduct_Z_REAL(skp_h[1]));

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_minproduct_sync_wtime( queue );
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
    
    tempo2 = magma_minproduct_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    double residual;
    CHECK(  magma_minproduct_zresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->iter_res = betanom;
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
        solver_par->info = MAGMA_minproduct_DIVERGENCE;
    }
    
cleanup:
    magma_minproduct_zmfree(&r, queue );
    magma_minproduct_zmfree(&z, queue );
    magma_minproduct_zmfree(&d, queue );

    magma_minproduct_free( d1 );
    magma_minproduct_free( d2 );
    magma_minproduct_free( skp );
    magma_minproduct_free_cpu( skp_h );

    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_zcg_merge */


