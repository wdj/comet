/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zfgmres.cpp normal z -> c, Sun May  3 11:22:59 2015
*/
#include "common_magma_tally4sparse.h"

#define PRECISION_c

/*
#define  q(i)     (q.dval + (i)*dofs)
#define  z(i)     (z.dval + (i)*dofs)
#define  H(i,j)  H[(i)   + (j)*(1+ldh)]
#define HH(i,j) HH[(i)   + (j)*ldh]
#define dH(i,j) dH[(i)   + (j)*(1+ldh)]
*/

// simulate 2-D arrays at the cost of some arithmetic
#define V(i) (V.dval+(i)*dofs)
#define W(i) (W.dval+(i)*dofs)
//#define Vv(i) (&V.dval[(i)*n])
//#define Wv(i) (&W.dval[(i)*n])
#define H(i,j) (H[(j)*m1+(i)])
#define ABS(x)   ((x)<0 ? (-(x)) : (x))




#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


static void
GeneratePlaneRotation(magma_tally4FloatComplex dx, magma_tally4FloatComplex dy, magma_tally4FloatComplex *cs, magma_tally4FloatComplex *sn)
{
    
    if (dy == MAGMA_tally4_C_ZERO) {
        *cs = MAGMA_tally4_C_ONE;
        *sn = MAGMA_tally4_C_ZERO;
    } else if (MAGMA_tally4_C_ABS((dy)) > MAGMA_tally4_C_ABS((dx))) {
        magma_tally4FloatComplex temp = dx / dy;
        *sn = MAGMA_tally4_C_ONE / magma_tally4_csqrt( ( MAGMA_tally4_C_ONE + temp*temp)) ;
        *cs = temp * *sn;
    } else {
        magma_tally4FloatComplex temp = dy / dx;
        *cs = MAGMA_tally4_C_ONE / magma_tally4_csqrt( ( MAGMA_tally4_C_ONE + temp*temp )) ;
        *sn = temp * *cs;
    }

  //  real_Double_t rho = sqrt(MAGMA_tally4_C_REAL(MAGMA_tally4_C_CNJG(dx)*dx + MAGMA_tally4_C_CNJG(dy)*dy));
  //  *cs = dx / rho;
  //  *sn = dy / rho;
}

static void ApplyPlaneRotation(magma_tally4FloatComplex *dx, magma_tally4FloatComplex *dy, magma_tally4FloatComplex cs, magma_tally4FloatComplex sn)
{

    magma_tally4FloatComplex temp = *dx;
    *dx =  cs * *dx + sn * *dy;
    *dy = -sn * temp + cs * *dy;

 //   magma_tally4FloatComplex temp  =  MAGMA_tally4_C_CNJG(cs) * *dx +  MAGMA_tally4_C_CNJG(sn) * *dy;
 //   *dy = -(sn) * *dx + cs * *dy;
 //   *dx = temp;
}



/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a complex sparse matrix stored in the GPU memory.
    X and B are complex vectors stored on the GPU memory.
    This is a GPU implementation of the right-preconditioned flexible GMRES.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_c_matrix
                descriptor for matrix A

    @param[in]
    b           magma_tally4_c_matrix
                RHS b vector

    @param[in,out]
    x           magma_tally4_c_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally4_c_solver_par*
                solver parameters

    @param[in]
    precond_par magma_tally4_c_preconditioner*
                preconditioner
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_cgesv
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cfgmres(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, magma_tally4_c_matrix *x,
    magma_tally4_c_solver_par *solver_par,
    magma_tally4_c_preconditioner *precond_par,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_int_t dofs = A.num_rows;

    // prepare solver feedback
    solver_par->solver = Magma_tally4_PGMRES;
    solver_par->numiter = 0;
    solver_par->info = MAGMA_tally4_SUCCESS;
    
    //Chronometry
    real_Double_t tempo1, tempo2;

    magma_tally4_int_t dim = solver_par->restart;
    magma_tally4_int_t m1 = dim+1; // used inside H macro
    magma_tally4_int_t i, j, k;
    magma_tally4FloatComplex beta;
    
    float rel_resid, resid0, r0=0.0, betanom = 0.0, nom;
    
    magma_tally4_c_matrix v_t={Magma_tally4_CSR}, w_t={Magma_tally4_CSR}, t={Magma_tally4_CSR}, t2={Magma_tally4_CSR}, V={Magma_tally4_CSR}, W={Magma_tally4_CSR};
    v_t.memory_location = Magma_tally4_DEV;
    v_t.num_rows = dofs;
    v_t.num_cols = 1;
    v_t.dval = NULL;
    v_t.storage_type = Magma_tally4_DENSE;

    w_t.memory_location = Magma_tally4_DEV;
    w_t.num_rows = dofs;
    w_t.num_cols = 1;
    w_t.dval = NULL;
    w_t.storage_type = Magma_tally4_DENSE;
    
    magma_tally4FloatComplex temp;
    
    magma_tally4FloatComplex *H={0}, *s={0}, *cs={0}, *sn={0};

    CHECK( magma_tally4_cvinit( &t, Magma_tally4_DEV, dofs, 1, MAGMA_tally4_C_ZERO, queue ));
    CHECK( magma_tally4_cvinit( &t2, Magma_tally4_DEV, dofs, 1, MAGMA_tally4_C_ZERO, queue ));
    
    CHECK( magma_tally4_cmalloc_pinned( &H, (dim+1)*dim ));
    CHECK( magma_tally4_cmalloc_pinned( &s,  dim+1 ));
    CHECK( magma_tally4_cmalloc_pinned( &cs, dim ));
    CHECK( magma_tally4_cmalloc_pinned( &sn, dim ));
    
    
    CHECK( magma_tally4_cvinit( &V, Magma_tally4_DEV, dofs*(dim+1), 1, MAGMA_tally4_C_ZERO, queue ));
    CHECK( magma_tally4_cvinit( &W, Magma_tally4_DEV, dofs*dim, 1, MAGMA_tally4_C_ZERO, queue ));
    
    CHECK(  magma_tally4_cresidual( A, b, *x, &nom, queue));

    solver_par->init_res = nom;
    
    if ( ( nom * solver_par->epsilon) < ATOLERANCE )
        r0 = ATOLERANCE;
    
    solver_par->numiter = 0;
    

    tempo1 = magma_tally4_sync_wtime( queue );
    do
    {
        // compute initial residual and its norm
        // A.mult(n, 1, x, n, V(0), n);                        // V(0) = A*x
        CHECK( magma_tally4_c_spmv( MAGMA_tally4_C_ONE, A, *x, MAGMA_tally4_C_ZERO, t, queue ));
        magma_tally4_ccopy( dofs, t.dval, 1, V(0), 1 );
        
        temp = MAGMA_tally4_C_MAKE(-1.0, 0.0);
        magma_tally4_caxpy(dofs,temp, b.dval, 1, V(0), 1);           // V(0) = V(0) - b
        beta = MAGMA_tally4_C_MAKE( magma_tally4_scnrm2( dofs, V(0), 1 ), 0.0); // beta = norm(V(0))
        if (solver_par->numiter == 0){
            solver_par->init_res = MAGMA_tally4_C_REAL( beta );
            resid0 = MAGMA_tally4_C_REAL( beta );
        
            if ( (r0 = resid0 * solver_par->epsilon) < ATOLERANCE )
                r0 = ATOLERANCE;
            if ( resid0 < r0 ) {
                solver_par->final_res = solver_par->init_res;
                solver_par->iter_res = solver_par->init_res;
                goto cleanup;
            }
        }
        if ( solver_par->verbose > 0 ) {
            solver_par->res_vec[0] = resid0;
            solver_par->timing[0] = 0.0;
        }
        temp = -1.0/beta;
        magma_tally4_cscal( dofs, temp, V(0), 1 );                 // V(0) = -V(0)/beta

        // save very first residual norm
        if (solver_par->numiter == 0)
            solver_par->init_res = MAGMA_tally4_C_REAL( beta );

        for (i = 1; i < dim+1; i++)
            s[i] = MAGMA_tally4_C_ZERO;
        s[0] = beta;

        i = -1;
        do
        {
            solver_par->numiter++;
            i++;
            
            // M.apply(n, 1, V(i), n, W(i), n);
            v_t.dval = V(i);
            CHECK( magma_tally4_c_applyprecond_left( A, v_t, &t, precond_par, queue ));
            CHECK( magma_tally4_c_applyprecond_right( A, t, &t2, precond_par, queue ));
            magma_tally4_ccopy( dofs, t2.dval, 1, W(i), 1 );

            // A.mult(n, 1, W(i), n, V(i+1), n);
            w_t.dval = W(i);
            CHECK( magma_tally4_c_spmv( MAGMA_tally4_C_ONE, A, w_t, MAGMA_tally4_C_ZERO, t, queue ));
            magma_tally4_ccopy( dofs, t.dval, 1, V(i+1), 1 );
            
            for (k = 0; k <= i; k++)
            {
                H(k, i) = magma_tally4_cdotc(dofs, V(k), 1, V(i+1), 1);
                temp = -H(k,i);
                // V(i+1) -= H(k, i) * V(k);
                magma_tally4_caxpy(dofs,-H(k,i), V(k), 1, V(i+1), 1);
            }

            H(i+1, i) = MAGMA_tally4_C_MAKE( magma_tally4_scnrm2(dofs, V(i+1), 1), 0. ); // H(i+1,i) = ||r||
            temp = 1.0 / H(i+1, i);
            // V(i+1) = V(i+1) / H(i+1, i)
            magma_tally4_cscal(dofs, temp, V(i+1), 1);    //  (to be fused)
    
            for (k = 0; k < i; k++)
                ApplyPlaneRotation(&H(k,i), &H(k+1,i), cs[k], sn[k]);
          
            GeneratePlaneRotation(H(i,i), H(i+1,i), &cs[i], &sn[i]);
            ApplyPlaneRotation(&H(i,i), &H(i+1,i), cs[i], sn[i]);
            ApplyPlaneRotation(&s[i], &s[i+1], cs[i], sn[i]);
            
            betanom = MAGMA_tally4_C_ABS( s[i+1] );
            rel_resid = betanom / resid0;
            if ( solver_par->verbose > 0 ) {
                tempo2 = magma_tally4_sync_wtime( queue );
                if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                    solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                            = (real_Double_t) betanom;
                    solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                            = (real_Double_t) tempo2-tempo1;
                }
            }
            if (rel_resid <= solver_par->epsilon){
                break;
            }
        }
        while (i+1 < dim && solver_par->numiter+1 <= solver_par->maxiter);

        // solve upper triangular system in place
        for (j = i; j >= 0; j--)
        {
            s[j] /= H(j,j);
            for (k = j-1; k >= 0; k--)
                s[k] -= H(k,j) * s[j];
        }

        // update the solution
        for (j = 0; j <= i; j++)
        {
            // x = x + s[j] * W(j)
            magma_tally4_caxpy(dofs, s[j], W(j), 1, x->dval, 1);
        }
    }
    while (rel_resid > solver_par->epsilon
                && solver_par->numiter+1 <= solver_par->maxiter);

    tempo2 = magma_tally4_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    float residual;
    CHECK( magma_tally4_cresidual( A, b, *x, &residual, queue ));
    solver_par->iter_res = betanom;
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
    // free pinned memory
    magma_tally4_free_pinned(s);
    magma_tally4_free_pinned(cs);
    magma_tally4_free_pinned(sn);
    magma_tally4_free_pinned(H);

    //free DEV memory
    magma_tally4_cmfree( &V, queue);
    magma_tally4_cmfree( &W, queue);
    magma_tally4_cmfree( &t, queue);
    magma_tally4_cmfree( &t2, queue);

    solver_par->info = info;
    return info;

} /* magma_tally4_cfgmres */



