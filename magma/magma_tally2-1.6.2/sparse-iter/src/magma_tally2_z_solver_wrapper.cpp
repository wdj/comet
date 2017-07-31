/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"




/**
    Purpose
    -------

    ALlows the user to choose a solver.

    Arguments
    ---------

    @param[in]
    A           magma_tally2_z_matrix
                sparse matrix A

    @param[in]
    b           magma_tally2_z_matrix
                input vector b

    @param[in]
    x           magma_tally2_z_matrix*
                output vector x

    @param[in]
    zopts     magma_tally2_zopts
              options for solver and preconditioner
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zaux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_z_solver(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b,
    magma_tally2_z_matrix *x, magma_tally2_zopts *zopts,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_tally2_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_tally2_ERR_NOT_SUPPORTED;
    }
    if( b.num_cols == 1 ){
    // preconditioner
        if ( zopts->solver_par.solver != Magma_tally2_ITERREF ) {
            int stat = magma_tally2_z_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally2_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally2_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally2_CG:
                    CHECK( magma_tally2_zcg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_CGMERGE:
                    CHECK( magma_tally2_zcg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_PCG:
                    CHECK( magma_tally2_zpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_BICGSTAB:
                    CHECK( magma_tally2_zbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_BICGSTABMERGE: 
                    CHECK( magma_tally2_zbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_PBICGSTAB:
                    CHECK( magma_tally2_zpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_GMRES:
                    CHECK( magma_tally2_zfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_PGMRES:
                    CHECK( magma_tally2_zfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_LOBPCG:
                    CHECK( magma_tally2_zlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_ITERREF:
                    CHECK( magma_tally2_ziterref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_JACOBI:
                    CHECK( magma_tally2_zjacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_BAITER:
                    CHECK( magma_tally2_zbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_tally2_ITERREF ) {
            int stat = magma_tally2_z_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally2_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally2_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally2_CG:
                    CHECK( magma_tally2_zbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_PCG:
                    CHECK( magma_tally2_zbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_LOBPCG:
                    CHECK( magma_tally2_zlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


