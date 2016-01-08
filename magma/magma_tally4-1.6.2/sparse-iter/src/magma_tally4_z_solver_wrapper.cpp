/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s
       @author Hartwig Anzt

*/
#include "common_magma_tally4sparse.h"




/**
    Purpose
    -------

    ALlows the user to choose a solver.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_z_matrix
                sparse matrix A

    @param[in]
    b           magma_tally4_z_matrix
                input vector b

    @param[in]
    x           magma_tally4_z_matrix*
                output vector x

    @param[in]
    zopts     magma_tally4_zopts
              options for solver and preconditioner
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_z_solver(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b,
    magma_tally4_z_matrix *x, magma_tally4_zopts *zopts,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_tally4_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_tally4_ERR_NOT_SUPPORTED;
    }
    if( b.num_cols == 1 ){
    // preconditioner
        if ( zopts->solver_par.solver != Magma_tally4_ITERREF ) {
            int stat = magma_tally4_z_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally4_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally4_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally4_CG:
                    CHECK( magma_tally4_zcg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_CGMERGE:
                    CHECK( magma_tally4_zcg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_PCG:
                    CHECK( magma_tally4_zpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_BICGSTAB:
                    CHECK( magma_tally4_zbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_BICGSTABMERGE: 
                    CHECK( magma_tally4_zbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_PBICGSTAB:
                    CHECK( magma_tally4_zpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_GMRES:
                    CHECK( magma_tally4_zfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_PGMRES:
                    CHECK( magma_tally4_zfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_LOBPCG:
                    CHECK( magma_tally4_zlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_ITERREF:
                    CHECK( magma_tally4_ziterref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_JACOBI:
                    CHECK( magma_tally4_zjacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_BAITER:
                    CHECK( magma_tally4_zbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_tally4_ITERREF ) {
            int stat = magma_tally4_z_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally4_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally4_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally4_CG:
                    CHECK( magma_tally4_zbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_PCG:
                    CHECK( magma_tally4_zbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_LOBPCG:
                    CHECK( magma_tally4_zlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


