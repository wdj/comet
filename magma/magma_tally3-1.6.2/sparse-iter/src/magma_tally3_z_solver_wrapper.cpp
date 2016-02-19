/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"




/**
    Purpose
    -------

    ALlows the user to choose a solver.

    Arguments
    ---------

    @param[in]
    A           magma_tally3_z_matrix
                sparse matrix A

    @param[in]
    b           magma_tally3_z_matrix
                input vector b

    @param[in]
    x           magma_tally3_z_matrix*
                output vector x

    @param[in]
    zopts     magma_tally3_zopts
              options for solver and preconditioner
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_z_solver(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b,
    magma_tally3_z_matrix *x, magma_tally3_zopts *zopts,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_tally3_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
    if( b.num_cols == 1 ){
    // preconditioner
        if ( zopts->solver_par.solver != Magma_tally3_ITERREF ) {
            int stat = magma_tally3_z_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally3_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally3_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally3_CG:
                    CHECK( magma_tally3_zcg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_CGMERGE:
                    CHECK( magma_tally3_zcg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_PCG:
                    CHECK( magma_tally3_zpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_BICGSTAB:
                    CHECK( magma_tally3_zbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_BICGSTABMERGE: 
                    CHECK( magma_tally3_zbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_PBICGSTAB:
                    CHECK( magma_tally3_zpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_GMRES:
                    CHECK( magma_tally3_zfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_PGMRES:
                    CHECK( magma_tally3_zfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_LOBPCG:
                    CHECK( magma_tally3_zlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_ITERREF:
                    CHECK( magma_tally3_ziterref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_JACOBI:
                    CHECK( magma_tally3_zjacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_BAITER:
                    CHECK( magma_tally3_zbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_tally3_ITERREF ) {
            int stat = magma_tally3_z_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally3_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally3_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally3_CG:
                    CHECK( magma_tally3_zbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_PCG:
                    CHECK( magma_tally3_zbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_LOBPCG:
                    CHECK( magma_tally3_zlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


