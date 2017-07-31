/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_z_solver_wrapper.cpp normal z -> d, Sun May  3 11:22:59 2015
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
    A           magma_tally2_d_matrix
                sparse matrix A

    @param[in]
    b           magma_tally2_d_matrix
                input vector b

    @param[in]
    x           magma_tally2_d_matrix*
                output vector x

    @param[in]
    zopts     magma_tally2_dopts
              options for solver and preconditioner
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_daux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_d_solver(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b,
    magma_tally2_d_matrix *x, magma_tally2_dopts *zopts,
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
            int stat = magma_tally2_d_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally2_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally2_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally2_CG:
                    CHECK( magma_tally2_dcg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_CGMERGE:
                    CHECK( magma_tally2_dcg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_PCG:
                    CHECK( magma_tally2_dpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_BICGSTAB:
                    CHECK( magma_tally2_dbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_BICGSTABMERGE: 
                    CHECK( magma_tally2_dbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_PBICGSTAB:
                    CHECK( magma_tally2_dpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_GMRES:
                    CHECK( magma_tally2_dfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_PGMRES:
                    CHECK( magma_tally2_dfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_LOBPCG:
                    CHECK( magma_tally2_dlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_ITERREF:
                    CHECK( magma_tally2_diterref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_JACOBI:
                    CHECK( magma_tally2_djacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally2_BAITER:
                    CHECK( magma_tally2_dbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_tally2_ITERREF ) {
            int stat = magma_tally2_d_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally2_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally2_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally2_CG:
                    CHECK( magma_tally2_dbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_PCG:
                    CHECK( magma_tally2_dbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally2_LOBPCG:
                    CHECK( magma_tally2_dlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


