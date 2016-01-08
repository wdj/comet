/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_z_solver_wrapper.cpp normal z -> c, Sun May  3 11:22:59 2015
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
    A           magma_tally4_c_matrix
                sparse matrix A

    @param[in]
    b           magma_tally4_c_matrix
                input vector b

    @param[in]
    x           magma_tally4_c_matrix*
                output vector x

    @param[in]
    zopts     magma_tally4_copts
              options for solver and preconditioner
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_caux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_c_solver(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b,
    magma_tally4_c_matrix *x, magma_tally4_copts *zopts,
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
            int stat = magma_tally4_c_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally4_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally4_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally4_CG:
                    CHECK( magma_tally4_ccg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_CGMERGE:
                    CHECK( magma_tally4_ccg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_PCG:
                    CHECK( magma_tally4_cpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_BICGSTAB:
                    CHECK( magma_tally4_cbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_BICGSTABMERGE: 
                    CHECK( magma_tally4_cbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_PBICGSTAB:
                    CHECK( magma_tally4_cpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_GMRES:
                    CHECK( magma_tally4_cfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_PGMRES:
                    CHECK( magma_tally4_cfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_LOBPCG:
                    CHECK( magma_tally4_clobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_ITERREF:
                    CHECK( magma_tally4_citerref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_JACOBI:
                    CHECK( magma_tally4_cjacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally4_BAITER:
                    CHECK( magma_tally4_cbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_tally4_ITERREF ) {
            int stat = magma_tally4_c_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally4_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally4_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally4_CG:
                    CHECK( magma_tally4_cbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_PCG:
                    CHECK( magma_tally4_cbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally4_LOBPCG:
                    CHECK( magma_tally4_clobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


