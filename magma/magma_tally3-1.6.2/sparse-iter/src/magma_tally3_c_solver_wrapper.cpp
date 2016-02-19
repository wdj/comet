/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_z_solver_wrapper.cpp normal z -> c, Sun May  3 11:22:59 2015
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
    A           magma_tally3_c_matrix
                sparse matrix A

    @param[in]
    b           magma_tally3_c_matrix
                input vector b

    @param[in]
    x           magma_tally3_c_matrix*
                output vector x

    @param[in]
    zopts     magma_tally3_copts
              options for solver and preconditioner
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_solver(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b,
    magma_tally3_c_matrix *x, magma_tally3_copts *zopts,
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
            int stat = magma_tally3_c_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally3_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally3_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally3_CG:
                    CHECK( magma_tally3_ccg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_CGMERGE:
                    CHECK( magma_tally3_ccg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_PCG:
                    CHECK( magma_tally3_cpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_BICGSTAB:
                    CHECK( magma_tally3_cbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_BICGSTABMERGE: 
                    CHECK( magma_tally3_cbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_PBICGSTAB:
                    CHECK( magma_tally3_cpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_GMRES:
                    CHECK( magma_tally3_cfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_PGMRES:
                    CHECK( magma_tally3_cfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_LOBPCG:
                    CHECK( magma_tally3_clobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_ITERREF:
                    CHECK( magma_tally3_citerref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_JACOBI:
                    CHECK( magma_tally3_cjacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_tally3_BAITER:
                    CHECK( magma_tally3_cbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_tally3_ITERREF ) {
            int stat = magma_tally3_c_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_tally3_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_tally3_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_tally3_CG:
                    CHECK( magma_tally3_cbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_PCG:
                    CHECK( magma_tally3_cbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_tally3_LOBPCG:
                    CHECK( magma_tally3_clobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


