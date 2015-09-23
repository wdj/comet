/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_z_solver_wrapper.cpp normal z -> c, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"




/**
    Purpose
    -------

    ALlows the user to choose a solver.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_c_matrix
                sparse matrix A

    @param[in]
    b           magma_minproduct_c_matrix
                input vector b

    @param[in]
    x           magma_minproduct_c_matrix*
                output vector x

    @param[in]
    zopts     magma_minproduct_copts
              options for solver and preconditioner
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_c_solver(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b,
    magma_minproduct_c_matrix *x, magma_minproduct_copts *zopts,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_minproduct_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
    if( b.num_cols == 1 ){
    // preconditioner
        if ( zopts->solver_par.solver != Magma_minproduct_ITERREF ) {
            int stat = magma_minproduct_c_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_minproduct_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_minproduct_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_minproduct_CG:
                    CHECK( magma_minproduct_ccg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_CGMERGE:
                    CHECK( magma_minproduct_ccg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_PCG:
                    CHECK( magma_minproduct_cpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_BICGSTAB:
                    CHECK( magma_minproduct_cbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_BICGSTABMERGE: 
                    CHECK( magma_minproduct_cbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_PBICGSTAB:
                    CHECK( magma_minproduct_cpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_GMRES:
                    CHECK( magma_minproduct_cfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_PGMRES:
                    CHECK( magma_minproduct_cfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_LOBPCG:
                    CHECK( magma_minproduct_clobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_ITERREF:
                    CHECK( magma_minproduct_citerref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_JACOBI:
                    CHECK( magma_minproduct_cjacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_BAITER:
                    CHECK( magma_minproduct_cbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_minproduct_ITERREF ) {
            int stat = magma_minproduct_c_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_minproduct_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_minproduct_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_minproduct_CG:
                    CHECK( magma_minproduct_cbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_PCG:
                    CHECK( magma_minproduct_cbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_LOBPCG:
                    CHECK( magma_minproduct_clobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


