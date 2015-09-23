/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s
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
    A           magma_minproduct_z_matrix
                sparse matrix A

    @param[in]
    b           magma_minproduct_z_matrix
                input vector b

    @param[in]
    x           magma_minproduct_z_matrix*
                output vector x

    @param[in]
    zopts     magma_minproduct_zopts
              options for solver and preconditioner
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_z_solver(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b,
    magma_minproduct_z_matrix *x, magma_minproduct_zopts *zopts,
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
            int stat = magma_minproduct_z_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_minproduct_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_minproduct_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_minproduct_CG:
                    CHECK( magma_minproduct_zcg_res( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_CGMERGE:
                    CHECK( magma_minproduct_zcg_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_PCG:
                    CHECK( magma_minproduct_zpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_BICGSTAB:
                    CHECK( magma_minproduct_zbicgstab( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_BICGSTABMERGE: 
                    CHECK( magma_minproduct_zbicgstab_merge( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_PBICGSTAB:
                    CHECK( magma_minproduct_zpbicgstab( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_GMRES:
                    CHECK( magma_minproduct_zfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_PGMRES:
                    CHECK( magma_minproduct_zfgmres( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_LOBPCG:
                    CHECK( magma_minproduct_zlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_ITERREF:
                    CHECK( magma_minproduct_ziterref( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_JACOBI:
                    CHECK( magma_minproduct_zjacobi( A, b, x, &zopts->solver_par, queue )); break;
            case  Magma_minproduct_BAITER:
                    CHECK( magma_minproduct_zbaiter( A, b, x, &zopts->solver_par, queue ) ); break;
            default:
                    printf("error: solver class not supported.\n"); break;
        }
    }
    else{
  // preconditioner
        if ( zopts->solver_par.solver != Magma_minproduct_ITERREF ) {
            int stat = magma_minproduct_z_precondsetup( A, b, &zopts->precond_par, queue );
            if (  stat != MAGMA_minproduct_SUCCESS ){
                printf("error: bad preconditioner.\n");
                return MAGMA_minproduct_ERR_BADPRECOND; 
            }
        }
        switch( zopts->solver_par.solver ) {
            case  Magma_minproduct_CG:
                    CHECK( magma_minproduct_zbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_PCG:
                    CHECK( magma_minproduct_zbpcg( A, b, x, &zopts->solver_par, &zopts->precond_par, queue )); break;
            case  Magma_minproduct_LOBPCG:
                    CHECK( magma_minproduct_zlobpcg( A, &zopts->solver_par, &zopts->precond_par, queue )); break;
            default:
                    printf("error: only 1 RHS supported for this solver class.\n"); break;
        }
    }
cleanup:
    return info; 
}


