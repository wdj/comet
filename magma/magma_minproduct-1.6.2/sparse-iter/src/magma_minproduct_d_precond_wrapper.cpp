/*
    -- MAGMA_minproduct (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @generated from magma_minproduct_z_precond_wrapper.cpp normal z -> d, Mon May  4 11:57:23 2015
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"




/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and the
    preconditioner parameters, the respective preconditioner
    is chosen. It approximates x for A x = y.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                sparse matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                input vector b

    @param[in]
    x           magma_minproduct_d_matrix*
                output vector x

    @param[in,out]
    precond     magma_minproduct_d_preconditioner
                preconditioner

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_d_precond(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set up precond parameters as solver parameters
    magma_minproduct_d_solver_par psolver_par;
    psolver_par.epsilon = precond->epsilon;
    psolver_par.maxiter = precond->maxiter;
    psolver_par.restart = precond->restart;
    psolver_par.verbose = 0;
    magma_minproduct_d_preconditioner pprecond;
    pprecond.solver = Magma_minproduct_NONE;

    switch( precond->solver ) {
        case  Magma_minproduct_CG:
                CHECK( magma_minproduct_dcg_res( A, b, x, &psolver_par, queue )); break;
        case  Magma_minproduct_BICGSTAB:
                CHECK( magma_minproduct_dbicgstab( A, b, x, &psolver_par, queue )); break;
        case  Magma_minproduct_GMRES:
                CHECK( magma_minproduct_dfgmres( A, b, x, &psolver_par, &pprecond, queue )); break;
        case  Magma_minproduct_JACOBI:
                CHECK( magma_minproduct_djacobi( A, b, x, &psolver_par, queue )); break;
        case  Magma_minproduct_BAITER:
                CHECK( magma_minproduct_dbaiter( A, b, x, &psolver_par, queue )); break;
        default:
                CHECK( magma_minproduct_dcg_res( A, b, x, &psolver_par, queue )); break;

    }
cleanup:
    return info;
}



/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and the
    preconditioner parameters, the respective preconditioner
    is preprocessed.
    E.g. for Jacobi: the scaling-vetor, for ILU the factorization.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                sparse matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                input vector y

    @param[in,out]
    precond     magma_minproduct_d_preconditioner
                preconditioner
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_d_precondsetup(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b,
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_minproduct_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }

    if ( precond->solver == Magma_minproduct_JACOBI ) {
        return magma_minproduct_djacobisetup_diagscal( A, &(precond->d), queue );
    }
    else if ( precond->solver == Magma_minproduct_PASTIX ) {
        //return magma_minproduct_dpastixsetup( A, b, precond, queue );
        return MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
    else if ( precond->solver == Magma_minproduct_ILU ) {
        return magma_minproduct_dcumilusetup( A, precond, queue );
    }
    else if ( precond->solver == Magma_minproduct_ICC ) {
        return magma_minproduct_dcumiccsetup( A, precond, queue );
    }
    else if ( precond->solver == Magma_minproduct_AICC ) {
        return magma_minproduct_ditericsetup( A, b, precond, queue );
    }
    else if ( precond->solver == Magma_minproduct_AILU ) {
        return magma_minproduct_diterilusetup( A, b, precond, queue );
    }
    else if ( precond->solver == Magma_minproduct_NONE ) {
        return MAGMA_minproduct_SUCCESS;
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        return MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
}



/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and the
    preconditioner parameters, the respective preconditioner
    is applied.
    E.g. for Jacobi: the scaling-vetor, for ILU the triangular solves.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                sparse matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                input vector b

    @param[in,out]
    x           magma_minproduct_d_matrix*
                output vector x

    @param[in]
    precond     magma_minproduct_d_preconditioner
                preconditioner

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_d_applyprecond(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_d_matrix tmp={Magma_minproduct_CSR};
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_minproduct_JACOBI ) {
        CHECK( magma_minproduct_djacobi_diagscal( A.num_rows, precond->d, b, x, queue ));
    }
    else if ( precond->solver == Magma_minproduct_PASTIX ) {
        //CHECK( magma_minproduct_dapplypastix( b, x, precond, queue ));
        return MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
    else if ( precond->solver == Magma_minproduct_ILU ) {
        CHECK( magma_minproduct_dvinit( &tmp, Magma_minproduct_DEV, A.num_rows, b.num_cols, MAGMA_minproduct_D_ZERO, queue ));
    }
    else if ( precond->solver == Magma_minproduct_ICC ) {
        CHECK( magma_minproduct_dvinit( &tmp, Magma_minproduct_DEV, A.num_rows, b.num_cols, MAGMA_minproduct_D_ZERO, queue ));
    }
    else if ( precond->solver == Magma_minproduct_NONE ) {
        magma_minproduct_dcopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_minproductblasSetKernelStream( orig_queue );
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_minproduct_dmfree( &tmp, queue );
    magma_minproductblasSetKernelStream( orig_queue );
    return info;
}


/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and the
    preconditioner parameters, the respective left preconditioner
    is applied.
    E.g. for Jacobi: the scaling-vetor, for ILU the left triangular solve.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                sparse matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                input vector b

    @param[in,out]
    x           magma_minproduct_d_matrix*
                output vector x

    @param[in]
    precond     magma_minproduct_d_preconditioner
                preconditioner

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_d_applyprecond_left(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_minproduct_JACOBI ) {
        CHECK( magma_minproduct_djacobi_diagscal( A.num_rows, precond->d, b, x, queue ));
    }
    else if ( ( precond->solver == Magma_minproduct_ILU ||
                precond->solver == Magma_minproduct_AILU ) && precond->maxiter >= 50 ) {
        CHECK( magma_minproduct_dapplycumilu_l( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_minproduct_ILU ||
                precond->solver == Magma_minproduct_AILU ) && precond->maxiter < 50 ) {
        magma_minproduct_dcopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_minproduct_d_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_minproduct_djacobiiter_sys( precond->L, b, precond->d, precond->work1, x, &solver_par, queue );
        CHECK( magma_minproduct_djacobispmvupdate(precond->maxiter, precond->L, precond->work1, b, precond->d, x, queue ));
    }
    else if ( ( precond->solver == Magma_minproduct_ICC ||
                precond->solver == Magma_minproduct_AICC ) && precond->maxiter >= 50 )  {
        CHECK( magma_minproduct_dapplycumicc_l( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_minproduct_ICC ||
                precond->solver == Magma_minproduct_AICC ) && precond->maxiter < 50 )  {
        magma_minproduct_dcopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_minproduct_d_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_minproduct_djacobiiter_sys( precond->L, b, precond->d, precond->work1, x, &solver_par, queue );
        CHECK( magma_minproduct_djacobispmvupdate(precond->maxiter, precond->L, precond->work1, b, precond->d, x, queue ));
    }
    else if ( precond->solver == Magma_minproduct_NONE ) {
        magma_minproduct_dcopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else if ( precond->solver == Magma_minproduct_FUNCTION ) {
        CHECK( magma_minproduct_dapplycustomprecond_l( b, x, precond, queue ));
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_minproductblasSetKernelStream( orig_queue );
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED; 
    }
cleanup:
    magma_minproductblasSetKernelStream( orig_queue );
    return info;
}


/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and the
    preconditioner parameters, the respective right-preconditioner
    is applied.
    E.g. for Jacobi: the scaling-vetor, for ILU the right triangular solve.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                sparse matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                input vector b

    @param[in,out]
    x           magma_minproduct_d_matrix*
                output vector x

    @param[in]
    precond     magma_minproduct_d_preconditioner
                preconditioner

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_d_applyprecond_right(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_minproduct_JACOBI ) {
        magma_minproduct_dcopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );    // x = b
    }
    else if ( ( precond->solver == Magma_minproduct_ILU ||
                precond->solver == Magma_minproduct_AILU ) && precond->maxiter >= 50 ) {
        CHECK( magma_minproduct_dapplycumilu_r( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_minproduct_ILU ||
                precond->solver == Magma_minproduct_AILU ) && precond->maxiter < 50 ) {
        magma_minproduct_dcopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_minproduct_d_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_minproduct_djacobiiter_sys( precond->U, b, precond->d2, precond->work2, x, &solver_par, queue );
        CHECK( magma_minproduct_djacobispmvupdate(precond->maxiter, precond->U, precond->work2, b, precond->d2, x, queue ));
    }

    else if ( ( precond->solver == Magma_minproduct_ICC ||
                precond->solver == Magma_minproduct_AICC ) && precond->maxiter >= 50 ) {
        CHECK( magma_minproduct_dapplycumicc_r( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_minproduct_ICC ||
               precond->solver == Magma_minproduct_AICC ) && precond->maxiter < 50 ) {
        magma_minproduct_dcopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_minproduct_d_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_minproduct_djacobiiter_sys( precond->U, b, precond->d2, precond->work2, x, &solver_par, queue );
        CHECK( magma_minproduct_djacobispmvupdate(precond->maxiter, precond->U, precond->work2, b, precond->d2, x, queue ));
    }
    else if ( precond->solver == Magma_minproduct_NONE ) {
        magma_minproduct_dcopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else if ( precond->solver == Magma_minproduct_FUNCTION ) {
        CHECK( magma_minproduct_dapplycustomprecond_r( b, x, precond, queue ));
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_minproductblasSetKernelStream( orig_queue );
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_minproductblasSetKernelStream( orig_queue );
    return info;
}


