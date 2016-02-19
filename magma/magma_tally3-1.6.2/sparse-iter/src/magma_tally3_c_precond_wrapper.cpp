/*
    -- MAGMA_tally3 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @generated from magma_tally3_z_precond_wrapper.cpp normal z -> c, Mon May  4 11:57:23 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"




/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and the
    preconditioner parameters, the respective preconditioner
    is chosen. It approximates x for A x = y.

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

    @param[in,out]
    precond     magma_tally3_c_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_precond(
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix b,
    magma_tally3_c_matrix *x,
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set up precond parameters as solver parameters
    magma_tally3_c_solver_par psolver_par;
    psolver_par.epsilon = precond->epsilon;
    psolver_par.maxiter = precond->maxiter;
    psolver_par.restart = precond->restart;
    psolver_par.verbose = 0;
    magma_tally3_c_preconditioner pprecond;
    pprecond.solver = Magma_tally3_NONE;

    switch( precond->solver ) {
        case  Magma_tally3_CG:
                CHECK( magma_tally3_ccg_res( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally3_BICGSTAB:
                CHECK( magma_tally3_cbicgstab( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally3_GMRES:
                CHECK( magma_tally3_cfgmres( A, b, x, &psolver_par, &pprecond, queue )); break;
        case  Magma_tally3_JACOBI:
                CHECK( magma_tally3_cjacobi( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally3_BAITER:
                CHECK( magma_tally3_cbaiter( A, b, x, &psolver_par, queue )); break;
        default:
                CHECK( magma_tally3_ccg_res( A, b, x, &psolver_par, queue )); break;

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
    A           magma_tally3_c_matrix
                sparse matrix A

    @param[in]
    b           magma_tally3_c_matrix
                input vector y

    @param[in,out]
    precond     magma_tally3_c_preconditioner
                preconditioner
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_precondsetup(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b,
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_tally3_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_tally3_ERR_NOT_SUPPORTED;
    }

    if ( precond->solver == Magma_tally3_JACOBI ) {
        return magma_tally3_cjacobisetup_diagscal( A, &(precond->d), queue );
    }
    else if ( precond->solver == Magma_tally3_PASTIX ) {
        //return magma_tally3_cpastixsetup( A, b, precond, queue );
        return MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
    else if ( precond->solver == Magma_tally3_ILU ) {
        return magma_tally3_ccumilusetup( A, precond, queue );
    }
    else if ( precond->solver == Magma_tally3_ICC ) {
        return magma_tally3_ccumiccsetup( A, precond, queue );
    }
    else if ( precond->solver == Magma_tally3_AICC ) {
        return magma_tally3_citericsetup( A, b, precond, queue );
    }
    else if ( precond->solver == Magma_tally3_AILU ) {
        return magma_tally3_citerilusetup( A, b, precond, queue );
    }
    else if ( precond->solver == Magma_tally3_NONE ) {
        return MAGMA_tally3_SUCCESS;
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        return MAGMA_tally3_ERR_NOT_SUPPORTED;
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
    A           magma_tally3_c_matrix
                sparse matrix A

    @param[in]
    b           magma_tally3_c_matrix
                input vector b

    @param[in,out]
    x           magma_tally3_c_matrix*
                output vector x

    @param[in]
    precond     magma_tally3_c_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_applyprecond(
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix b,
    magma_tally3_c_matrix *x,
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_c_matrix tmp={Magma_tally3_CSR};
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally3_JACOBI ) {
        CHECK( magma_tally3_cjacobi_diagscal( A.num_rows, precond->d, b, x, queue ));
    }
    else if ( precond->solver == Magma_tally3_PASTIX ) {
        //CHECK( magma_tally3_capplypastix( b, x, precond, queue ));
        return MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
    else if ( precond->solver == Magma_tally3_ILU ) {
        CHECK( magma_tally3_cvinit( &tmp, Magma_tally3_DEV, A.num_rows, b.num_cols, MAGMA_tally3_C_ZERO, queue ));
    }
    else if ( precond->solver == Magma_tally3_ICC ) {
        CHECK( magma_tally3_cvinit( &tmp, Magma_tally3_DEV, A.num_rows, b.num_cols, MAGMA_tally3_C_ZERO, queue ));
    }
    else if ( precond->solver == Magma_tally3_NONE ) {
        magma_tally3_ccopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally3blasSetKernelStream( orig_queue );
        info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_tally3_cmfree( &tmp, queue );
    magma_tally3blasSetKernelStream( orig_queue );
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
    A           magma_tally3_c_matrix
                sparse matrix A

    @param[in]
    b           magma_tally3_c_matrix
                input vector b

    @param[in,out]
    x           magma_tally3_c_matrix*
                output vector x

    @param[in]
    precond     magma_tally3_c_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_applyprecond_left(
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix b,
    magma_tally3_c_matrix *x,
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally3_JACOBI ) {
        CHECK( magma_tally3_cjacobi_diagscal( A.num_rows, precond->d, b, x, queue ));
    }
    else if ( ( precond->solver == Magma_tally3_ILU ||
                precond->solver == Magma_tally3_AILU ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally3_capplycumilu_l( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally3_ILU ||
                precond->solver == Magma_tally3_AILU ) && precond->maxiter < 50 ) {
        magma_tally3_ccopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally3_c_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally3_cjacobiiter_sys( precond->L, b, precond->d, precond->work1, x, &solver_par, queue );
        CHECK( magma_tally3_cjacobispmvupdate(precond->maxiter, precond->L, precond->work1, b, precond->d, x, queue ));
    }
    else if ( ( precond->solver == Magma_tally3_ICC ||
                precond->solver == Magma_tally3_AICC ) && precond->maxiter >= 50 )  {
        CHECK( magma_tally3_capplycumicc_l( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally3_ICC ||
                precond->solver == Magma_tally3_AICC ) && precond->maxiter < 50 )  {
        magma_tally3_ccopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally3_c_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally3_cjacobiiter_sys( precond->L, b, precond->d, precond->work1, x, &solver_par, queue );
        CHECK( magma_tally3_cjacobispmvupdate(precond->maxiter, precond->L, precond->work1, b, precond->d, x, queue ));
    }
    else if ( precond->solver == Magma_tally3_NONE ) {
        magma_tally3_ccopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else if ( precond->solver == Magma_tally3_FUNCTION ) {
        CHECK( magma_tally3_capplycustomprecond_l( b, x, precond, queue ));
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally3blasSetKernelStream( orig_queue );
        info = MAGMA_tally3_ERR_NOT_SUPPORTED; 
    }
cleanup:
    magma_tally3blasSetKernelStream( orig_queue );
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
    A           magma_tally3_c_matrix
                sparse matrix A

    @param[in]
    b           magma_tally3_c_matrix
                input vector b

    @param[in,out]
    x           magma_tally3_c_matrix*
                output vector x

    @param[in]
    precond     magma_tally3_c_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_applyprecond_right(
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix b,
    magma_tally3_c_matrix *x,
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally3_JACOBI ) {
        magma_tally3_ccopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );    // x = b
    }
    else if ( ( precond->solver == Magma_tally3_ILU ||
                precond->solver == Magma_tally3_AILU ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally3_capplycumilu_r( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally3_ILU ||
                precond->solver == Magma_tally3_AILU ) && precond->maxiter < 50 ) {
        magma_tally3_ccopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally3_c_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally3_cjacobiiter_sys( precond->U, b, precond->d2, precond->work2, x, &solver_par, queue );
        CHECK( magma_tally3_cjacobispmvupdate(precond->maxiter, precond->U, precond->work2, b, precond->d2, x, queue ));
    }

    else if ( ( precond->solver == Magma_tally3_ICC ||
                precond->solver == Magma_tally3_AICC ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally3_capplycumicc_r( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally3_ICC ||
               precond->solver == Magma_tally3_AICC ) && precond->maxiter < 50 ) {
        magma_tally3_ccopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally3_c_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally3_cjacobiiter_sys( precond->U, b, precond->d2, precond->work2, x, &solver_par, queue );
        CHECK( magma_tally3_cjacobispmvupdate(precond->maxiter, precond->U, precond->work2, b, precond->d2, x, queue ));
    }
    else if ( precond->solver == Magma_tally3_NONE ) {
        magma_tally3_ccopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else if ( precond->solver == Magma_tally3_FUNCTION ) {
        CHECK( magma_tally3_capplycustomprecond_r( b, x, precond, queue ));
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally3blasSetKernelStream( orig_queue );
        info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_tally3blasSetKernelStream( orig_queue );
    return info;
}


