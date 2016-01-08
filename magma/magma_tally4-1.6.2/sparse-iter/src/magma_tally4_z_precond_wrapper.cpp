/*
    -- MAGMA_tally4 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @precisions normal z -> c d s
       @author Hartwig Anzt

*/
#include "common_magma_tally4sparse.h"




/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and the
    preconditioner parameters, the respective preconditioner
    is chosen. It approximates x for A x = y.

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

    @param[in,out]
    precond     magma_tally4_z_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_z_precond(
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix b,
    magma_tally4_z_matrix *x,
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // set up precond parameters as solver parameters
    magma_tally4_z_solver_par psolver_par;
    psolver_par.epsilon = precond->epsilon;
    psolver_par.maxiter = precond->maxiter;
    psolver_par.restart = precond->restart;
    psolver_par.verbose = 0;
    magma_tally4_z_preconditioner pprecond;
    pprecond.solver = Magma_tally4_NONE;

    switch( precond->solver ) {
        case  Magma_tally4_CG:
                CHECK( magma_tally4_zcg_res( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally4_BICGSTAB:
                CHECK( magma_tally4_zbicgstab( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally4_GMRES:
                CHECK( magma_tally4_zfgmres( A, b, x, &psolver_par, &pprecond, queue )); break;
        case  Magma_tally4_JACOBI:
                CHECK( magma_tally4_zjacobi( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally4_BAITER:
                CHECK( magma_tally4_zbaiter( A, b, x, &psolver_par, queue )); break;
        default:
                CHECK( magma_tally4_zcg_res( A, b, x, &psolver_par, queue )); break;

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
    A           magma_tally4_z_matrix
                sparse matrix A

    @param[in]
    b           magma_tally4_z_matrix
                input vector y

    @param[in,out]
    precond     magma_tally4_z_preconditioner
                preconditioner
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_z_precondsetup(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b,
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_tally4_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_tally4_ERR_NOT_SUPPORTED;
    }

    if ( precond->solver == Magma_tally4_JACOBI ) {
        return magma_tally4_zjacobisetup_diagscal( A, &(precond->d), queue );
    }
    else if ( precond->solver == Magma_tally4_PASTIX ) {
        //return magma_tally4_zpastixsetup( A, b, precond, queue );
        return MAGMA_tally4_ERR_NOT_SUPPORTED;
    }
    else if ( precond->solver == Magma_tally4_ILU ) {
        return magma_tally4_zcumilusetup( A, precond, queue );
    }
    else if ( precond->solver == Magma_tally4_ICC ) {
        return magma_tally4_zcumiccsetup( A, precond, queue );
    }
    else if ( precond->solver == Magma_tally4_AICC ) {
        return magma_tally4_zitericsetup( A, b, precond, queue );
    }
    else if ( precond->solver == Magma_tally4_AILU ) {
        return magma_tally4_ziterilusetup( A, b, precond, queue );
    }
    else if ( precond->solver == Magma_tally4_NONE ) {
        return MAGMA_tally4_SUCCESS;
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        return MAGMA_tally4_ERR_NOT_SUPPORTED;
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
    A           magma_tally4_z_matrix
                sparse matrix A

    @param[in]
    b           magma_tally4_z_matrix
                input vector b

    @param[in,out]
    x           magma_tally4_z_matrix*
                output vector x

    @param[in]
    precond     magma_tally4_z_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_z_applyprecond(
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix b,
    magma_tally4_z_matrix *x,
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_z_matrix tmp={Magma_tally4_CSR};
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally4_JACOBI ) {
        CHECK( magma_tally4_zjacobi_diagscal( A.num_rows, precond->d, b, x, queue ));
    }
    else if ( precond->solver == Magma_tally4_PASTIX ) {
        //CHECK( magma_tally4_zapplypastix( b, x, precond, queue ));
        return MAGMA_tally4_ERR_NOT_SUPPORTED;
    }
    else if ( precond->solver == Magma_tally4_ILU ) {
        CHECK( magma_tally4_zvinit( &tmp, Magma_tally4_DEV, A.num_rows, b.num_cols, MAGMA_tally4_Z_ZERO, queue ));
    }
    else if ( precond->solver == Magma_tally4_ICC ) {
        CHECK( magma_tally4_zvinit( &tmp, Magma_tally4_DEV, A.num_rows, b.num_cols, MAGMA_tally4_Z_ZERO, queue ));
    }
    else if ( precond->solver == Magma_tally4_NONE ) {
        magma_tally4_zcopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally4blasSetKernelStream( orig_queue );
        info = MAGMA_tally4_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_tally4_zmfree( &tmp, queue );
    magma_tally4blasSetKernelStream( orig_queue );
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
    A           magma_tally4_z_matrix
                sparse matrix A

    @param[in]
    b           magma_tally4_z_matrix
                input vector b

    @param[in,out]
    x           magma_tally4_z_matrix*
                output vector x

    @param[in]
    precond     magma_tally4_z_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_z_applyprecond_left(
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix b,
    magma_tally4_z_matrix *x,
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally4_JACOBI ) {
        CHECK( magma_tally4_zjacobi_diagscal( A.num_rows, precond->d, b, x, queue ));
    }
    else if ( ( precond->solver == Magma_tally4_ILU ||
                precond->solver == Magma_tally4_AILU ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally4_zapplycumilu_l( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally4_ILU ||
                precond->solver == Magma_tally4_AILU ) && precond->maxiter < 50 ) {
        magma_tally4_zcopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally4_z_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally4_zjacobiiter_sys( precond->L, b, precond->d, precond->work1, x, &solver_par, queue );
        CHECK( magma_tally4_zjacobispmvupdate(precond->maxiter, precond->L, precond->work1, b, precond->d, x, queue ));
    }
    else if ( ( precond->solver == Magma_tally4_ICC ||
                precond->solver == Magma_tally4_AICC ) && precond->maxiter >= 50 )  {
        CHECK( magma_tally4_zapplycumicc_l( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally4_ICC ||
                precond->solver == Magma_tally4_AICC ) && precond->maxiter < 50 )  {
        magma_tally4_zcopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally4_z_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally4_zjacobiiter_sys( precond->L, b, precond->d, precond->work1, x, &solver_par, queue );
        CHECK( magma_tally4_zjacobispmvupdate(precond->maxiter, precond->L, precond->work1, b, precond->d, x, queue ));
    }
    else if ( precond->solver == Magma_tally4_NONE ) {
        magma_tally4_zcopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else if ( precond->solver == Magma_tally4_FUNCTION ) {
        CHECK( magma_tally4_zapplycustomprecond_l( b, x, precond, queue ));
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally4blasSetKernelStream( orig_queue );
        info = MAGMA_tally4_ERR_NOT_SUPPORTED; 
    }
cleanup:
    magma_tally4blasSetKernelStream( orig_queue );
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
    A           magma_tally4_z_matrix
                sparse matrix A

    @param[in]
    b           magma_tally4_z_matrix
                input vector b

    @param[in,out]
    x           magma_tally4_z_matrix*
                output vector x

    @param[in]
    precond     magma_tally4_z_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_z_applyprecond_right(
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix b,
    magma_tally4_z_matrix *x,
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally4_JACOBI ) {
        magma_tally4_zcopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );    // x = b
    }
    else if ( ( precond->solver == Magma_tally4_ILU ||
                precond->solver == Magma_tally4_AILU ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally4_zapplycumilu_r( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally4_ILU ||
                precond->solver == Magma_tally4_AILU ) && precond->maxiter < 50 ) {
        magma_tally4_zcopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally4_z_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally4_zjacobiiter_sys( precond->U, b, precond->d2, precond->work2, x, &solver_par, queue );
        CHECK( magma_tally4_zjacobispmvupdate(precond->maxiter, precond->U, precond->work2, b, precond->d2, x, queue ));
    }

    else if ( ( precond->solver == Magma_tally4_ICC ||
                precond->solver == Magma_tally4_AICC ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally4_zapplycumicc_r( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally4_ICC ||
               precond->solver == Magma_tally4_AICC ) && precond->maxiter < 50 ) {
        magma_tally4_zcopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally4_z_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally4_zjacobiiter_sys( precond->U, b, precond->d2, precond->work2, x, &solver_par, queue );
        CHECK( magma_tally4_zjacobispmvupdate(precond->maxiter, precond->U, precond->work2, b, precond->d2, x, queue ));
    }
    else if ( precond->solver == Magma_tally4_NONE ) {
        magma_tally4_zcopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else if ( precond->solver == Magma_tally4_FUNCTION ) {
        CHECK( magma_tally4_zapplycustomprecond_r( b, x, precond, queue ));
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally4blasSetKernelStream( orig_queue );
        info = MAGMA_tally4_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_tally4blasSetKernelStream( orig_queue );
    return info;
}


