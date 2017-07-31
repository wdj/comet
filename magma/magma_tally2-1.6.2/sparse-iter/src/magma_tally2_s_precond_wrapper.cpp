/*
    -- MAGMA_tally2 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @generated from magma_tally2_z_precond_wrapper.cpp normal z -> s, Mon May  4 11:57:23 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"




/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and the
    preconditioner parameters, the respective preconditioner
    is chosen. It approximates x for A x = y.

    Arguments
    ---------

    @param[in]
    A           magma_tally2_s_matrix
                sparse matrix A

    @param[in]
    b           magma_tally2_s_matrix
                input vector b

    @param[in]
    x           magma_tally2_s_matrix*
                output vector x

    @param[in,out]
    precond     magma_tally2_s_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_s_precond(
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    // set up precond parameters as solver parameters
    magma_tally2_s_solver_par psolver_par;
    psolver_par.epsilon = precond->epsilon;
    psolver_par.maxiter = precond->maxiter;
    psolver_par.restart = precond->restart;
    psolver_par.verbose = 0;
    magma_tally2_s_preconditioner pprecond;
    pprecond.solver = Magma_tally2_NONE;

    switch( precond->solver ) {
        case  Magma_tally2_CG:
                CHECK( magma_tally2_scg_res( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally2_BICGSTAB:
                CHECK( magma_tally2_sbicgstab( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally2_GMRES:
                CHECK( magma_tally2_sfgmres( A, b, x, &psolver_par, &pprecond, queue )); break;
        case  Magma_tally2_JACOBI:
                CHECK( magma_tally2_sjacobi( A, b, x, &psolver_par, queue )); break;
        case  Magma_tally2_BAITER:
                CHECK( magma_tally2_sbaiter( A, b, x, &psolver_par, queue )); break;
        default:
                CHECK( magma_tally2_scg_res( A, b, x, &psolver_par, queue )); break;

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
    A           magma_tally2_s_matrix
                sparse matrix A

    @param[in]
    b           magma_tally2_s_matrix
                input vector y

    @param[in,out]
    precond     magma_tally2_s_preconditioner
                preconditioner
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_s_precondsetup(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    
    // make sure RHS is a dense matrix
    if ( b.storage_type != Magma_tally2_DENSE ) {
        printf( "error: sparse RHS not yet supported.\n" );
        return MAGMA_tally2_ERR_NOT_SUPPORTED;
    }

    if ( precond->solver == Magma_tally2_JACOBI ) {
        return magma_tally2_sjacobisetup_diagscal( A, &(precond->d), queue );
    }
    else if ( precond->solver == Magma_tally2_PASTIX ) {
        //return magma_tally2_spastixsetup( A, b, precond, queue );
        return MAGMA_tally2_ERR_NOT_SUPPORTED;
    }
    else if ( precond->solver == Magma_tally2_ILU ) {
        return magma_tally2_scumilusetup( A, precond, queue );
    }
    else if ( precond->solver == Magma_tally2_ICC ) {
        return magma_tally2_scumiccsetup( A, precond, queue );
    }
    else if ( precond->solver == Magma_tally2_AICC ) {
        return magma_tally2_sitericsetup( A, b, precond, queue );
    }
    else if ( precond->solver == Magma_tally2_AILU ) {
        return magma_tally2_siterilusetup( A, b, precond, queue );
    }
    else if ( precond->solver == Magma_tally2_NONE ) {
        return MAGMA_tally2_SUCCESS;
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        return MAGMA_tally2_ERR_NOT_SUPPORTED;
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
    A           magma_tally2_s_matrix
                sparse matrix A

    @param[in]
    b           magma_tally2_s_matrix
                input vector b

    @param[in,out]
    x           magma_tally2_s_matrix*
                output vector x

    @param[in]
    precond     magma_tally2_s_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_s_applyprecond(
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_s_matrix tmp={Magma_tally2_CSR};
    
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue=NULL;
    magma_tally2blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally2_JACOBI ) {
        CHECK( magma_tally2_sjacobi_diagscal( A.num_rows, precond->d, b, x, queue ));
    }
    else if ( precond->solver == Magma_tally2_PASTIX ) {
        //CHECK( magma_tally2_sapplypastix( b, x, precond, queue ));
        return MAGMA_tally2_ERR_NOT_SUPPORTED;
    }
    else if ( precond->solver == Magma_tally2_ILU ) {
        CHECK( magma_tally2_svinit( &tmp, Magma_tally2_DEV, A.num_rows, b.num_cols, MAGMA_tally2_S_ZERO, queue ));
    }
    else if ( precond->solver == Magma_tally2_ICC ) {
        CHECK( magma_tally2_svinit( &tmp, Magma_tally2_DEV, A.num_rows, b.num_cols, MAGMA_tally2_S_ZERO, queue ));
    }
    else if ( precond->solver == Magma_tally2_NONE ) {
        magma_tally2_scopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally2blasSetKernelStream( orig_queue );
        info = MAGMA_tally2_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_tally2_smfree( &tmp, queue );
    magma_tally2blasSetKernelStream( orig_queue );
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
    A           magma_tally2_s_matrix
                sparse matrix A

    @param[in]
    b           magma_tally2_s_matrix
                input vector b

    @param[in,out]
    x           magma_tally2_s_matrix*
                output vector x

    @param[in]
    precond     magma_tally2_s_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_s_applyprecond_left(
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue=NULL;
    magma_tally2blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally2_JACOBI ) {
        CHECK( magma_tally2_sjacobi_diagscal( A.num_rows, precond->d, b, x, queue ));
    }
    else if ( ( precond->solver == Magma_tally2_ILU ||
                precond->solver == Magma_tally2_AILU ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally2_sapplycumilu_l( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally2_ILU ||
                precond->solver == Magma_tally2_AILU ) && precond->maxiter < 50 ) {
        magma_tally2_scopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally2_s_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally2_sjacobiiter_sys( precond->L, b, precond->d, precond->work1, x, &solver_par, queue );
        CHECK( magma_tally2_sjacobispmvupdate(precond->maxiter, precond->L, precond->work1, b, precond->d, x, queue ));
    }
    else if ( ( precond->solver == Magma_tally2_ICC ||
                precond->solver == Magma_tally2_AICC ) && precond->maxiter >= 50 )  {
        CHECK( magma_tally2_sapplycumicc_l( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally2_ICC ||
                precond->solver == Magma_tally2_AICC ) && precond->maxiter < 50 )  {
        magma_tally2_scopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally2_s_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally2_sjacobiiter_sys( precond->L, b, precond->d, precond->work1, x, &solver_par, queue );
        CHECK( magma_tally2_sjacobispmvupdate(precond->maxiter, precond->L, precond->work1, b, precond->d, x, queue ));
    }
    else if ( precond->solver == Magma_tally2_NONE ) {
        magma_tally2_scopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else if ( precond->solver == Magma_tally2_FUNCTION ) {
        CHECK( magma_tally2_sapplycustomprecond_l( b, x, precond, queue ));
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally2blasSetKernelStream( orig_queue );
        info = MAGMA_tally2_ERR_NOT_SUPPORTED; 
    }
cleanup:
    magma_tally2blasSetKernelStream( orig_queue );
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
    A           magma_tally2_s_matrix
                sparse matrix A

    @param[in]
    b           magma_tally2_s_matrix
                input vector b

    @param[in,out]
    x           magma_tally2_s_matrix*
                output vector x

    @param[in]
    precond     magma_tally2_s_preconditioner
                preconditioner

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_saux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_s_applyprecond_right(
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue=NULL;
    magma_tally2blasGetKernelStream( &orig_queue );

    if ( precond->solver == Magma_tally2_JACOBI ) {
        magma_tally2_scopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );    // x = b
    }
    else if ( ( precond->solver == Magma_tally2_ILU ||
                precond->solver == Magma_tally2_AILU ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally2_sapplycumilu_r( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally2_ILU ||
                precond->solver == Magma_tally2_AILU ) && precond->maxiter < 50 ) {
        magma_tally2_scopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally2_s_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally2_sjacobiiter_sys( precond->U, b, precond->d2, precond->work2, x, &solver_par, queue );
        CHECK( magma_tally2_sjacobispmvupdate(precond->maxiter, precond->U, precond->work2, b, precond->d2, x, queue ));
    }

    else if ( ( precond->solver == Magma_tally2_ICC ||
                precond->solver == Magma_tally2_AICC ) && precond->maxiter >= 50 ) {
        CHECK( magma_tally2_sapplycumicc_r( b, x, precond, queue ));
    }
    else if ( ( precond->solver == Magma_tally2_ICC ||
               precond->solver == Magma_tally2_AICC ) && precond->maxiter < 50 ) {
        magma_tally2_scopy( b.num_rows*b.num_cols, b.dval, b.num_cols, x->dval, b.num_cols );
        magma_tally2_s_solver_par solver_par;
        solver_par.maxiter = precond->maxiter;
        //magma_tally2_sjacobiiter_sys( precond->U, b, precond->d2, precond->work2, x, &solver_par, queue );
        CHECK( magma_tally2_sjacobispmvupdate(precond->maxiter, precond->U, precond->work2, b, precond->d2, x, queue ));
    }
    else if ( precond->solver == Magma_tally2_NONE ) {
        magma_tally2_scopy( b.num_rows*b.num_cols, b.dval, 1, x->dval, 1 );      //  x = b
    }
    else if ( precond->solver == Magma_tally2_FUNCTION ) {
        CHECK( magma_tally2_sapplycustomprecond_r( b, x, precond, queue ));
    }
    else {
        printf( "error: preconditioner type not yet supported.\n" );
        magma_tally2blasSetKernelStream( orig_queue );
        info = MAGMA_tally2_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_tally2blasSetKernelStream( orig_queue );
    return info;
}


