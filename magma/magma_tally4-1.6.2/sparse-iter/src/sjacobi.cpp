/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zjacobi.cpp normal z -> s, Sun May  3 11:22:59 2015
*/

#include "common_magma_tally4sparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a real symmetric N-by-N positive definite matrix A.
    This is a GPU implementation of the Jacobi method.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_s_matrix
                input matrix A

    @param[in]
    b           magma_tally4_s_matrix
                RHS b

    @param[in,out]
    x           magma_tally4_s_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_tally4_s_solver_par*
                solver parameters

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_sgesv
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sjacobi(
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix b,
    magma_tally4_s_matrix *x,
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    

    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );
    
    // some useful variables
    float c_zero = MAGMA_tally4_S_ZERO;
    
    float nom0 = 0.0;

    magma_tally4_s_matrix r={Magma_tally4_CSR}, d={Magma_tally4_CSR}, ACSR={Magma_tally4_CSR} ;
    
    CHECK( magma_tally4_smconvert(A, &ACSR, A.storage_type, Magma_tally4_CSR, queue ) );

    // prepare solver feedback
    solver_par->solver = Magma_tally4_JACOBI;
    solver_par->info = MAGMA_tally4_SUCCESS;

    real_Double_t tempo1, tempo2;
    float residual;
    // solver setup
    CHECK( magma_tally4_svinit( &r, Magma_tally4_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK(  magma_tally4_sresidualvec( ACSR, b, *x, &r, &residual, queue));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;
    nom0 = residual;

    // Jacobi setup
    CHECK( magma_tally4_sjacobisetup_diagscal( ACSR, &d, queue ));
    magma_tally4_s_solver_par jacobiiter_par;
    jacobiiter_par.maxiter = solver_par->maxiter;

    tempo1 = magma_tally4_sync_wtime( queue );

    // Jacobi iterator
    CHECK( magma_tally4_sjacobispmvupdate(jacobiiter_par.maxiter, ACSR, r, b, d, x, queue ));

    tempo2 = magma_tally4_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_tally4_sresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res )
        info = MAGMA_tally4_SUCCESS;
    else
        info = MAGMA_tally4_DIVERGENCE;
    
cleanup:
    magma_tally4_smfree( &r, queue );
    magma_tally4_smfree( &d, queue );
    magma_tally4_smfree( &ACSR, queue );

    magma_tally4blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally4_sjacobi */






/**
    Purpose
    -------

    Prepares the Matrix M for the Jacobi Iteration according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    It returns the preconditioner Matrix M and a vector d
    containing the diagonal elements.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_s_matrix
                input matrix A

    @param[in]
    M           magma_tally4_s_matrix*
                M = D^(-1) * (L+U)

    @param[in,out]
    d           magma_tally4_s_matrix*
                vector with diagonal elements of A
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_s
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sjacobisetup_matrix(
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix *M, magma_tally4_s_matrix *d,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_int_t i;

    magma_tally4_s_matrix A_h1={Magma_tally4_CSR}, A_h2={Magma_tally4_CSR}, B={Magma_tally4_CSR}, C={Magma_tally4_CSR};
    magma_tally4_s_matrix diag={Magma_tally4_CSR};
    CHECK( magma_tally4_svinit( &diag, Magma_tally4_CPU, A.num_rows, 1, MAGMA_tally4_S_ZERO, queue ));

    if ( A.storage_type != Magma_tally4_CSR) {
        CHECK( magma_tally4_smtransfer( A, &A_h1, A.memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_smconvert( A_h1, &B, A_h1.storage_type, Magma_tally4_CSR, queue ));
    }
    else {
        CHECK( magma_tally4_smtransfer( A, &B, A.memory_location, Magma_tally4_CPU, queue ));
    }
    for( magma_tally4_int_t rowindex=0; rowindex<B.num_rows; rowindex++ ) {
        magma_tally4_int_t start = (B.drow[rowindex]);
        magma_tally4_int_t end = (B.drow[rowindex+1]);
        for( i=start; i<end; i++ ) {
            if ( B.dcol[i]==rowindex ) {
                diag.val[rowindex] = B.val[i];
                if ( MAGMA_tally4_S_REAL( diag.val[rowindex]) == 0 )
                    printf(" error: zero diagonal element in row %d!\n",
                                                                (int) rowindex);
            }
        }
        for( i=start; i<end; i++ ) {
            B.val[i] = B.val[i] / diag.val[rowindex];
            if ( B.dcol[i]==rowindex ) {
                B.val[i] = MAGMA_tally4_S_MAKE( 0., 0. );
            }
        }
    }
    CHECK( magma_tally4_s_csr_compressor(&B.val, &B.drow, &B.dcol,
                           &C.val, &C.drow, &C.dcol, &B.num_rows, queue ));
    C.num_rows = B.num_rows;
    C.num_cols = B.num_cols;
    C.memory_location = B.memory_location;
    C.nnz = C.drow[B.num_rows];
    C.storage_type = B.storage_type;
    C.memory_location = B.memory_location;
    if ( A.storage_type != Magma_tally4_CSR) {
        CHECK( magma_tally4_smconvert( C, &A_h2, Magma_tally4_CSR, A_h1.storage_type, queue ));
        CHECK( magma_tally4_smtransfer( A_h2, M, Magma_tally4_CPU, A.memory_location, queue ));
    }
    else {
        CHECK( magma_tally4_smtransfer( C, M, Magma_tally4_CPU, A.memory_location, queue ));
    }
    CHECK( magma_tally4_smtransfer( diag, d, Magma_tally4_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_tally4_CSR) {
        magma_tally4_smfree( &A_h1, queue );
        magma_tally4_smfree( &A_h2, queue );
    }
    
cleanup:
    magma_tally4_smfree( &B, queue );
    magma_tally4_smfree( &C, queue );

    magma_tally4_smfree( &diag, queue );
 
    return info;
}


/**
    Purpose
    -------


    It returns a vector d
    containing the inverse diagonal elements.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_s_matrix
                input matrix A

    @param[in,out]
    d           magma_tally4_s_matrix*
                vector with diagonal elements
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_s
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sjacobisetup_diagscal(
    magma_tally4_s_matrix A, magma_tally4_s_matrix *d,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_int_t i;

    magma_tally4_s_matrix A_h1={Magma_tally4_CSR}, B={Magma_tally4_CSR};
    magma_tally4_s_matrix diag={Magma_tally4_CSR};
    CHECK( magma_tally4_svinit( &diag, Magma_tally4_CPU, A.num_rows, 1, MAGMA_tally4_S_ZERO, queue ));

    if ( A.storage_type != Magma_tally4_CSR) {
        CHECK( magma_tally4_smtransfer( A, &A_h1, A.memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_smconvert( A_h1, &B, A_h1.storage_type, Magma_tally4_CSR, queue ));
    }
    else {
        CHECK( magma_tally4_smtransfer( A, &B, A.memory_location, Magma_tally4_CPU, queue ));
    }
    for( magma_tally4_int_t rowindex=0; rowindex<B.num_rows; rowindex++ ) {
        magma_tally4_int_t start = (B.drow[rowindex]);
        magma_tally4_int_t end = (B.drow[rowindex+1]);
        for( i=start; i<end; i++ ) {
            if ( B.dcol[i]==rowindex ) {
                diag.val[rowindex] = 1.0/B.val[i];
                break;
            }
        }
        if ( diag.val[rowindex] == MAGMA_tally4_S_ZERO ){
            printf(" error: zero diagonal element in row %d!\n",
                                                        (int) rowindex);
            
            if ( A.storage_type != Magma_tally4_CSR) {
                magma_tally4_smfree( &A_h1, queue );
            }
            magma_tally4_smfree( &B, queue );
            magma_tally4_smfree( &diag, queue );
            info = MAGMA_tally4_ERR_BADPRECOND;
            goto cleanup;
        }
    }
    CHECK( magma_tally4_smtransfer( diag, d, Magma_tally4_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_tally4_CSR) {
        magma_tally4_smfree( &A_h1, queue );
    }
    
cleanup:
    magma_tally4_smfree( &B, queue );
    magma_tally4_smfree( &diag, queue );
 
    return info;
}



/**
    Purpose
    -------

    Prepares the Jacobi Iteration according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    Returns the vector c

    Arguments
    ---------

    @param[in]
    b           magma_tally4_s_matrix
                RHS b

    @param[in]
    d           magma_tally4_s_matrix
                vector with diagonal entries

    @param[in]
    c           magma_tally4_s_matrix*
                c = D^(-1) * b
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_s
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sjacobisetup_vector(
    magma_tally4_s_matrix b, magma_tally4_s_matrix d,
    magma_tally4_s_matrix *c,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_s_matrix diag={Magma_tally4_CSR}, c_t={Magma_tally4_CSR}, b_h={Magma_tally4_CSR}, tmp={Magma_tally4_CSR};
    
    if ( b.memory_location == Magma_tally4_CPU ) {

        CHECK( magma_tally4_svinit( &c_t, Magma_tally4_CPU, b.num_rows, b.num_cols, MAGMA_tally4_S_ZERO, queue ));

        CHECK( magma_tally4_smtransfer( b, &b_h, b.memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_smtransfer( d, &diag, b.memory_location, Magma_tally4_CPU, queue ));

        for( magma_tally4_int_t rowindex=0; rowindex<b.num_rows; rowindex++ ) {
            c_t.val[rowindex] = b_h.val[rowindex] / diag.val[rowindex];

        }
        CHECK( magma_tally4_smtransfer( c_t, c, Magma_tally4_CPU, b.memory_location, queue ));
    }
    else if ( b.memory_location == Magma_tally4_DEV ) {
        // fill vector
        CHECK( magma_tally4_svinit( &tmp, Magma_tally4_DEV, b.num_rows, b.num_cols, MAGMA_tally4_S_ZERO, queue ));
        CHECK( magma_tally4_sjacobisetup_vector_gpu(
                    b.num_rows, b, d, *c, &tmp, queue ));
        goto cleanup;
    }

cleanup:
    magma_tally4_smfree( &tmp, queue );
    magma_tally4_smfree( &diag, queue );
    magma_tally4_smfree( &c_t, queue );
    magma_tally4_smfree( &b_h, queue );
    
    return info;
}


/**
    Purpose
    -------

    Prepares the Jacobi Iteration according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_s_matrix
                input matrix A

    @param[in]
    b           magma_tally4_s_matrix
                RHS b

    @param[in]
    M           magma_tally4_s_matrix*
                M = D^(-1) * (L+U)

    @param[in]
    c           magma_tally4_s_matrix*
                c = D^(-1) * b
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_s
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sjacobisetup(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b,
    magma_tally4_s_matrix *M, magma_tally4_s_matrix *c,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_int_t i;

    magma_tally4_s_matrix A_h1={Magma_tally4_CSR}, A_h2={Magma_tally4_CSR}, B={Magma_tally4_CSR}, C={Magma_tally4_CSR};
    magma_tally4_s_matrix diag={Magma_tally4_CSR}, c_t={Magma_tally4_CSR}, b_h={Magma_tally4_CSR};
    CHECK( magma_tally4_svinit( &c_t, Magma_tally4_CPU, A.num_rows, b.num_cols, MAGMA_tally4_S_ZERO, queue ));
    CHECK( magma_tally4_svinit( &diag, Magma_tally4_CPU, A.num_rows, b.num_cols, MAGMA_tally4_S_ZERO, queue ));
    CHECK( magma_tally4_smtransfer( b, &b_h, A.memory_location, Magma_tally4_CPU, queue ));

    if ( A.storage_type != Magma_tally4_CSR ) {
        CHECK( magma_tally4_smtransfer( A, &A_h1, A.memory_location, Magma_tally4_CPU, queue ));
        CHECK( magma_tally4_smconvert( A_h1, &B, A_h1.storage_type, Magma_tally4_CSR, queue ));
    }
    else {
        CHECK( magma_tally4_smtransfer( A, &B, A.memory_location, Magma_tally4_CPU, queue ));
    }
    for( magma_tally4_int_t rowindex=0; rowindex<B.num_rows; rowindex++ ) {
        magma_tally4_int_t start = (B.drow[rowindex]);
        magma_tally4_int_t end = (B.drow[rowindex+1]);
        for( i=start; i<end; i++ ) {
            if ( B.dcol[i]==rowindex ) {
                diag.val[rowindex] = B.val[i];
                if ( MAGMA_tally4_S_REAL( diag.val[rowindex]) == 0 )
                    printf(" error: zero diagonal element in row %d!\n",
                                                               (int) rowindex);
            }
        }
        for( i=start; i<end; i++ ) {
            B.val[i] = B.val[i] / diag.val[rowindex];
            if ( B.dcol[i]==rowindex ) {
                B.val[i] = MAGMA_tally4_S_MAKE( 0., 0. );
            }
        }
        c_t.val[rowindex] = b_h.val[rowindex] / diag.val[rowindex];

    }

    CHECK( magma_tally4_s_csr_compressor(&B.val, &B.drow, &B.dcol,
                           &C.val, &C.drow, &C.dcol, &B.num_rows, queue ));

    C.num_rows = B.num_rows;
    C.num_cols = B.num_cols;
    C.memory_location = B.memory_location;
    C.nnz = C.drow[B.num_rows];
    C.storage_type = B.storage_type;
    C.memory_location = B.memory_location;
    if ( A.storage_type != Magma_tally4_CSR) {
        A_h2.alignment = A.alignment;
        A_h2.blocksize = A.blocksize;
        CHECK( magma_tally4_smconvert( C, &A_h2, Magma_tally4_CSR, A_h1.storage_type, queue ));
        CHECK( magma_tally4_smtransfer( A_h2, M, Magma_tally4_CPU, A.memory_location, queue ));
    }
    else {
        CHECK( magma_tally4_smtransfer( C, M, Magma_tally4_CPU, A.memory_location, queue ));
    }
    CHECK( magma_tally4_smtransfer( c_t, c, Magma_tally4_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_tally4_CSR) {
        magma_tally4_smfree( &A_h1, queue );
        magma_tally4_smfree( &A_h2, queue );
    }
    
cleanup:
    magma_tally4_smfree( &B, queue );
    magma_tally4_smfree( &C, queue );
    magma_tally4_smfree( &diag, queue );
    magma_tally4_smfree( &c_t, queue );
    magma_tally4_smfree( &b_h, queue );

    return info;
}



/**
    Purpose
    -------

    Iterates the solution approximation according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    This routine takes the iteration matrix M as input.

    Arguments
    ---------

    @param[in]
    M           magma_tally4_s_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_tally4_s_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_tally4_s_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_tally4_s_solver_par*
                solver parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_s
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sjacobiiter(
    magma_tally4_s_matrix M, magma_tally4_s_matrix c, magma_tally4_s_matrix *x,
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    // local variables
    float c_zero = MAGMA_tally4_S_ZERO, c_one = MAGMA_tally4_S_ONE,
                                            c_mone = MAGMA_tally4_S_NEG_ONE;
    magma_tally4_int_t dofs = M.num_rows*x->num_cols;
    magma_tally4_s_matrix t={Magma_tally4_CSR}, swap={Magma_tally4_CSR};
    CHECK( magma_tally4_svinit( &t, Magma_tally4_DEV, M.num_rows, x->num_cols, c_zero, queue ));


    for( magma_tally4_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_tally4_s_spmv( c_mone, M, *x, c_zero, t, queue ));        // t = - M * x
        magma_tally4_saxpy( dofs, c_one , c.dval, 1 , t.dval, 1 );        // t = t + c

        // swap so that x again contains solution, and y is ready to be used
        swap = *x;
        *x = t;
        t = swap;
        //magma_tally4_scopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }

cleanup:
    magma_tally4_smfree( &t, queue );

    magma_tally4blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally4_sjacobiiter */



/**
    Purpose
    -------

    Iterates the solution approximation according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    Arguments
    ---------

    @param[in]
    M           magma_tally4_s_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_tally4_s_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_tally4_s_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_tally4_s_solver_par*
                solver parameters

    @param[in]
    solver_par  magma_tally4_s_precond_par*
                precond parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_s
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sjacobiiter_precond(
    magma_tally4_s_matrix M, magma_tally4_s_matrix *x,
    magma_tally4_s_solver_par *solver_par, magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    // local variables
    float c_zero = MAGMA_tally4_S_ZERO, c_one = MAGMA_tally4_S_ONE,
                                            c_mone = MAGMA_tally4_S_NEG_ONE;
    magma_tally4_int_t dofs = M.num_rows;
    magma_tally4_int_t num_vecs = x->num_rows / dofs;
    magma_tally4_s_matrix swap={Magma_tally4_CSR};

    for( magma_tally4_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_tally4_s_spmv( c_mone, M, *x, c_zero, precond->work2, queue )); // t = - M * x

        magma_tally4_saxpy( num_vecs*dofs, c_one ,
                precond->work1.dval, 1 , precond->work2.dval, 1 ); // t = t + c

        // swap so that x again contains solution, and y is ready to be used
        swap = *x;
        *x = precond->work2;
        precond->work2 = swap;
        //magma_tally4_scopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }
    
cleanup:
    magma_tally4blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally4_sjacobiiter */



    /**
    Purpose
    -------

    Iterates the solution approximation according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    This routine takes the system matrix A and the RHS b as input.

    Arguments
    ---------

    @param[in]
    M           magma_tally4_s_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_tally4_s_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_tally4_s_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_tally4_s_solver_par*
                solver parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_s
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sjacobiiter_sys(
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix b,
    magma_tally4_s_matrix d,
    magma_tally4_s_matrix t,
    magma_tally4_s_matrix *x,
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    // local variables
    float c_zero = MAGMA_tally4_S_ZERO, c_one = MAGMA_tally4_S_ONE;

    for( magma_tally4_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_tally4_s_spmv( c_one, A, *x, c_zero, t, queue ));        // t =  A * x
        CHECK( magma_tally4_sjacobiupdate( t, b, d, x, queue ));
        // swap so that x again contains solution, and y is ready to be used
        //swap = *x;
        //*x = t;
        //t = swap;
        //magma_tally4_scopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }
    
cleanup:
    magma_tally4blasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_tally4_sjacobiiter_sys */
