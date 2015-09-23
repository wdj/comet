/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zjacobi.cpp normal z -> c, Sun May  3 11:22:59 2015
*/

#include "common_magma_minproductsparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a complex Hermitian N-by-N positive definite matrix A.
    This is a GPU implementation of the Jacobi method.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_c_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_c_matrix
                RHS b

    @param[in,out]
    x           magma_minproduct_c_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_minproduct_c_solver_par*
                solver parameters

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_cgesv
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cjacobi(
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix b,
    magma_minproduct_c_matrix *x,
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    

    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );
    
    // some useful variables
    magma_minproductFloatComplex c_zero = MAGMA_minproduct_C_ZERO;
    
    float nom0 = 0.0;

    magma_minproduct_c_matrix r={Magma_minproduct_CSR}, d={Magma_minproduct_CSR}, ACSR={Magma_minproduct_CSR} ;
    
    CHECK( magma_minproduct_cmconvert(A, &ACSR, A.storage_type, Magma_minproduct_CSR, queue ) );

    // prepare solver feedback
    solver_par->solver = Magma_minproduct_JACOBI;
    solver_par->info = MAGMA_minproduct_SUCCESS;

    real_Double_t tempo1, tempo2;
    float residual;
    // solver setup
    CHECK( magma_minproduct_cvinit( &r, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK(  magma_minproduct_cresidualvec( ACSR, b, *x, &r, &residual, queue));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;
    nom0 = residual;

    // Jacobi setup
    CHECK( magma_minproduct_cjacobisetup_diagscal( ACSR, &d, queue ));
    magma_minproduct_c_solver_par jacobiiter_par;
    jacobiiter_par.maxiter = solver_par->maxiter;

    tempo1 = magma_minproduct_sync_wtime( queue );

    // Jacobi iterator
    CHECK( magma_minproduct_cjacobispmvupdate(jacobiiter_par.maxiter, ACSR, r, b, d, x, queue ));

    tempo2 = magma_minproduct_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_minproduct_cresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res )
        info = MAGMA_minproduct_SUCCESS;
    else
        info = MAGMA_minproduct_DIVERGENCE;
    
cleanup:
    magma_minproduct_cmfree( &r, queue );
    magma_minproduct_cmfree( &d, queue );
    magma_minproduct_cmfree( &ACSR, queue );

    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_cjacobi */






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
    A           magma_minproduct_c_matrix
                input matrix A

    @param[in]
    M           magma_minproduct_c_matrix*
                M = D^(-1) * (L+U)

    @param[in,out]
    d           magma_minproduct_c_matrix*
                vector with diagonal elements of A
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cjacobisetup_matrix(
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix *M, magma_minproduct_c_matrix *d,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_int_t i;

    magma_minproduct_c_matrix A_h1={Magma_minproduct_CSR}, A_h2={Magma_minproduct_CSR}, B={Magma_minproduct_CSR}, C={Magma_minproduct_CSR};
    magma_minproduct_c_matrix diag={Magma_minproduct_CSR};
    CHECK( magma_minproduct_cvinit( &diag, Magma_minproduct_CPU, A.num_rows, 1, MAGMA_minproduct_C_ZERO, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        CHECK( magma_minproduct_cmtransfer( A, &A_h1, A.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_cmconvert( A_h1, &B, A_h1.storage_type, Magma_minproduct_CSR, queue ));
    }
    else {
        CHECK( magma_minproduct_cmtransfer( A, &B, A.memory_location, Magma_minproduct_CPU, queue ));
    }
    for( magma_minproduct_int_t rowindex=0; rowindex<B.num_rows; rowindex++ ) {
        magma_minproduct_int_t start = (B.drow[rowindex]);
        magma_minproduct_int_t end = (B.drow[rowindex+1]);
        for( i=start; i<end; i++ ) {
            if ( B.dcol[i]==rowindex ) {
                diag.val[rowindex] = B.val[i];
                if ( MAGMA_minproduct_C_REAL( diag.val[rowindex]) == 0 )
                    printf(" error: zero diagonal element in row %d!\n",
                                                                (int) rowindex);
            }
        }
        for( i=start; i<end; i++ ) {
            B.val[i] = B.val[i] / diag.val[rowindex];
            if ( B.dcol[i]==rowindex ) {
                B.val[i] = MAGMA_minproduct_C_MAKE( 0., 0. );
            }
        }
    }
    CHECK( magma_minproduct_c_csr_compressor(&B.val, &B.drow, &B.dcol,
                           &C.val, &C.drow, &C.dcol, &B.num_rows, queue ));
    C.num_rows = B.num_rows;
    C.num_cols = B.num_cols;
    C.memory_location = B.memory_location;
    C.nnz = C.drow[B.num_rows];
    C.storage_type = B.storage_type;
    C.memory_location = B.memory_location;
    if ( A.storage_type != Magma_minproduct_CSR) {
        CHECK( magma_minproduct_cmconvert( C, &A_h2, Magma_minproduct_CSR, A_h1.storage_type, queue ));
        CHECK( magma_minproduct_cmtransfer( A_h2, M, Magma_minproduct_CPU, A.memory_location, queue ));
    }
    else {
        CHECK( magma_minproduct_cmtransfer( C, M, Magma_minproduct_CPU, A.memory_location, queue ));
    }
    CHECK( magma_minproduct_cmtransfer( diag, d, Magma_minproduct_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        magma_minproduct_cmfree( &A_h1, queue );
        magma_minproduct_cmfree( &A_h2, queue );
    }
    
cleanup:
    magma_minproduct_cmfree( &B, queue );
    magma_minproduct_cmfree( &C, queue );

    magma_minproduct_cmfree( &diag, queue );
 
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
    A           magma_minproduct_c_matrix
                input matrix A

    @param[in,out]
    d           magma_minproduct_c_matrix*
                vector with diagonal elements
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cjacobisetup_diagscal(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix *d,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_int_t i;

    magma_minproduct_c_matrix A_h1={Magma_minproduct_CSR}, B={Magma_minproduct_CSR};
    magma_minproduct_c_matrix diag={Magma_minproduct_CSR};
    CHECK( magma_minproduct_cvinit( &diag, Magma_minproduct_CPU, A.num_rows, 1, MAGMA_minproduct_C_ZERO, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        CHECK( magma_minproduct_cmtransfer( A, &A_h1, A.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_cmconvert( A_h1, &B, A_h1.storage_type, Magma_minproduct_CSR, queue ));
    }
    else {
        CHECK( magma_minproduct_cmtransfer( A, &B, A.memory_location, Magma_minproduct_CPU, queue ));
    }
    for( magma_minproduct_int_t rowindex=0; rowindex<B.num_rows; rowindex++ ) {
        magma_minproduct_int_t start = (B.drow[rowindex]);
        magma_minproduct_int_t end = (B.drow[rowindex+1]);
        for( i=start; i<end; i++ ) {
            if ( B.dcol[i]==rowindex ) {
                diag.val[rowindex] = 1.0/B.val[i];
                break;
            }
        }
        if ( diag.val[rowindex] == MAGMA_minproduct_C_ZERO ){
            printf(" error: zero diagonal element in row %d!\n",
                                                        (int) rowindex);
            
            if ( A.storage_type != Magma_minproduct_CSR) {
                magma_minproduct_cmfree( &A_h1, queue );
            }
            magma_minproduct_cmfree( &B, queue );
            magma_minproduct_cmfree( &diag, queue );
            info = MAGMA_minproduct_ERR_BADPRECOND;
            goto cleanup;
        }
    }
    CHECK( magma_minproduct_cmtransfer( diag, d, Magma_minproduct_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        magma_minproduct_cmfree( &A_h1, queue );
    }
    
cleanup:
    magma_minproduct_cmfree( &B, queue );
    magma_minproduct_cmfree( &diag, queue );
 
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
    b           magma_minproduct_c_matrix
                RHS b

    @param[in]
    d           magma_minproduct_c_matrix
                vector with diagonal entries

    @param[in]
    c           magma_minproduct_c_matrix*
                c = D^(-1) * b
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cjacobisetup_vector(
    magma_minproduct_c_matrix b, magma_minproduct_c_matrix d,
    magma_minproduct_c_matrix *c,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_c_matrix diag={Magma_minproduct_CSR}, c_t={Magma_minproduct_CSR}, b_h={Magma_minproduct_CSR}, tmp={Magma_minproduct_CSR};
    
    if ( b.memory_location == Magma_minproduct_CPU ) {

        CHECK( magma_minproduct_cvinit( &c_t, Magma_minproduct_CPU, b.num_rows, b.num_cols, MAGMA_minproduct_C_ZERO, queue ));

        CHECK( magma_minproduct_cmtransfer( b, &b_h, b.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_cmtransfer( d, &diag, b.memory_location, Magma_minproduct_CPU, queue ));

        for( magma_minproduct_int_t rowindex=0; rowindex<b.num_rows; rowindex++ ) {
            c_t.val[rowindex] = b_h.val[rowindex] / diag.val[rowindex];

        }
        CHECK( magma_minproduct_cmtransfer( c_t, c, Magma_minproduct_CPU, b.memory_location, queue ));
    }
    else if ( b.memory_location == Magma_minproduct_DEV ) {
        // fill vector
        CHECK( magma_minproduct_cvinit( &tmp, Magma_minproduct_DEV, b.num_rows, b.num_cols, MAGMA_minproduct_C_ZERO, queue ));
        CHECK( magma_minproduct_cjacobisetup_vector_gpu(
                    b.num_rows, b, d, *c, &tmp, queue ));
        goto cleanup;
    }

cleanup:
    magma_minproduct_cmfree( &tmp, queue );
    magma_minproduct_cmfree( &diag, queue );
    magma_minproduct_cmfree( &c_t, queue );
    magma_minproduct_cmfree( &b_h, queue );
    
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
    A           magma_minproduct_c_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_c_matrix
                RHS b

    @param[in]
    M           magma_minproduct_c_matrix*
                M = D^(-1) * (L+U)

    @param[in]
    c           magma_minproduct_c_matrix*
                c = D^(-1) * b
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cjacobisetup(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b,
    magma_minproduct_c_matrix *M, magma_minproduct_c_matrix *c,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_int_t i;

    magma_minproduct_c_matrix A_h1={Magma_minproduct_CSR}, A_h2={Magma_minproduct_CSR}, B={Magma_minproduct_CSR}, C={Magma_minproduct_CSR};
    magma_minproduct_c_matrix diag={Magma_minproduct_CSR}, c_t={Magma_minproduct_CSR}, b_h={Magma_minproduct_CSR};
    CHECK( magma_minproduct_cvinit( &c_t, Magma_minproduct_CPU, A.num_rows, b.num_cols, MAGMA_minproduct_C_ZERO, queue ));
    CHECK( magma_minproduct_cvinit( &diag, Magma_minproduct_CPU, A.num_rows, b.num_cols, MAGMA_minproduct_C_ZERO, queue ));
    CHECK( magma_minproduct_cmtransfer( b, &b_h, A.memory_location, Magma_minproduct_CPU, queue ));

    if ( A.storage_type != Magma_minproduct_CSR ) {
        CHECK( magma_minproduct_cmtransfer( A, &A_h1, A.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_cmconvert( A_h1, &B, A_h1.storage_type, Magma_minproduct_CSR, queue ));
    }
    else {
        CHECK( magma_minproduct_cmtransfer( A, &B, A.memory_location, Magma_minproduct_CPU, queue ));
    }
    for( magma_minproduct_int_t rowindex=0; rowindex<B.num_rows; rowindex++ ) {
        magma_minproduct_int_t start = (B.drow[rowindex]);
        magma_minproduct_int_t end = (B.drow[rowindex+1]);
        for( i=start; i<end; i++ ) {
            if ( B.dcol[i]==rowindex ) {
                diag.val[rowindex] = B.val[i];
                if ( MAGMA_minproduct_C_REAL( diag.val[rowindex]) == 0 )
                    printf(" error: zero diagonal element in row %d!\n",
                                                               (int) rowindex);
            }
        }
        for( i=start; i<end; i++ ) {
            B.val[i] = B.val[i] / diag.val[rowindex];
            if ( B.dcol[i]==rowindex ) {
                B.val[i] = MAGMA_minproduct_C_MAKE( 0., 0. );
            }
        }
        c_t.val[rowindex] = b_h.val[rowindex] / diag.val[rowindex];

    }

    CHECK( magma_minproduct_c_csr_compressor(&B.val, &B.drow, &B.dcol,
                           &C.val, &C.drow, &C.dcol, &B.num_rows, queue ));

    C.num_rows = B.num_rows;
    C.num_cols = B.num_cols;
    C.memory_location = B.memory_location;
    C.nnz = C.drow[B.num_rows];
    C.storage_type = B.storage_type;
    C.memory_location = B.memory_location;
    if ( A.storage_type != Magma_minproduct_CSR) {
        A_h2.alignment = A.alignment;
        A_h2.blocksize = A.blocksize;
        CHECK( magma_minproduct_cmconvert( C, &A_h2, Magma_minproduct_CSR, A_h1.storage_type, queue ));
        CHECK( magma_minproduct_cmtransfer( A_h2, M, Magma_minproduct_CPU, A.memory_location, queue ));
    }
    else {
        CHECK( magma_minproduct_cmtransfer( C, M, Magma_minproduct_CPU, A.memory_location, queue ));
    }
    CHECK( magma_minproduct_cmtransfer( c_t, c, Magma_minproduct_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        magma_minproduct_cmfree( &A_h1, queue );
        magma_minproduct_cmfree( &A_h2, queue );
    }
    
cleanup:
    magma_minproduct_cmfree( &B, queue );
    magma_minproduct_cmfree( &C, queue );
    magma_minproduct_cmfree( &diag, queue );
    magma_minproduct_cmfree( &c_t, queue );
    magma_minproduct_cmfree( &b_h, queue );

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
    M           magma_minproduct_c_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_minproduct_c_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_minproduct_c_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_minproduct_c_solver_par*
                solver parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cjacobiiter(
    magma_minproduct_c_matrix M, magma_minproduct_c_matrix c, magma_minproduct_c_matrix *x,
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // local variables
    magma_minproductFloatComplex c_zero = MAGMA_minproduct_C_ZERO, c_one = MAGMA_minproduct_C_ONE,
                                            c_mone = MAGMA_minproduct_C_NEG_ONE;
    magma_minproduct_int_t dofs = M.num_rows*x->num_cols;
    magma_minproduct_c_matrix t={Magma_minproduct_CSR}, swap={Magma_minproduct_CSR};
    CHECK( magma_minproduct_cvinit( &t, Magma_minproduct_DEV, M.num_rows, x->num_cols, c_zero, queue ));


    for( magma_minproduct_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_minproduct_c_spmv( c_mone, M, *x, c_zero, t, queue ));        // t = - M * x
        magma_minproduct_caxpy( dofs, c_one , c.dval, 1 , t.dval, 1 );        // t = t + c

        // swap so that x again contains solution, and y is ready to be used
        swap = *x;
        *x = t;
        t = swap;
        //magma_minproduct_ccopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }

cleanup:
    magma_minproduct_cmfree( &t, queue );

    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_cjacobiiter */



/**
    Purpose
    -------

    Iterates the solution approximation according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    Arguments
    ---------

    @param[in]
    M           magma_minproduct_c_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_minproduct_c_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_minproduct_c_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_minproduct_c_solver_par*
                solver parameters

    @param[in]
    solver_par  magma_minproduct_c_precond_par*
                precond parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cjacobiiter_precond(
    magma_minproduct_c_matrix M, magma_minproduct_c_matrix *x,
    magma_minproduct_c_solver_par *solver_par, magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // local variables
    magma_minproductFloatComplex c_zero = MAGMA_minproduct_C_ZERO, c_one = MAGMA_minproduct_C_ONE,
                                            c_mone = MAGMA_minproduct_C_NEG_ONE;
    magma_minproduct_int_t dofs = M.num_rows;
    magma_minproduct_int_t num_vecs = x->num_rows / dofs;
    magma_minproduct_c_matrix swap={Magma_minproduct_CSR};

    for( magma_minproduct_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_minproduct_c_spmv( c_mone, M, *x, c_zero, precond->work2, queue )); // t = - M * x

        magma_minproduct_caxpy( num_vecs*dofs, c_one ,
                precond->work1.dval, 1 , precond->work2.dval, 1 ); // t = t + c

        // swap so that x again contains solution, and y is ready to be used
        swap = *x;
        *x = precond->work2;
        precond->work2 = swap;
        //magma_minproduct_ccopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }
    
cleanup:
    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_cjacobiiter */



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
    M           magma_minproduct_c_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_minproduct_c_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_minproduct_c_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_minproduct_c_solver_par*
                solver parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cjacobiiter_sys(
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix b,
    magma_minproduct_c_matrix d,
    magma_minproduct_c_matrix t,
    magma_minproduct_c_matrix *x,
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // local variables
    magma_minproductFloatComplex c_zero = MAGMA_minproduct_C_ZERO, c_one = MAGMA_minproduct_C_ONE;

    for( magma_minproduct_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_minproduct_c_spmv( c_one, A, *x, c_zero, t, queue ));        // t =  A * x
        CHECK( magma_minproduct_cjacobiupdate( t, b, d, x, queue ));
        // swap so that x again contains solution, and y is ready to be used
        //swap = *x;
        //*x = t;
        //t = swap;
        //magma_minproduct_ccopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }
    
cleanup:
    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_cjacobiiter_sys */
