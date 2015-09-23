/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zjacobi.cpp normal z -> d, Sun May  3 11:22:59 2015
*/

#include "common_magma_minproductsparse.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


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
    A           magma_minproduct_d_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                RHS b

    @param[in,out]
    x           magma_minproduct_d_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_minproduct_d_solver_par*
                solver parameters

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_dgesv
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_djacobi(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    

    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );
    
    // some useful variables
    double c_zero = MAGMA_minproduct_D_ZERO;
    
    double nom0 = 0.0;

    magma_minproduct_d_matrix r={Magma_minproduct_CSR}, d={Magma_minproduct_CSR}, ACSR={Magma_minproduct_CSR} ;
    
    CHECK( magma_minproduct_dmconvert(A, &ACSR, A.storage_type, Magma_minproduct_CSR, queue ) );

    // prepare solver feedback
    solver_par->solver = Magma_minproduct_JACOBI;
    solver_par->info = MAGMA_minproduct_SUCCESS;

    real_Double_t tempo1, tempo2;
    double residual;
    // solver setup
    CHECK( magma_minproduct_dvinit( &r, Magma_minproduct_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK(  magma_minproduct_dresidualvec( ACSR, b, *x, &r, &residual, queue));
    solver_par->init_res = residual;
    solver_par->res_vec = NULL;
    solver_par->timing = NULL;
    nom0 = residual;

    // Jacobi setup
    CHECK( magma_minproduct_djacobisetup_diagscal( ACSR, &d, queue ));
    magma_minproduct_d_solver_par jacobiiter_par;
    jacobiiter_par.maxiter = solver_par->maxiter;

    tempo1 = magma_minproduct_sync_wtime( queue );

    // Jacobi iterator
    CHECK( magma_minproduct_djacobispmvupdate(jacobiiter_par.maxiter, ACSR, r, b, d, x, queue ));

    tempo2 = magma_minproduct_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    CHECK(  magma_minproduct_dresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->final_res = residual;
    solver_par->numiter = solver_par->maxiter;

    if ( solver_par->init_res > solver_par->final_res )
        info = MAGMA_minproduct_SUCCESS;
    else
        info = MAGMA_minproduct_DIVERGENCE;
    
cleanup:
    magma_minproduct_dmfree( &r, queue );
    magma_minproduct_dmfree( &d, queue );
    magma_minproduct_dmfree( &ACSR, queue );

    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_djacobi */






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
    A           magma_minproduct_d_matrix
                input matrix A

    @param[in]
    M           magma_minproduct_d_matrix*
                M = D^(-1) * (L+U)

    @param[in,out]
    d           magma_minproduct_d_matrix*
                vector with diagonal elements of A
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_d
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_djacobisetup_matrix(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix *M, magma_minproduct_d_matrix *d,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_int_t i;

    magma_minproduct_d_matrix A_h1={Magma_minproduct_CSR}, A_h2={Magma_minproduct_CSR}, B={Magma_minproduct_CSR}, C={Magma_minproduct_CSR};
    magma_minproduct_d_matrix diag={Magma_minproduct_CSR};
    CHECK( magma_minproduct_dvinit( &diag, Magma_minproduct_CPU, A.num_rows, 1, MAGMA_minproduct_D_ZERO, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        CHECK( magma_minproduct_dmtransfer( A, &A_h1, A.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_dmconvert( A_h1, &B, A_h1.storage_type, Magma_minproduct_CSR, queue ));
    }
    else {
        CHECK( magma_minproduct_dmtransfer( A, &B, A.memory_location, Magma_minproduct_CPU, queue ));
    }
    for( magma_minproduct_int_t rowindex=0; rowindex<B.num_rows; rowindex++ ) {
        magma_minproduct_int_t start = (B.drow[rowindex]);
        magma_minproduct_int_t end = (B.drow[rowindex+1]);
        for( i=start; i<end; i++ ) {
            if ( B.dcol[i]==rowindex ) {
                diag.val[rowindex] = B.val[i];
                if ( MAGMA_minproduct_D_REAL( diag.val[rowindex]) == 0 )
                    printf(" error: zero diagonal element in row %d!\n",
                                                                (int) rowindex);
            }
        }
        for( i=start; i<end; i++ ) {
            B.val[i] = B.val[i] / diag.val[rowindex];
            if ( B.dcol[i]==rowindex ) {
                B.val[i] = MAGMA_minproduct_D_MAKE( 0., 0. );
            }
        }
    }
    CHECK( magma_minproduct_d_csr_compressor(&B.val, &B.drow, &B.dcol,
                           &C.val, &C.drow, &C.dcol, &B.num_rows, queue ));
    C.num_rows = B.num_rows;
    C.num_cols = B.num_cols;
    C.memory_location = B.memory_location;
    C.nnz = C.drow[B.num_rows];
    C.storage_type = B.storage_type;
    C.memory_location = B.memory_location;
    if ( A.storage_type != Magma_minproduct_CSR) {
        CHECK( magma_minproduct_dmconvert( C, &A_h2, Magma_minproduct_CSR, A_h1.storage_type, queue ));
        CHECK( magma_minproduct_dmtransfer( A_h2, M, Magma_minproduct_CPU, A.memory_location, queue ));
    }
    else {
        CHECK( magma_minproduct_dmtransfer( C, M, Magma_minproduct_CPU, A.memory_location, queue ));
    }
    CHECK( magma_minproduct_dmtransfer( diag, d, Magma_minproduct_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        magma_minproduct_dmfree( &A_h1, queue );
        magma_minproduct_dmfree( &A_h2, queue );
    }
    
cleanup:
    magma_minproduct_dmfree( &B, queue );
    magma_minproduct_dmfree( &C, queue );

    magma_minproduct_dmfree( &diag, queue );
 
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
    A           magma_minproduct_d_matrix
                input matrix A

    @param[in,out]
    d           magma_minproduct_d_matrix*
                vector with diagonal elements
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_d
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_djacobisetup_diagscal(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix *d,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_int_t i;

    magma_minproduct_d_matrix A_h1={Magma_minproduct_CSR}, B={Magma_minproduct_CSR};
    magma_minproduct_d_matrix diag={Magma_minproduct_CSR};
    CHECK( magma_minproduct_dvinit( &diag, Magma_minproduct_CPU, A.num_rows, 1, MAGMA_minproduct_D_ZERO, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        CHECK( magma_minproduct_dmtransfer( A, &A_h1, A.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_dmconvert( A_h1, &B, A_h1.storage_type, Magma_minproduct_CSR, queue ));
    }
    else {
        CHECK( magma_minproduct_dmtransfer( A, &B, A.memory_location, Magma_minproduct_CPU, queue ));
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
        if ( diag.val[rowindex] == MAGMA_minproduct_D_ZERO ){
            printf(" error: zero diagonal element in row %d!\n",
                                                        (int) rowindex);
            
            if ( A.storage_type != Magma_minproduct_CSR) {
                magma_minproduct_dmfree( &A_h1, queue );
            }
            magma_minproduct_dmfree( &B, queue );
            magma_minproduct_dmfree( &diag, queue );
            info = MAGMA_minproduct_ERR_BADPRECOND;
            goto cleanup;
        }
    }
    CHECK( magma_minproduct_dmtransfer( diag, d, Magma_minproduct_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        magma_minproduct_dmfree( &A_h1, queue );
    }
    
cleanup:
    magma_minproduct_dmfree( &B, queue );
    magma_minproduct_dmfree( &diag, queue );
 
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
    b           magma_minproduct_d_matrix
                RHS b

    @param[in]
    d           magma_minproduct_d_matrix
                vector with diagonal entries

    @param[in]
    c           magma_minproduct_d_matrix*
                c = D^(-1) * b
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_d
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_djacobisetup_vector(
    magma_minproduct_d_matrix b, magma_minproduct_d_matrix d,
    magma_minproduct_d_matrix *c,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_d_matrix diag={Magma_minproduct_CSR}, c_t={Magma_minproduct_CSR}, b_h={Magma_minproduct_CSR}, tmp={Magma_minproduct_CSR};
    
    if ( b.memory_location == Magma_minproduct_CPU ) {

        CHECK( magma_minproduct_dvinit( &c_t, Magma_minproduct_CPU, b.num_rows, b.num_cols, MAGMA_minproduct_D_ZERO, queue ));

        CHECK( magma_minproduct_dmtransfer( b, &b_h, b.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_dmtransfer( d, &diag, b.memory_location, Magma_minproduct_CPU, queue ));

        for( magma_minproduct_int_t rowindex=0; rowindex<b.num_rows; rowindex++ ) {
            c_t.val[rowindex] = b_h.val[rowindex] / diag.val[rowindex];

        }
        CHECK( magma_minproduct_dmtransfer( c_t, c, Magma_minproduct_CPU, b.memory_location, queue ));
    }
    else if ( b.memory_location == Magma_minproduct_DEV ) {
        // fill vector
        CHECK( magma_minproduct_dvinit( &tmp, Magma_minproduct_DEV, b.num_rows, b.num_cols, MAGMA_minproduct_D_ZERO, queue ));
        CHECK( magma_minproduct_djacobisetup_vector_gpu(
                    b.num_rows, b, d, *c, &tmp, queue ));
        goto cleanup;
    }

cleanup:
    magma_minproduct_dmfree( &tmp, queue );
    magma_minproduct_dmfree( &diag, queue );
    magma_minproduct_dmfree( &c_t, queue );
    magma_minproduct_dmfree( &b_h, queue );
    
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
    A           magma_minproduct_d_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                RHS b

    @param[in]
    M           magma_minproduct_d_matrix*
                M = D^(-1) * (L+U)

    @param[in]
    c           magma_minproduct_d_matrix*
                c = D^(-1) * b
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_d
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_djacobisetup(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix *M, magma_minproduct_d_matrix *c,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_int_t i;

    magma_minproduct_d_matrix A_h1={Magma_minproduct_CSR}, A_h2={Magma_minproduct_CSR}, B={Magma_minproduct_CSR}, C={Magma_minproduct_CSR};
    magma_minproduct_d_matrix diag={Magma_minproduct_CSR}, c_t={Magma_minproduct_CSR}, b_h={Magma_minproduct_CSR};
    CHECK( magma_minproduct_dvinit( &c_t, Magma_minproduct_CPU, A.num_rows, b.num_cols, MAGMA_minproduct_D_ZERO, queue ));
    CHECK( magma_minproduct_dvinit( &diag, Magma_minproduct_CPU, A.num_rows, b.num_cols, MAGMA_minproduct_D_ZERO, queue ));
    CHECK( magma_minproduct_dmtransfer( b, &b_h, A.memory_location, Magma_minproduct_CPU, queue ));

    if ( A.storage_type != Magma_minproduct_CSR ) {
        CHECK( magma_minproduct_dmtransfer( A, &A_h1, A.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_dmconvert( A_h1, &B, A_h1.storage_type, Magma_minproduct_CSR, queue ));
    }
    else {
        CHECK( magma_minproduct_dmtransfer( A, &B, A.memory_location, Magma_minproduct_CPU, queue ));
    }
    for( magma_minproduct_int_t rowindex=0; rowindex<B.num_rows; rowindex++ ) {
        magma_minproduct_int_t start = (B.drow[rowindex]);
        magma_minproduct_int_t end = (B.drow[rowindex+1]);
        for( i=start; i<end; i++ ) {
            if ( B.dcol[i]==rowindex ) {
                diag.val[rowindex] = B.val[i];
                if ( MAGMA_minproduct_D_REAL( diag.val[rowindex]) == 0 )
                    printf(" error: zero diagonal element in row %d!\n",
                                                               (int) rowindex);
            }
        }
        for( i=start; i<end; i++ ) {
            B.val[i] = B.val[i] / diag.val[rowindex];
            if ( B.dcol[i]==rowindex ) {
                B.val[i] = MAGMA_minproduct_D_MAKE( 0., 0. );
            }
        }
        c_t.val[rowindex] = b_h.val[rowindex] / diag.val[rowindex];

    }

    CHECK( magma_minproduct_d_csr_compressor(&B.val, &B.drow, &B.dcol,
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
        CHECK( magma_minproduct_dmconvert( C, &A_h2, Magma_minproduct_CSR, A_h1.storage_type, queue ));
        CHECK( magma_minproduct_dmtransfer( A_h2, M, Magma_minproduct_CPU, A.memory_location, queue ));
    }
    else {
        CHECK( magma_minproduct_dmtransfer( C, M, Magma_minproduct_CPU, A.memory_location, queue ));
    }
    CHECK( magma_minproduct_dmtransfer( c_t, c, Magma_minproduct_CPU, A.memory_location, queue ));

    if ( A.storage_type != Magma_minproduct_CSR) {
        magma_minproduct_dmfree( &A_h1, queue );
        magma_minproduct_dmfree( &A_h2, queue );
    }
    
cleanup:
    magma_minproduct_dmfree( &B, queue );
    magma_minproduct_dmfree( &C, queue );
    magma_minproduct_dmfree( &diag, queue );
    magma_minproduct_dmfree( &c_t, queue );
    magma_minproduct_dmfree( &b_h, queue );

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
    M           magma_minproduct_d_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_minproduct_d_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_minproduct_d_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_minproduct_d_solver_par*
                solver parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_d
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_djacobiiter(
    magma_minproduct_d_matrix M, magma_minproduct_d_matrix c, magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // local variables
    double c_zero = MAGMA_minproduct_D_ZERO, c_one = MAGMA_minproduct_D_ONE,
                                            c_mone = MAGMA_minproduct_D_NEG_ONE;
    magma_minproduct_int_t dofs = M.num_rows*x->num_cols;
    magma_minproduct_d_matrix t={Magma_minproduct_CSR}, swap={Magma_minproduct_CSR};
    CHECK( magma_minproduct_dvinit( &t, Magma_minproduct_DEV, M.num_rows, x->num_cols, c_zero, queue ));


    for( magma_minproduct_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_minproduct_d_spmv( c_mone, M, *x, c_zero, t, queue ));        // t = - M * x
        magma_minproduct_daxpy( dofs, c_one , c.dval, 1 , t.dval, 1 );        // t = t + c

        // swap so that x again contains solution, and y is ready to be used
        swap = *x;
        *x = t;
        t = swap;
        //magma_minproduct_dcopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }

cleanup:
    magma_minproduct_dmfree( &t, queue );

    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_djacobiiter */



/**
    Purpose
    -------

    Iterates the solution approximation according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    Arguments
    ---------

    @param[in]
    M           magma_minproduct_d_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_minproduct_d_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_minproduct_d_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_minproduct_d_solver_par*
                solver parameters

    @param[in]
    solver_par  magma_minproduct_d_precond_par*
                precond parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_d
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_djacobiiter_precond(
    magma_minproduct_d_matrix M, magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par, magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // local variables
    double c_zero = MAGMA_minproduct_D_ZERO, c_one = MAGMA_minproduct_D_ONE,
                                            c_mone = MAGMA_minproduct_D_NEG_ONE;
    magma_minproduct_int_t dofs = M.num_rows;
    magma_minproduct_int_t num_vecs = x->num_rows / dofs;
    magma_minproduct_d_matrix swap={Magma_minproduct_CSR};

    for( magma_minproduct_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_minproduct_d_spmv( c_mone, M, *x, c_zero, precond->work2, queue )); // t = - M * x

        magma_minproduct_daxpy( num_vecs*dofs, c_one ,
                precond->work1.dval, 1 , precond->work2.dval, 1 ); // t = t + c

        // swap so that x again contains solution, and y is ready to be used
        swap = *x;
        *x = precond->work2;
        precond->work2 = swap;
        //magma_minproduct_dcopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }
    
cleanup:
    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_djacobiiter */



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
    M           magma_minproduct_d_matrix
                input matrix M = D^(-1) * (L+U)

    @param[in]
    c           magma_minproduct_d_matrix
                c = D^(-1) * b

    @param[in,out]
    x           magma_minproduct_d_matrix*
                iteration vector x

    @param[in,out]
    solver_par  magma_minproduct_d_solver_par*
                solver parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_d
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_djacobiiter_sys(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix d,
    magma_minproduct_d_matrix t,
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // local variables
    double c_zero = MAGMA_minproduct_D_ZERO, c_one = MAGMA_minproduct_D_ONE;

    for( magma_minproduct_int_t i=0; i<solver_par->maxiter; i++ ) {
        CHECK( magma_minproduct_d_spmv( c_one, A, *x, c_zero, t, queue ));        // t =  A * x
        CHECK( magma_minproduct_djacobiupdate( t, b, d, x, queue ));
        // swap so that x again contains solution, and y is ready to be used
        //swap = *x;
        //*x = t;
        //t = swap;
        //magma_minproduct_dcopy( dofs, t.dval, 1 , x->dval, 1 );               // x = t
    }
    
cleanup:
    magma_minproductblasSetKernelStream( orig_queue );
    solver_par->info = info;
    return info;
}   /* magma_minproduct_djacobiiter_sys */
