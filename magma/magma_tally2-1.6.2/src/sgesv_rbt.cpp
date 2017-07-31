/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgesv_rbt.cpp normal z -> s, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_tally2.h"

#define BWDMAX 1.0
#define ITERMAX 30

/**
    Purpose
    -------
    SGERFS  improve the computed solution to a system of linear
          equations.

        
    The iterative refinement process is stopped if
        ITER > ITERMAX
    or for all the RHS we have:
        RNRM < SQRT(n)*XNRM*ANRM*EPS*BWDMAX
    where
        o ITER is the number of the current iteration in the iterative
          refinement process
        o RNRM is the infinity-norm of the residual
        o XNRM is the infinity-norm of the solution
        o ANRM is the infinity-operator-norm of the matrix A
        o EPS is the machine epsilon returned by SLAMCH('Epsilon')
    The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00 respectively.

    Arguments
    ---------
    @param[in]
    trans   magma_tally2_trans_t
            Specifies the form of the system of equations:
      -     = Magma_tally2NoTrans:    A    * X = B  (No transpose)
      -     = Magma_tally2Trans:      A**T * X = B  (Transpose)
      -     = Magma_tally2ConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The number of linear equations, i.e., the order of the
            matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      REAL array on the GPU, dimension (ldda,N)
            the N-by-N coefficient matrix A.
            
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  ldda >= max(1,N).

    @param[in]
    dB      REAL array on the GPU, dimension (lddb,NRHS)
            The N-by-NRHS right hand side matrix B.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  lddb >= max(1,N).

    @param[in, out]
    dX      REAL array on the GPU, dimension (lddx,NRHS)
            On entry, the solution matrix X, as computed by
            SGETRS_NOPIV.  On exit, the improved solution matrix X.

    @param[in]
    lddx    INTEGER
            The leading dimension of the array dX.  lddx >= max(1,N).

    @param
    dworkd  (workspace) REAL array on the GPU, dimension (N*NRHS)
            This array is used to hold the residual vectors.

    @param
    dAF     REAL array on the GPU, dimension (ldda,n)
            The factors L and U from the factorization A = L*U
            as computed by SGETRF_NOPIV.

    @param[out]
    iter    INTEGER
      -     < 0: iterative refinement has failed, real
                 factorization has been performed
        +        -1 : the routine fell back to full precision for
                      implementation- or machine-specific reasons
        +        -2 : narrowing the precision induced an overflow,
                      the routine fell back to full precision
        +        -3 : failure of SGETRF
        +        -31: stop the iterative refinement after the 30th iteration
      -     > 0: iterative refinement has been successfully used.
                 Returns the number of iterations
 
    @param[out]
    info   INTEGER
      -     = 0:  successful exit
      -     < 0:  if info = -i, the i-th argument had an illegal value
      -     > 0:  if info = i, U(i,i) computed in REAL is
                  exactly zero.  The factorization has been completed,
                  but the factor U is exactly singular, so the solution
                  could not be computed.

    @ingroup magma_tally2_sgesv_driver
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_sgerfs_nopiv_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    magma_tally2Float_ptr dX, magma_tally2_int_t lddx,
    magma_tally2Float_ptr dworkd, magma_tally2Float_ptr dAF,
    magma_tally2_int_t *iter,
    magma_tally2_int_t *info)
{
    #define dB(i,j)     (dB + (i) + (j)*lddb)
    #define dX(i,j)     (dX + (i) + (j)*lddx)
    #define dR(i,j)     (dR + (i) + (j)*lddr)
    
    float c_neg_one = MAGMA_tally2_S_NEG_ONE;
    float c_one     = MAGMA_tally2_S_ONE;
    magma_tally2_int_t     ione  = 1;
    magma_tally2Float_ptr dR;
    float Xnrmv, Rnrmv;
    float          Anrm, Xnrm, Rnrm, cte, eps;
    magma_tally2_int_t     i, j, iiter, lddsa, lddr;
    
    /* Check arguments */
    *iter = 0;
    *info = 0;
    if ( n < 0 )
        *info = -1;
    else if ( nrhs < 0 )
        *info = -2;
    else if ( ldda < max(1,n))
        *info = -4;
    else if ( lddb < max(1,n))
        *info = -8;
    else if ( lddx < max(1,n))
        *info = -10;
    
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }
    
    if ( n == 0 || nrhs == 0 )
        return *info;

    lddsa = n;
    lddr  = n;
    
    dR  = dworkd;
    
    eps  = lapackf77_slamch("Epsilon");
    Anrm = magma_tally2blas_slange(Magma_tally2InfNorm, n, n, dA, ldda, (float*)dworkd );
    cte  = Anrm * eps * pow( (float)n, (float)0.5 ) * BWDMAX;
    
    // residual dR = dB - dA*dX in real
    magma_tally2blas_slacpy( Magma_tally2UpperLower, n, nrhs, dB, lddb, dR, lddr );
    if ( nrhs == 1 ) {
        magma_tally2_sgemv( trans, n, n,
                     c_neg_one, dA, ldda,
                                dX, 1,
                     c_one,     dR, 1 );
    }
    else {
        magma_tally2_sgemm( trans, Magma_tally2NoTrans, n, nrhs, n,
                     c_neg_one, dA, ldda,
                                dX, lddx,
                     c_one,     dR, lddr );
    }
    
    // TODO: use MAGMA_tally2_S_ABS( dX(i,j) ) instead of slange?
    for( j=0; j < nrhs; j++ ) {
        i = magma_tally2_isamax( n, dX(0,j), 1) - 1;
        magma_tally2_sgetmatrix( 1, 1, dX(i,j), 1, &Xnrmv, 1 );
        Xnrm = lapackf77_slange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
        
        i = magma_tally2_isamax ( n, dR(0,j), 1 ) - 1;
        magma_tally2_sgetmatrix( 1, 1, dR(i,j), 1, &Rnrmv, 1 );
        Rnrm = lapackf77_slange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
       


 //       printf("Rnrm : %e, Xnrm*cte : %e\n", Rnrm, Xnrm*cte);



        if ( Rnrm >  Xnrm*cte ) {
            goto REFINEMENT;
        }
    }
    
    *iter = 0;
    return *info;

REFINEMENT:
    for( iiter=1; iiter < ITERMAX; ) {
        *info = 0;
        // solve dAF*dX = dR 
        // it's okay that dR is used for both dB input and dX output.
        magma_tally2_sgetrs_nopiv_gpu( trans, n, nrhs, dAF, lddsa, dR, lddr, info );
        if (*info != 0) {
            *iter = -3;
            goto FALLBACK;
        }
        
        // Add correction and setup residual
        // dX += dR  --and--
        // dR = dB
        // This saves going through dR a second time (if done with one more kernel).
        // -- not really: first time is read, second time is write.
        for( j=0; j < nrhs; j++ ) {
            magma_tally2blas_saxpycp2( n, dR(0,j), dX(0,j), dB(0,j) );
        }
        
        // residual dR = dB - dA*dX in real
        if ( nrhs == 1 ) {
            magma_tally2_sgemv( trans, n, n,
                         c_neg_one, dA, ldda,
                                    dX, 1,
                         c_one,     dR, 1 );
        }
        else {
            magma_tally2_sgemm( trans, Magma_tally2NoTrans, n, nrhs, n,
                         c_neg_one, dA, ldda,
                                    dX, lddx,
                         c_one,     dR, lddr );
        }
        
        /*  Check whether the nrhs normwise backward errors satisfy the
         *  stopping criterion. If yes, set ITER=IITER > 0 and return. */
        for( j=0; j < nrhs; j++ ) {
            i = magma_tally2_isamax( n, dX(0,j), 1) - 1;
            magma_tally2_sgetmatrix( 1, 1, dX(i,j), 1, &Xnrmv, 1 );
            Xnrm = lapackf77_slange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
            
            i = magma_tally2_isamax ( n, dR(0,j), 1 ) - 1;
            magma_tally2_sgetmatrix( 1, 1, dR(i,j), 1, &Rnrmv, 1 );
            Rnrm = lapackf77_slange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
            
            if ( Rnrm >  Xnrm*cte ) {
                goto L20;
            }
        }
        
        /*  If we are here, the nrhs normwise backward errors satisfy
         *  the stopping criterion, we are good to exit. */
        *iter = iiter;
        return *info;
        
      L20:
        iiter++;
    }

    
    /* If we are at this place of the code, this is because we have
     * performed ITER=ITERMAX iterations and never satisified the
     * stopping criterion. Set up the ITER flag accordingly. */
    *iter = -ITERMAX - 1;
    
FALLBACK:
    /* Iterative refinement failed to converge to a
     * satisfactory solution. */
    
    return *info;
}











/**
    Purpose
    -------
    Solves a system of linear equations
       A * X = B
    where A is a general N-by-N matrix and X and B are N-by-NRHS matrices.
    Random Butterfly Tranformation is applied on A and B, then 
    the LU decomposition with no pivoting is
    used to factor A as
       A = L * U,
    where L is unit lower triangular, and U is
    upper triangular.  The factored form of A is then used to solve the
    system of equations A * X = B.
    The solution can then be improved using iterative refinement.


    Arguments
    ---------
    @param[in]
    ref  magma_tally2_bool_t
    Specifies if iterative refinement have to be applied to improve the solution.
    -     = Magma_tally2True:   Iterative refinement is applied.
    -     = Magma_tally2False:  Iterative refinement is not applied.

    @param[in]
    n       INTEGER
    The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
    The number of right hand sides, i.e., the number of columns
    of the matrix B.  NRHS >= 0.

    @param[in,out]
    A       REAL array, dimension (LDA,N).
    On entry, the M-by-N matrix to be factored.
    On exit, the factors L and U from the factorization
    A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    lda     INTEGER
    The leading dimension of the array A.  LDA >= max(1,N).


    @param[in,out]
    B       REAL array, dimension (LDB,NRHS)
    On entry, the right hand side matrix B.
    On exit, the solution matrix X.

    @param[in]
    ldb     INTEGER
    The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
    -     = 0:  successful exit
    -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_sgesv_driver
 ********************************************************************/


extern "C" magma_tally2_int_t 
magma_tally2_sgesv_rbt(
    magma_tally2_bool_t ref, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    float *A, magma_tally2_int_t lda, 
    float *B, magma_tally2_int_t ldb, 
    magma_tally2_int_t *info)
{

    /* Function Body */
    *info = 0;
    if ( ! (ref == Magma_tally2True) &&
         ! (ref == Magma_tally2False) ) {
        *info = -1;
    }
    else if (n < 0) {
        *info = -2;
    } else if (nrhs < 0) {
        *info = -3;
    } else if (lda < max(1,n)) {
        *info = -5;
    } else if (ldb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (nrhs == 0 || n == 0)
        return *info;


    magma_tally2_int_t nn = n + ((4-(n % 4))%4);
    float *dA, *hu, *hv, *db, *dAo, *dBo, *dwork;
    magma_tally2_int_t n2;

    magma_tally2_int_t iter;
    n2 = nn*nn;

    if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc( &dA, n2 )) {
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }
    if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc( &db, nn*nrhs )) {
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }

    if (ref == Magma_tally2True) {
        if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc( &dAo, n2 )) {
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }
        if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc( &dwork, nn*nrhs )) {
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }
        if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc( &dBo, nn*nrhs )) {
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }
    }

    if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc_cpu( &hu, 2*nn )) {
        *info = MAGMA_tally2_ERR_HOST_ALLOC;
        return *info;
    }
    if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc_cpu( &hv, 2*nn )) {
        *info = MAGMA_tally2_ERR_HOST_ALLOC;
        return *info;
    }

    magma_tally2blas_slaset(Magma_tally2Full, nn, nn, MAGMA_tally2_S_ZERO, MAGMA_tally2_S_ONE, dA, nn);

    /* Send matrix on the GPU*/
    magma_tally2_ssetmatrix(n, n, A, lda, dA, nn);

    /* Send b on the GPU*/
    magma_tally2_ssetmatrix(n, nrhs, B, ldb, db, nn);

    *info = magma_tally2_sgerbt_gpu(Magma_tally2True, nn, nrhs, dA, nn, db, nn, hu, hv, info);
    if (*info != MAGMA_tally2_SUCCESS)  {
        return *info;
    }

    if (ref == Magma_tally2True) {
        magma_tally2_scopymatrix(nn, nn, dA, nn, dAo, nn);
        magma_tally2_scopymatrix(nn, nrhs, db, nn, dBo, nn);
    }
    /* Solve the system U^TAV.y = U^T.b on the GPU*/ 
    magma_tally2_sgesv_nopiv_gpu( nn, nrhs, dA, nn, db, nn, info);


    /* Iterative refinement */
    if (ref == Magma_tally2True) {
        magma_tally2_sgerfs_nopiv_gpu(Magma_tally2NoTrans, nn, nrhs, dAo, nn, dBo, nn, db, nn, dwork, dA, &iter, info);
    }
    //printf("iter = %d\n", iter);

    /* The solution of A.x = b is Vy computed on the GPU */
    float *dv;

    if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc( &dv, 2*nn )) {
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }

    magma_tally2_ssetvector(2*nn, hv, 1, dv, 1);
    
    for(int i = 0; i < nrhs; i++) {
        magma_tally2blas_sprbt_mv(nn, dv, db+(i*nn));
    }

    magma_tally2_sgetmatrix(n, nrhs, db, nn, B, ldb);

    magma_tally2_free_cpu( hu);
    magma_tally2_free_cpu( hv);

    magma_tally2_free( dA );
    magma_tally2_free( dv );
    magma_tally2_free( db );
    
    if (ref == Magma_tally2True) {    
        magma_tally2_free( dAo );
        magma_tally2_free( dBo );
        magma_tally2_free( dwork );
    }
    return *info;
}
