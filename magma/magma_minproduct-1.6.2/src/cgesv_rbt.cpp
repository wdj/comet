/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgesv_rbt.cpp normal z -> c, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_minproduct.h"

#define BWDMAX 1.0
#define ITERMAX 30

/**
    Purpose
    -------
    CGERFS  improve the computed solution to a system of linear
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
    trans   magma_minproduct_trans_t
            Specifies the form of the system of equations:
      -     = Magma_minproductNoTrans:    A    * X = B  (No transpose)
      -     = Magma_minproductTrans:      A**T * X = B  (Transpose)
      -     = Magma_minproductConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The number of linear equations, i.e., the order of the
            matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX array on the GPU, dimension (ldda,N)
            the N-by-N coefficient matrix A.
            
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  ldda >= max(1,N).

    @param[in]
    dB      COMPLEX array on the GPU, dimension (lddb,NRHS)
            The N-by-NRHS right hand side matrix B.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  lddb >= max(1,N).

    @param[in, out]
    dX      COMPLEX array on the GPU, dimension (lddx,NRHS)
            On entry, the solution matrix X, as computed by
            CGETRS_NOPIV.  On exit, the improved solution matrix X.

    @param[in]
    lddx    INTEGER
            The leading dimension of the array dX.  lddx >= max(1,N).

    @param
    dworkd  (workspace) COMPLEX array on the GPU, dimension (N*NRHS)
            This array is used to hold the residual vectors.

    @param
    dAF     COMPLEX array on the GPU, dimension (ldda,n)
            The factors L and U from the factorization A = L*U
            as computed by CGETRF_NOPIV.

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

    @ingroup magma_minproduct_cgesv_driver
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_cgerfs_nopiv_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductFloatComplex_ptr dworkd, magma_minproductFloatComplex_ptr dAF,
    magma_minproduct_int_t *iter,
    magma_minproduct_int_t *info)
{
    #define dB(i,j)     (dB + (i) + (j)*lddb)
    #define dX(i,j)     (dX + (i) + (j)*lddx)
    #define dR(i,j)     (dR + (i) + (j)*lddr)
    
    magma_minproductFloatComplex c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex c_one     = MAGMA_minproduct_C_ONE;
    magma_minproduct_int_t     ione  = 1;
    magma_minproductFloatComplex_ptr dR;
    magma_minproductFloatComplex Xnrmv, Rnrmv;
    float          Anrm, Xnrm, Rnrm, cte, eps;
    magma_minproduct_int_t     i, j, iiter, lddsa, lddr;
    
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
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }
    
    if ( n == 0 || nrhs == 0 )
        return *info;

    lddsa = n;
    lddr  = n;
    
    dR  = dworkd;
    
    eps  = lapackf77_slamch("Epsilon");
    Anrm = magma_minproductblas_clange(Magma_minproductInfNorm, n, n, dA, ldda, (float*)dworkd );
    cte  = Anrm * eps * pow( (float)n, (float)0.5 ) * BWDMAX;
    
    // residual dR = dB - dA*dX in real
    magma_minproductblas_clacpy( Magma_minproductUpperLower, n, nrhs, dB, lddb, dR, lddr );
    if ( nrhs == 1 ) {
        magma_minproduct_cgemv( trans, n, n,
                     c_neg_one, dA, ldda,
                                dX, 1,
                     c_one,     dR, 1 );
    }
    else {
        magma_minproduct_cgemm( trans, Magma_minproductNoTrans, n, nrhs, n,
                     c_neg_one, dA, ldda,
                                dX, lddx,
                     c_one,     dR, lddr );
    }
    
    // TODO: use MAGMA_minproduct_C_ABS( dX(i,j) ) instead of clange?
    for( j=0; j < nrhs; j++ ) {
        i = magma_minproduct_icamax( n, dX(0,j), 1) - 1;
        magma_minproduct_cgetmatrix( 1, 1, dX(i,j), 1, &Xnrmv, 1 );
        Xnrm = lapackf77_clange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
        
        i = magma_minproduct_icamax ( n, dR(0,j), 1 ) - 1;
        magma_minproduct_cgetmatrix( 1, 1, dR(i,j), 1, &Rnrmv, 1 );
        Rnrm = lapackf77_clange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
       


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
        magma_minproduct_cgetrs_nopiv_gpu( trans, n, nrhs, dAF, lddsa, dR, lddr, info );
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
            magma_minproductblas_caxpycp2( n, dR(0,j), dX(0,j), dB(0,j) );
        }
        
        // residual dR = dB - dA*dX in real
        if ( nrhs == 1 ) {
            magma_minproduct_cgemv( trans, n, n,
                         c_neg_one, dA, ldda,
                                    dX, 1,
                         c_one,     dR, 1 );
        }
        else {
            magma_minproduct_cgemm( trans, Magma_minproductNoTrans, n, nrhs, n,
                         c_neg_one, dA, ldda,
                                    dX, lddx,
                         c_one,     dR, lddr );
        }
        
        /*  Check whether the nrhs normwise backward errors satisfy the
         *  stopping criterion. If yes, set ITER=IITER > 0 and return. */
        for( j=0; j < nrhs; j++ ) {
            i = magma_minproduct_icamax( n, dX(0,j), 1) - 1;
            magma_minproduct_cgetmatrix( 1, 1, dX(i,j), 1, &Xnrmv, 1 );
            Xnrm = lapackf77_clange( "F", &ione, &ione, &Xnrmv, &ione, NULL );
            
            i = magma_minproduct_icamax ( n, dR(0,j), 1 ) - 1;
            magma_minproduct_cgetmatrix( 1, 1, dR(i,j), 1, &Rnrmv, 1 );
            Rnrm = lapackf77_clange( "F", &ione, &ione, &Rnrmv, &ione, NULL );
            
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
    ref  magma_minproduct_bool_t
    Specifies if iterative refinement have to be applied to improve the solution.
    -     = Magma_minproductTrue:   Iterative refinement is applied.
    -     = Magma_minproductFalse:  Iterative refinement is not applied.

    @param[in]
    n       INTEGER
    The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
    The number of right hand sides, i.e., the number of columns
    of the matrix B.  NRHS >= 0.

    @param[in,out]
    A       COMPLEX array, dimension (LDA,N).
    On entry, the M-by-N matrix to be factored.
    On exit, the factors L and U from the factorization
    A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    lda     INTEGER
    The leading dimension of the array A.  LDA >= max(1,N).


    @param[in,out]
    B       COMPLEX array, dimension (LDB,NRHS)
    On entry, the right hand side matrix B.
    On exit, the solution matrix X.

    @param[in]
    ldb     INTEGER
    The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
    -     = 0:  successful exit
    -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_cgesv_driver
 ********************************************************************/


extern "C" magma_minproduct_int_t 
magma_minproduct_cgesv_rbt(
    magma_minproduct_bool_t ref, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, 
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb, 
    magma_minproduct_int_t *info)
{

    /* Function Body */
    *info = 0;
    if ( ! (ref == Magma_minproductTrue) &&
         ! (ref == Magma_minproductFalse) ) {
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
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (nrhs == 0 || n == 0)
        return *info;


    magma_minproduct_int_t nn = n + ((4-(n % 4))%4);
    magma_minproductFloatComplex *dA, *hu, *hv, *db, *dAo, *dBo, *dwork;
    magma_minproduct_int_t n2;

    magma_minproduct_int_t iter;
    n2 = nn*nn;

    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dA, n2 )) {
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        return *info;
    }
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &db, nn*nrhs )) {
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        return *info;
    }

    if (ref == Magma_minproductTrue) {
        if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dAo, n2 )) {
            *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
            return *info;
        }
        if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dwork, nn*nrhs )) {
            *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
            return *info;
        }
        if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dBo, nn*nrhs )) {
            *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
            return *info;
        }
    }

    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc_cpu( &hu, 2*nn )) {
        *info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return *info;
    }
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc_cpu( &hv, 2*nn )) {
        *info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return *info;
    }

    magma_minproductblas_claset(Magma_minproductFull, nn, nn, MAGMA_minproduct_C_ZERO, MAGMA_minproduct_C_ONE, dA, nn);

    /* Send matrix on the GPU*/
    magma_minproduct_csetmatrix(n, n, A, lda, dA, nn);

    /* Send b on the GPU*/
    magma_minproduct_csetmatrix(n, nrhs, B, ldb, db, nn);

    *info = magma_minproduct_cgerbt_gpu(Magma_minproductTrue, nn, nrhs, dA, nn, db, nn, hu, hv, info);
    if (*info != MAGMA_minproduct_SUCCESS)  {
        return *info;
    }

    if (ref == Magma_minproductTrue) {
        magma_minproduct_ccopymatrix(nn, nn, dA, nn, dAo, nn);
        magma_minproduct_ccopymatrix(nn, nrhs, db, nn, dBo, nn);
    }
    /* Solve the system U^TAV.y = U^T.b on the GPU*/ 
    magma_minproduct_cgesv_nopiv_gpu( nn, nrhs, dA, nn, db, nn, info);


    /* Iterative refinement */
    if (ref == Magma_minproductTrue) {
        magma_minproduct_cgerfs_nopiv_gpu(Magma_minproductNoTrans, nn, nrhs, dAo, nn, dBo, nn, db, nn, dwork, dA, &iter, info);
    }
    //printf("iter = %d\n", iter);

    /* The solution of A.x = b is Vy computed on the GPU */
    magma_minproductFloatComplex *dv;

    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dv, 2*nn )) {
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        return *info;
    }

    magma_minproduct_csetvector(2*nn, hv, 1, dv, 1);
    
    for(int i = 0; i < nrhs; i++) {
        magma_minproductblas_cprbt_mv(nn, dv, db+(i*nn));
    }

    magma_minproduct_cgetmatrix(n, nrhs, db, nn, B, ldb);

    magma_minproduct_free_cpu( hu);
    magma_minproduct_free_cpu( hv);

    magma_minproduct_free( dA );
    magma_minproduct_free( dv );
    magma_minproduct_free( db );
    
    if (ref == Magma_minproductTrue) {    
        magma_minproduct_free( dAo );
        magma_minproduct_free( dBo );
        magma_minproduct_free( dwork );
    }
    return *info;
}
