/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds

*/
#include "common_magma_tally2.h"

#define BWDMAX 1.0
#define ITERMAX 30

/**
    Purpose
    -------
    ZCGEQRSV solves the least squares problem
       min || A*X - B ||,
    where A is an M-by-N matrix and X and B are M-by-NRHS matrices.

    ZCGEQRSV first attempts to factorize the matrix in complex SINGLE PRECISION
    and use this factorization within an iterative refinement procedure
    to produce a solution with complex DOUBLE PRECISION norm-wise backward error
    quality (see below). If the approach fails the method switches to a
    complex DOUBLE PRECISION factorization and solve.

    The iterative refinement is not going to be a winning strategy if
    the ratio complex SINGLE PRECISION performance over complex DOUBLE PRECISION
    performance is too small. A reasonable strategy should take the
    number of right-hand sides and the size of the matrix into account.
    This might be done with a call to ILAENV in the future. Up to now, we
    always try iterative refinement.
    
    The iterative refinement process is stopped if
        ITER > ITERMAX
    or for all the RHS we have:
        RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
    where
        o ITER is the number of the current iteration in the iterative
          refinement process
        o RNRM is the infinity-norm of the residual
        o XNRM is the infinity-norm of the solution
        o ANRM is the infinity-operator-norm of the matrix A
        o EPS is the machine epsilon returned by DLAMCH('Epsilon')
    The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00 respectively.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A. M >= N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the M-by-N coefficient matrix A.
            On exit, if iterative refinement has been successfully used
            (info.EQ.0 and ITER.GE.0, see description below), A is
            unchanged. If double precision factorization has been used
            (info.EQ.0 and ITER.LT.0, see description below), then the
            array dA contains the QR factorization of A as returned by
            function DGEQRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).

    @param[in,out]
    dB      COMPLEX_16 array on the GPU, dimension (LDDB,NRHS)
            The M-by-NRHS right hand side matrix B.
            May be overwritten (e.g., if refinement fails).

    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  LDDB >= max(1,M).

    @param[out]
    dX      COMPLEX_16 array on the GPU, dimension (LDDX,NRHS)
            If info = 0, the N-by-NRHS solution matrix X.

    @param[in]
    lddx    INTEGER
            The leading dimension of the array dX.  LDDX >= max(1,N).

    @param[out]
    iter    INTEGER
      -     < 0: iterative refinement has failed, double precision
                 factorization has been performed
        +        -1 : the routine fell back to full precision for
                      implementation- or machine-specific reasons
        +        -2 : narrowing the precision induced an overflow,
                      the routine fell back to full precision
        +        -3 : failure of SGEQRF
        +        -31: stop the iterative refinement after the 30th iteration
      -     > 0: iterative refinement has been successfully used.
                 Returns the number of iterations

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if info = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_zgels_driver
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zcgeqrsv_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB,  magma_tally2_int_t lddb,
    magma_tally2DoubleComplex_ptr dX,  magma_tally2_int_t lddx,
    magma_tally2_int_t *iter,
    magma_tally2_int_t *info)
{
    #define dB(i,j)     (dB + (i) + (j)*lddb)
    #define dX(i,j)     (dX + (i) + (j)*lddx)
    #define dR(i,j)     (dR + (i) + (j)*lddr)
    #define dSX(i,j)    (dSX + (i) + (j)*lddsx)
    
    magma_tally2DoubleComplex c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    magma_tally2DoubleComplex c_one     = MAGMA_tally2_Z_ONE;
    magma_tally2_int_t     ione  = 1;
    magma_tally2DoubleComplex *hworkd;
    magma_tally2FloatComplex  *hworks;
    magma_tally2DoubleComplex *tau;
    magma_tally2FloatComplex  *stau;
    magma_tally2DoubleComplex_ptr dworkd;
    magma_tally2FloatComplex_ptr  dworks;
    magma_tally2DoubleComplex_ptr dR, dT;
    magma_tally2FloatComplex_ptr  dSA, dSX, dST;
    magma_tally2DoubleComplex Xnrmv, Rnrmv;
    double          Anrm, Xnrm, Rnrm, cte, eps;
    magma_tally2_int_t     i, j, iiter, lddsa, lddsx, lddr, nb, lhwork, minmn, size, ldworkd;

    /* Check arguments */
    *iter = 0;
    *info = 0;
    if ( m < 0 )
        *info = -1;
    else if ( n < 0 || n > m )
        *info = -2;
    else if ( nrhs < 0 )
        *info = -3;
    else if ( ldda < max(1,m))
        *info = -5;
    else if ( lddb < max(1,m))
        *info = -7;
    else if ( lddx < max(1,n))
        *info = -9;

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    if ( m == 0 || n == 0 || nrhs == 0 )
        return *info;

    nb   = magma_tally2_get_cgeqrf_nb(m);
    minmn= min(m, n);
    
    /* dSX contains both B and X, so must be max(m or lddb,n). */
    lddsa = ldda;
    lddsx = max(lddb,n);
    lddr  = lddb;
    
    /*
     * Allocate temporary buffers
     */
    /* dworks(dSA + dSX + dST) */
    size = lddsa*n + lddsx*nrhs + ( 2*minmn + ((n+31)/32)*32 )*nb;
    if (MAGMA_tally2_SUCCESS != magma_tally2_cmalloc( &dworks, size )) {
        fprintf(stderr, "Allocation of dworks failed (%d)\n", (int) size);
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }
    dSA = dworks;
    dSX = dSA + lddsa*n;
    dST = dSX + lddsx*nrhs;

    /* dworkd(dR) = lddr*nrhs */
    ldworkd = lddr*nrhs;
    if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dworkd, ldworkd )) {
        magma_tally2_free( dworks );
        fprintf(stderr, "Allocation of dworkd failed\n");
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }
    dR = dworkd;

    /* hworks(workspace for cgeqrs + stau) = min(m,n) + lhworks */
    lhwork = (m - n + nb)*(nrhs + nb) + nrhs*nb;
    size = lhwork + minmn;
    magma_tally2_cmalloc_cpu( &hworks, size );
    if ( hworks == NULL ) {
        magma_tally2_free( dworks );
        magma_tally2_free( dworkd );
        fprintf(stderr, "Allocation of hworks failed\n");
        *info = MAGMA_tally2_ERR_HOST_ALLOC;
        return *info;
    }
    stau = hworks + lhwork;

    eps  = lapackf77_dlamch("Epsilon");
    Anrm = magma_tally2blas_zlange(Magma_tally2InfNorm, m, n, dA, ldda, (double*)dworkd );
    cte  = Anrm * eps * pow((double)n, 0.5) * BWDMAX;

    /*
     * Convert to single precision
     */
    magma_tally2blas_zlag2c( m, nrhs, dB, lddb, dSX, lddsx, info );
    if (*info != 0) {
        *iter = -2;
        goto FALLBACK;
    }

    magma_tally2blas_zlag2c( m, n, dA, ldda, dSA, lddsa, info );
    if (*info != 0) {
        *iter = -2;
        goto FALLBACK;
    }

    // factor dSA in single precision
    magma_tally2_cgeqrf_gpu( m, n, dSA, lddsa, stau, dST, info );
    if (*info != 0) {
        *iter = -3;
        goto FALLBACK;
    }

    // solve dSA*dSX = dB in single precision
    magma_tally2_cgeqrs_gpu( m, n, nrhs, dSA, lddsa, stau, dST, dSX, lddsx, hworks, lhwork, info );
    if (*info != 0) {
        *iter = -3;
        goto FALLBACK;
    }

    // residual dR = dB - dA*dX in double precision
    magma_tally2blas_clag2z( n, nrhs, dSX, lddsx, dX, lddx, info );
    magma_tally2blas_zlacpy( Magma_tally2UpperLower, m, nrhs, dB, lddb, dR, lddr );
    if ( nrhs == 1 ) {
        magma_tally2_zgemv( Magma_tally2NoTrans, m, n,
                     c_neg_one, dA, ldda,
                                dX, 1,
                     c_one,     dR, 1 );
    }
    else {
        magma_tally2_zgemm( Magma_tally2NoTrans, Magma_tally2NoTrans, m, nrhs, n,
                     c_neg_one, dA, ldda,
                                dX, lddx,
                     c_one,     dR, lddr );
    }

    // TODO: use MAGMA_tally2_Z_ABS( dX(i,j) ) instead of zlange?
    for( j=0; j < nrhs; j++ ) {
        i = magma_tally2_izamax( n, dX(0,j), 1) - 1;
        magma_tally2_zgetmatrix( 1, 1, dX(i,j), 1, &Xnrmv, 1 );
        Xnrm = lapackf77_zlange( "F", &ione, &ione, &Xnrmv, &ione, NULL );

        i = magma_tally2_izamax ( m, dR(0,j), 1 ) - 1;
        magma_tally2_zgetmatrix( 1, 1, dR(i,j), 1, &Rnrmv, 1 );
        Rnrm = lapackf77_zlange( "F", &ione, &ione, &Rnrmv, &ione, NULL );

        if ( Rnrm >  Xnrm*cte ) {
            goto REFINEMENT;
        }
    }

    *iter = 0;

    /* Free workspaces */
    magma_tally2_free( dworks );
    magma_tally2_free( dworkd );
    magma_tally2_free_cpu( hworks );
    return *info;

REFINEMENT:
    /* TODO: this iterative refinement algorithm works only for compatibile
     * systems (B in colspan of A).
     * See Matrix Computations (3rd ed) p. 267 for correct algorithm. */
    for( iiter=1; iiter < ITERMAX; ) {
        *info = 0;
        // convert residual dR to single precision dSX
        magma_tally2blas_zlag2c( m, nrhs, dR, lddr, dSX, lddsx, info );
        if (*info != 0) {
            *iter = -2;
            goto FALLBACK;
        }
        // solve dSA*dSX = R in single precision
        magma_tally2_cgeqrs_gpu( m, n, nrhs, dSA, lddsa, stau, dST, dSX, lddsx, hworks, lhwork, info );
        if (*info != 0) {
            *iter = -3;
            goto FALLBACK;
        }

        // Add correction and setup residual
        // dX += dSX [including conversion]  --and--
        // dR[1:n] = dB[1:n]   (only n rows, not whole m rows! -- useless if m > n)
        for( j=0; j < nrhs; j++ ) {
            magma_tally2blas_zcaxpycp( n, dSX(0,j), dX(0,j), dB(0,j), dR(0,j) );
        }
        // dR = dB  (whole m rows)
        magma_tally2blas_zlacpy( Magma_tally2UpperLower, m, nrhs, dB, lddb, dR, lddr );
        
        // residual dR = dB - dA*dX in double precision
        if ( nrhs == 1 ) {
            magma_tally2_zgemv( Magma_tally2NoTrans, m, n,
                         c_neg_one, dA, ldda,
                                    dX, 1,
                         c_one,     dR, 1 );
        }
        else {
            magma_tally2_zgemm( Magma_tally2NoTrans, Magma_tally2NoTrans, m, nrhs, n,
                         c_neg_one, dA, ldda,
                                    dX, lddx,
                         c_one,     dR, lddr );
        }

        /*  Check whether the nrhs normwise backward errors satisfy the
         *  stopping criterion. If yes, set ITER=IITER > 0 and return. */
        for( j=0; j < nrhs; j++ ) {
            i = magma_tally2_izamax( n, dX(0,j), 1) - 1;
            magma_tally2_zgetmatrix( 1, 1, dX(i,j), 1, &Xnrmv, 1 );
            Xnrm = lapackf77_zlange( "F", &ione, &ione, &Xnrmv, &ione, NULL );

            i = magma_tally2_izamax ( m, dR(0,j), 1 ) - 1;
            magma_tally2_zgetmatrix( 1, 1, dR(i,j), 1, &Rnrmv, 1 );
            Rnrm = lapackf77_zlange( "F", &ione, &ione, &Rnrmv, &ione, NULL );

            if ( Rnrm >  Xnrm*cte ) {
                goto L20;
            }
        }

        /*  If we are here, the nrhs normwise backward errors satisfy
         *  the stopping criterion, we are good to exit. */
        *iter = iiter;

        /* Free workspaces */
        magma_tally2_free( dworks );
        magma_tally2_free( dworkd );
        magma_tally2_free_cpu( hworks );
        return *info;
        
      L20:
        iiter++;
    }

    /* If we are at this place of the code, this is because we have
     * performed ITER=ITERMAX iterations and never satisified the
     * stopping criterion. Set up the ITER flag accordingly and follow
     * up on double precision routine. */
    *iter = -ITERMAX - 1;
    
FALLBACK:
    /* Single-precision iterative refinement failed to converge to a
     * satisfactory solution, so we resort to double precision. */
    magma_tally2_free( dworks );
    magma_tally2_free_cpu( hworks );

    /*
     * Allocate temporary buffers
     */
    /* dworkd = dT for zgeqrf */
    nb   = magma_tally2_get_zgeqrf_nb( m );
    size = (2*min(m, n) + (n+31)/32*32 )*nb;
    if ( size > ldworkd ) {
        magma_tally2_free( dworkd );
        if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dworkd, size )) {
            fprintf(stderr, "Allocation of dworkd2 failed\n");
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }
    }
    dT = dworkd;

    /* hworkd(dtau + workspace for zgeqrs) = min(m,n) + lhwork */
    size = lhwork + minmn;
    magma_tally2_zmalloc_cpu( &hworkd, size );
    if ( hworkd == NULL ) {
        magma_tally2_free( dworkd );
        fprintf(stderr, "Allocation of hworkd2 failed\n");
        *info = MAGMA_tally2_ERR_HOST_ALLOC;
        return *info;
    }
    tau = hworkd + lhwork;

    magma_tally2_zgeqrf_gpu( m, n, dA, ldda, tau, dT, info );
    if (*info == 0) {
        // if m > n, then dB won't fit in dX, so solve with dB and copy n rows to dX
        magma_tally2_zgeqrs_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hworkd, lhwork, info );
        magma_tally2blas_zlacpy( Magma_tally2UpperLower, n, nrhs, dB, lddb, dX, lddx );
    }

    magma_tally2_free( dworkd );
    magma_tally2_free_cpu( hworkd );
    return *info;
}
