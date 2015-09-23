/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar
       @author Mark Gates

       @precisions normal z -> c

*/
#include "common_magma_minproduct.h"
#include "magma_minproduct_timer.h"

#define COMPLEX

/**
    Purpose
    -------
    ZHEGVDX computes selected eigenvalues and, optionally, eigenvectors
    of a complex generalized Hermitian-definite eigenproblem, of the form
    A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
    B are assumed to be Hermitian and B is also positive definite.
    Eigenvalues and eigenvectors can be selected by specifying either a
    range of values or a range of indices for the desired eigenvalues.
    If eigenvectors are desired, it uses a divide and conquer algorithm.

    The divide and conquer algorithm makes very mild assumptions about
    floating point arithmetic. It will work on machines with a guard
    digit in add/subtract, or on those binary machines without guard
    digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    Cray-2. It could conceivably fail on hexadecimal or decimal machines
    without guard digits, but we know of none.

    Arguments
    ---------
    @param[in]
    itype   INTEGER
            Specifies the problem type to be solved:
            = 1:  A*x = (lambda)*B*x
            = 2:  A*B*x = (lambda)*x
            = 3:  B*A*x = (lambda)*x

    @param[in]
    jobz    magma_minproduct_vec_t
      -     = Magma_minproductNoVec:  Compute eigenvalues only;
      -     = Magma_minproductVec:    Compute eigenvalues and eigenvectors.

    @param[in]
    range   magma_minproduct_range_t
      -     = Magma_minproductRangeAll: all eigenvalues will be found.
      -     = Magma_minproductRangeV:   all eigenvalues in the half-open interval (VL,VU]
                   will be found.
      -     = Magma_minproductRangeI:   the IL-th through IU-th eigenvalues will be found.

    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  Upper triangles of A and B are stored;
      -     = Magma_minproductLower:  Lower triangles of A and B are stored.

    @param[in]
    n       INTEGER
            The order of the matrices A and B.  N >= 0.

    @param[in,out]
    A       COMPLEX_16 array, dimension (LDA, N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_minproductUpper, the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = Magma_minproductLower,
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.
    \n
            On exit, if JOBZ = Magma_minproductVec, then if INFO = 0, A contains the
            matrix Z of eigenvectors.  The eigenvectors are normalized
            as follows:
            if ITYPE = 1 or 2, Z**H*B*Z = I;
            if ITYPE = 3, Z**H*inv(B)*Z = I.
            If JOBZ = Magma_minproductNoVec, then on exit the upper triangle (if UPLO=Magma_minproductUpper)
            or the lower triangle (if UPLO=Magma_minproductLower) of A, including the
            diagonal, is destroyed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in,out]
    B       COMPLEX_16 array, dimension (LDB, N)
            On entry, the Hermitian matrix B.  If UPLO = Magma_minproductUpper, the
            leading N-by-N upper triangular part of B contains the
            upper triangular part of the matrix B.  If UPLO = Magma_minproductLower,
            the leading N-by-N lower triangular part of B contains
            the lower triangular part of the matrix B.
    \n
            On exit, if INFO <= N, the part of B containing the matrix is
            overwritten by the triangular factor U or L from the Cholesky
            factorization B = U**H*U or B = L*L**H.

    @param[in]
    ldb     INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[in]
    vl      DOUBLE PRECISION
    @param[in]
    vu      DOUBLE PRECISION
            If RANGE=Magma_minproductRangeV, the lower and upper bounds of the interval to
            be searched for eigenvalues. VL < VU.
            Not referenced if RANGE = Magma_minproductRangeAll or Magma_minproductRangeI.

    @param[in]
    il      INTEGER
    @param[in]
    iu      INTEGER
            If RANGE=Magma_minproductRangeI, the indices (in ascending order) of the
            smallest and largest eigenvalues to be returned.
            1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
            Not referenced if RANGE = Magma_minproductRangeAll or Magma_minproductRangeV.

    @param[out]
    mout    INTEGER
            The total number of eigenvalues found.  0 <= MOUT <= N.
            If RANGE = Magma_minproductRangeAll, MOUT = N, and if RANGE = Magma_minproductRangeI, MOUT = IU-IL+1.

    @param[out]
    w       DOUBLE PRECISION array, dimension (N)
            If INFO = 0, the eigenvalues in ascending order.

    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The length of the array WORK.
            If N <= 1,                      LWORK >= 1.
            If JOBZ = Magma_minproductNoVec and N > 1, LWORK >= N + N*NB.
            If JOBZ = Magma_minproductVec   and N > 1, LWORK >= max( N + N*NB, 2*N + N**2 ).
            NB can be obtained through magma_minproduct_get_zhetrd_nb(N).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal sizes of the WORK, RWORK and
            IWORK arrays, returns these values as the first entries of
            the WORK, RWORK and IWORK arrays, and no error message
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.

    @param[out]
    rwork   (workspace) DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
            On exit, if INFO = 0, RWORK[0] returns the optimal LRWORK.

    @param[in]
    lrwork  INTEGER
            The dimension of the array RWORK.
            If N <= 1,                      LRWORK >= 1.
            If JOBZ = Magma_minproductNoVec and N > 1, LRWORK >= N.
            If JOBZ = Magma_minproductVec   and N > 1, LRWORK >= 1 + 5*N + 2*N**2.
    \n
            If LRWORK = -1, then a workspace query is assumed; the
            routine only calculates the optimal sizes of the WORK, RWORK
            and IWORK arrays, returns these values as the first entries
            of the WORK, RWORK and IWORK arrays, and no error message
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.

    @param[out]
    iwork   (workspace) INTEGER array, dimension (MAX(1,LIWORK))
            On exit, if INFO = 0, IWORK[0] returns the optimal LIWORK.

    @param[in]
    liwork  INTEGER
            The dimension of the array IWORK.
            If N <= 1,                      LIWORK >= 1.
            If JOBZ = Magma_minproductNoVec and N > 1, LIWORK >= 1.
            If JOBZ = Magma_minproductVec   and N > 1, LIWORK >= 3 + 5*N.
    \n
            If LIWORK = -1, then a workspace query is assumed; the
            routine only calculates the optimal sizes of the WORK, RWORK
            and IWORK arrays, returns these values as the first entries
            of the WORK, RWORK and IWORK arrays, and no error message
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  ZPOTRF or ZHEEVD returned an error code:
               <= N:  if INFO = i and JOBZ = Magma_minproductNoVec, then the algorithm
                      failed to converge; i off-diagonal elements of an
                      intermediate tridiagonal form did not converge to
                      zero;
                      if INFO = i and JOBZ = Magma_minproductVec, then the algorithm
                      failed to compute an eigenvalue while working on
                      the submatrix lying in rows and columns INFO/(N+1)
                      through mod(INFO,N+1);
               > N:   if INFO = N + i, for 1 <= i <= N, then the leading
                      minor of order i of B is not positive definite.
                      The factorization of B could not be completed and
                      no eigenvalues or eigenvectors were computed.

    Further Details
    ---------------
    Based on contributions by
       Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA

    Modified so that no backsubstitution is performed if ZHEEVD fails to
    converge (NEIG in old code could be greater than N causing out of
    bounds reference to A - reported by Ralf Meyer).  Also corrected the
    description of INFO and the test on ITYPE. Sven, 16 Feb 05.

    @ingroup magma_minproduct_zhegv_driver
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_zhegvdx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *mout, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info)
{
    const char* uplo_  = lapack_uplo_const( uplo  );
    const char* jobz_  = lapack_vec_const( jobz  );

    magma_minproductDoubleComplex c_one = MAGMA_minproduct_Z_ONE;

    magma_minproductDoubleComplex *dA=NULL, *dB=NULL;
    magma_minproduct_int_t ldda = roundup( n, 32 );
    magma_minproduct_int_t lddb = ldda;

    magma_minproduct_int_t lower;
    magma_minproduct_trans_t trans;
    magma_minproduct_int_t wantz;
    magma_minproduct_int_t lquery;
    magma_minproduct_int_t alleig, valeig, indeig;

    magma_minproduct_int_t lwmin;
    magma_minproduct_int_t liwmin;
    magma_minproduct_int_t lrwmin;

    magma_minproduct_queue_t stream;
    magma_minproduct_queue_create( &stream );

    wantz  = (jobz == Magma_minproductVec);
    lower  = (uplo == Magma_minproductLower);
    alleig = (range == Magma_minproductRangeAll);
    valeig = (range == Magma_minproductRangeV);
    indeig = (range == Magma_minproductRangeI);
    lquery = (lwork == -1 || lrwork == -1 || liwork == -1);

    *info = 0;
    if (itype < 1 || itype > 3) {
        *info = -1;
    } else if (! (alleig || valeig || indeig)) {
        *info = -2;
    } else if (! (wantz || (jobz == Magma_minproductNoVec))) {
        *info = -3;
    } else if (! (lower || (uplo == Magma_minproductUpper))) {
        *info = -4;
    } else if (n < 0) {
        *info = -5;
    } else if (lda < max(1,n)) {
        *info = -7;
    } else if (ldb < max(1,n)) {
        *info = -9;
    } else {
        if (valeig) {
            if (n > 0 && vu <= vl) {
                *info = -11;
            }
        } else if (indeig) {
            if (il < 1 || il > max(1,n)) {
                *info = -12;
            } else if (iu < min(n,il) || iu > n) {
                *info = -13;
            }
        }
    }

    magma_minproduct_int_t nb = magma_minproduct_get_zhetrd_nb( n );
    if ( n <= 1 ) {
        lwmin  = 1;
        lrwmin = 1;
        liwmin = 1;
    }
    else if ( wantz ) {
        lwmin  = max( n + n*nb, 2*n + n*n );
        lrwmin = 1 + 5*n + 2*n*n;
        liwmin = 3 + 5*n;
    }
    else {
        lwmin  = n + n*nb;
        lrwmin = n;
        liwmin = 1;
    }

    // multiply by 1+eps (in Double!) to ensure length gets rounded up,
    // if it cannot be exactly represented in floating point.
    real_Double_t one_eps = 1. + lapackf77_dlamch("Epsilon");
    work[0]  = MAGMA_minproduct_Z_MAKE( lwmin * one_eps, 0 );  // round up
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    if (lwork < lwmin && ! lquery) {
        *info = -17;
    } else if (lrwork < lrwmin && ! lquery) {
        *info = -19;
    } else if (liwork < liwmin && ! lquery) {
        *info = -21;
    }

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (n == 0) {
        return *info;
    }

    /* If matrix is very small, then just call LAPACK on CPU, no need for GPU */
    if (n <= 128) {
        lapackf77_zhegvd( &itype, jobz_, uplo_,
                          &n, A, &lda, B, &ldb,
                          w, work, &lwork,
                          #ifdef COMPLEX
                          rwork, &lrwork,
                          #endif
                          iwork, &liwork, info );
        *mout = n;
        return *info;
    }

    if (MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc( &dA, n*ldda ) ||
        MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc( &dB, n*lddb )) {
        magma_minproduct_free( dA );
        magma_minproduct_free( dB );
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        return *info;
    }

    /* Form a Cholesky factorization of B. */
    magma_minproduct_zsetmatrix( n, n, B, ldb, dB, lddb );

    magma_minproduct_zsetmatrix_async( n, n,
                            A,  lda,
                            dA, ldda, stream );

    magma_minproduct_timer_t time=0;
    timer_start( time );
    magma_minproduct_zpotrf_gpu( uplo, n, dB, lddb, info );
    if (*info != 0) {
        *info = n + *info;
        return *info;
    }
    timer_stop( time );
    timer_printf( "time zpotrf_gpu = %6.2f\n", time );

    magma_minproduct_queue_sync( stream );
    magma_minproduct_zgetmatrix_async( n, n,
                            dB, lddb,
                            B,  ldb, stream );

    /* Transform problem to standard eigenvalue problem and solve. */
    timer_start( time );
    magma_minproduct_zhegst_gpu( itype, uplo, n, dA, ldda, dB, lddb, info );
    timer_stop( time );
    timer_printf( "time zhegst_gpu = %6.2f\n", time );

    /* simple fix to be able to run bigger size.
     * set dB=NULL so we know to re-allocate below
     * TODO: have dwork here that will be used as dB and then passed to  dsyevd.
     */
    if (n > 5000) {
        magma_minproduct_queue_sync( stream );
        magma_minproduct_free( dB );  dB=NULL;
    }

    timer_start( time );
    magma_minproduct_zheevdx_gpu( jobz, range, uplo, n, dA, ldda, vl, vu, il, iu, mout, w, A, lda,
                       work, lwork, rwork, lrwork, iwork, liwork, info );
    timer_stop( time );
    timer_printf( "time zheevdx_gpu = %6.2f\n", time );

    if (wantz && *info == 0) {
        timer_start( time );
        
        /* allocate and copy dB back */
        if (dB == NULL) {
            if (MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc( &dB, n*lddb ) ) {
                magma_minproduct_free( dA );  dA=NULL;
                *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
                return *info;
            }
            magma_minproduct_zsetmatrix( n, n, B, ldb, dB, lddb );
        }
        /* Backtransform eigenvectors to the original problem. */
        if (itype == 1 || itype == 2) {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
               backtransform eigenvectors: x = inv(L)'*y or inv(U)*y */
            if (lower) {
                trans = Magma_minproductConjTrans;
            } else {
                trans = Magma_minproductNoTrans;
            }
            magma_minproduct_ztrsm( Magma_minproductLeft, uplo, trans, Magma_minproductNonUnit,
                         n, *mout, c_one, dB, lddb, dA, ldda );
        }
        else if (itype == 3) {
            /* For B*A*x=(lambda)*x;
               backtransform eigenvectors: x = L*y or U'*y */
            if (lower) {
                trans = Magma_minproductNoTrans;
            } else {
                trans = Magma_minproductConjTrans;
            }
            magma_minproduct_ztrmm( Magma_minproductLeft, uplo, trans, Magma_minproductNonUnit,
                         n, *mout, c_one, dB, lddb, dA, ldda );
        }

        magma_minproduct_zgetmatrix( n, *mout, dA, ldda, A, lda );
        
        timer_stop( time );
        timer_printf( "time ztrsm/mm + getmatrix = %6.2f\n", time );
    }

    magma_minproduct_queue_sync( stream );
    magma_minproduct_queue_destroy( stream );

    work[0]  = MAGMA_minproduct_Z_MAKE( lwmin * one_eps, 0 );  // round up
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    magma_minproduct_free( dA );  dA=NULL;
    magma_minproduct_free( dB );  dB=NULL;

    return *info;
} /* magma_minproduct_zhegvdx */
