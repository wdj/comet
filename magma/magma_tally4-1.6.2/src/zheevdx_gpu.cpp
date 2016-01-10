/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
    
       @author Raffaele Solca
       @author Stan Tomov
       @author Mark Gates
       @author Azzam Haidar
    
       @precisions normal z -> c

*/
#include "common_magma_tally4.h"
#include "magma_tally4_timer.h"

#define COMPLEX
#define FAST_HEMV

/**
    Purpose
    -------
    ZHEEVDX_GPU computes selected eigenvalues and, optionally, eigenvectors
    of a complex Hermitian matrix A. Eigenvalues and eigenvectors can
    be selected by specifying either a range of values or a range of
    indices for the desired eigenvalues.
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
    jobz    magma_tally4_vec_t
      -     = Magma_tally4NoVec:  Compute eigenvalues only;
      -     = Magma_tally4Vec:    Compute eigenvalues and eigenvectors.

    @param[in]
    range   magma_tally4_range_t
      -     = Magma_tally4RangeAll: all eigenvalues will be found.
      -     = Magma_tally4RangeV:   all eigenvalues in the half-open interval (VL,VU]
                   will be found.
      -     = Magma_tally4RangeI:   the IL-th through IU-th eigenvalues will be found.

    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of A is stored;
      -     = Magma_tally4Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU,
            dimension (LDDA, N).
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally4Upper, the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = Magma_tally4Lower,
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.
            On exit, if JOBZ = Magma_tally4Vec, then if INFO = 0, the first mout columns
            of A contains the required
            orthonormal eigenvectors of the matrix A.
            If JOBZ = Magma_tally4NoVec, then on exit the lower triangle (if UPLO=Magma_tally4Lower)
            or the upper triangle (if UPLO=Magma_tally4Upper) of A, including the
            diagonal, is destroyed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.  LDDA >= max(1,N).

    @param[in]
    vl      DOUBLE PRECISION
    @param[in]
    vu      DOUBLE PRECISION
            If RANGE=Magma_tally4RangeV, the lower and upper bounds of the interval to
            be searched for eigenvalues. VL < VU.
            Not referenced if RANGE = Magma_tally4RangeAll or Magma_tally4RangeI.

    @param[in]
    il      INTEGER
    @param[in]
    iu      INTEGER
            If RANGE=Magma_tally4RangeI, the indices (in ascending order) of the
            smallest and largest eigenvalues to be returned.
            1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
            Not referenced if RANGE = Magma_tally4RangeAll or Magma_tally4RangeV.

    @param[out]
    mout    INTEGER
            The total number of eigenvalues found.  0 <= MOUT <= N.
            If RANGE = Magma_tally4RangeAll, MOUT = N, and if RANGE = Magma_tally4RangeI, MOUT = IU-IL+1.

    @param[out]
    w       DOUBLE PRECISION array, dimension (N)
            If INFO = 0, the required mout eigenvalues in ascending order.

    @param
    wA      (workspace) COMPLEX_16 array, dimension (LDWA, N)

    @param[in]
    ldwa    INTEGER
            The leading dimension of the array wA.  LDWA >= max(1,N).

    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The length of the array WORK.
            If N <= 1,                      LWORK >= 1.
            If JOBZ = Magma_tally4NoVec and N > 1, LWORK >= N + N*NB.
            If JOBZ = Magma_tally4Vec   and N > 1, LWORK >= max( N + N*NB, 2*N + N**2 ).
            NB can be obtained through magma_tally4_get_zhetrd_nb(N).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal sizes of the WORK, RWORK and
            IWORK arrays, returns these values as the first entries of
            the WORK, RWORK and IWORK arrays, and no error message
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.

    @param[out]
    rwork   (workspace) DOUBLE PRECISION array, dimension (LRWORK)
            On exit, if INFO = 0, RWORK[0] returns the optimal LRWORK.

    @param[in]
    lrwork  INTEGER
            The dimension of the array RWORK.
            If N <= 1,                      LRWORK >= 1.
            If JOBZ = Magma_tally4NoVec and N > 1, LRWORK >= N.
            If JOBZ = Magma_tally4Vec   and N > 1, LRWORK >= 1 + 5*N + 2*N**2.
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
            If JOBZ = Magma_tally4NoVec and N > 1, LIWORK >= 1.
            If JOBZ = Magma_tally4Vec   and N > 1, LIWORK >= 3 + 5*N.
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
      -     > 0:  if INFO = i and JOBZ = Magma_tally4NoVec, then the algorithm failed
                  to converge; i off-diagonal elements of an intermediate
                  tridiagonal form did not converge to zero;
                  if INFO = i and JOBZ = Magma_tally4Vec, then the algorithm failed
                  to compute an eigenvalue while working on the submatrix
                  lying in rows and columns INFO/(N+1) through
                  mod(INFO,N+1).

    Further Details
    ---------------
    Based on contributions by
       Jeff Rutter, Computer Science Division, University of California
       at Berkeley, USA

    Modified description of INFO. Sven, 16 Feb 05.

    @ingroup magma_tally4_zheev_driver
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zheevdx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *mout, double *w,
    magma_tally4DoubleComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info)
{
    const char* uplo_  = lapack_uplo_const_tally4( uplo  );
    const char* jobz_  = lapack_vec_const_tally4( jobz  );
    magma_tally4_int_t ione = 1;

    double d__1;

    double eps;
    magma_tally4_int_t inde;
    double anrm;
    magma_tally4_int_t imax;
    double rmin, rmax;
    double sigma;
    magma_tally4_int_t iinfo, lwmin;
    magma_tally4_int_t lower;
    magma_tally4_int_t llrwk;
    magma_tally4_int_t wantz;
    //magma_tally4_int_t indwk2;
    magma_tally4_int_t iscale;
    double safmin;
    double bignum;
    magma_tally4_int_t indtau;
    magma_tally4_int_t indrwk, indwrk, liwmin;
    magma_tally4_int_t lrwmin, llwork;
    double smlnum;
    magma_tally4_int_t lquery;
    magma_tally4_int_t alleig, valeig, indeig;

    magma_tally4Double_ptr dwork;
    magma_tally4DoubleComplex_ptr dC;
    magma_tally4_int_t lddc = ldda;

    wantz = (jobz == Magma_tally4Vec);
    lower = (uplo == Magma_tally4Lower);

    alleig = (range == Magma_tally4RangeAll);
    valeig = (range == Magma_tally4RangeV);
    indeig = (range == Magma_tally4RangeI);

    lquery = (lwork == -1 || lrwork == -1 || liwork == -1);

    *info = 0;
    if (! (wantz || (jobz == Magma_tally4NoVec))) {
        *info = -1;
    } else if (! (alleig || valeig || indeig)) {
        *info = -2;
    } else if (! (lower || (uplo == Magma_tally4Upper))) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (ldda < max(1,n)) {
        *info = -6;
    } else if (ldwa < max(1,n)) {
        *info = -14;
    } else {
        if (valeig) {
            if (n > 0 && vu <= vl) {
                *info = -8;
            }
        } else if (indeig) {
            if (il < 1 || il > max(1,n)) {
                *info = -9;
            } else if (iu < min(n,il) || iu > n) {
                *info = -10;
            }
        }
    }

    magma_tally4_int_t nb = magma_tally4_get_zhetrd_nb( n );
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
    work[0]  = MAGMA_tally4_Z_MAKE( lwmin * one_eps, 0 );
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    if ((lwork < lwmin) && !lquery) {
        *info = -16;
    } else if ((lrwork < lrwmin) && ! lquery) {
        *info = -18;
    } else if ((liwork < liwmin) && ! lquery) {
        *info = -20;
    }

    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* If matrix is very small, then just call LAPACK on CPU, no need for GPU */
    if (n <= 128) {
        magma_tally4_int_t lda = n;
        magma_tally4DoubleComplex *A;
        magma_tally4_zmalloc_cpu( &A, lda*n );
        magma_tally4_zgetmatrix( n, n, dA, ldda, A, lda );
        lapackf77_zheevd( jobz_, uplo_,
                          &n, A, &lda,
                          w, work, &lwork,
                          rwork, &lrwork,
                          iwork, &liwork, info );
        magma_tally4_zsetmatrix( n, n, A, lda, dA, ldda );
        magma_tally4_free_cpu( A );
        *mout = n;
        return *info;
    }

    magma_tally4_queue_t stream;
    magma_tally4_queue_create( &stream );

    // dC and dwork are never used together, so use one buffer for both;
    // unfortunately they're different types (complex and double).
    // (this is easier in dsyevd_gpu where everything is double.)
    // zhetrd2_gpu requires ldda*ceildiv(n,64) + 2*ldda*nb, in double-complex.
    // zunmtr_gpu  requires lddc*n,                         in double-complex.
    // zlanhe      requires n, in double.
    magma_tally4_int_t ldwork = max( ldda*ceildiv(n,64) + 2*ldda*nb, lddc*n );
    magma_tally4_int_t ldwork_real = max( ldwork*2, n );
    if ( wantz ) {
        // zstedx requrise 3n^2/2, in double
        ldwork_real = max( ldwork_real, 3*n*(n/2 + 1) );
    }
    if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc( &dwork, ldwork_real )) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }
    dC = (magma_tally4DoubleComplex*) dwork;

    /* Get machine constants. */
    safmin = lapackf77_dlamch("Safe minimum");
    eps    = lapackf77_dlamch("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = magma_tally4_dsqrt( smlnum );
    rmax = magma_tally4_dsqrt( bignum );

    /* Scale matrix to allowable range, if necessary. */
    anrm = magma_tally4blas_zlanhe( Magma_tally4MaxNorm, uplo, n, dA, ldda, dwork );
    iscale = 0;
    sigma  = 1;
    if (anrm > 0. && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        magma_tally4blas_zlascl( uplo, 0, 0, 1., sigma, n, n, dA, ldda, info );
    }

    /* Call ZHETRD to reduce Hermitian matrix to tridiagonal form. */
    // zhetrd rwork: e (n)
    // zstedx rwork: e (n) + llrwk (1 + 4*N + 2*N**2)  ==>  1 + 5n + 2n^2
    inde   = 0;
    indrwk = inde + n;
    llrwk  = lrwork - indrwk;

    // zhetrd work: tau (n) + llwork (n*nb)  ==>  n + n*nb
    // zstedx work: tau (n) + z (n^2)
    // zunmtr work: tau (n) + z (n^2) + llwrk2 (n or n*nb)  ==>  2n + n^2, or n + n*nb + n^2
    indtau = 0;
    indwrk = indtau + n;
    //indwk2 = indwrk + n*n;
    llwork = lwork - indwrk;
    //llwrk2 = lwork - indwk2;

    magma_tally4_timer_t time=0;
    timer_start( time );

#ifdef FAST_HEMV
    magma_tally4_zhetrd2_gpu( uplo, n, dA, ldda, w, &rwork[inde],
                       &work[indtau], wA, ldwa, &work[indwrk], llwork,
                       dC, ldwork, &iinfo );
#else
    magma_tally4_zhetrd_gpu ( uplo, n, dA, ldda, w, &rwork[inde],
                       &work[indtau], wA, ldwa, &work[indwrk], llwork,
                       &iinfo );
#endif

    timer_stop( time );
    timer_printf( "time zhetrd_gpu = %6.2f\n", time );

    /* For eigenvalues only, call DSTERF.  For eigenvectors, first call
       ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
       tridiagonal matrix, then call ZUNMTR to multiply it to the Householder
       transformations represented as Householder vectors in A. */
    if (! wantz) {
        lapackf77_dsterf( &n, w, &rwork[inde], info );

        magma_tally4_dmove_eig( range, n, w, &il, &iu, vl, vu, mout );
    }
    else {
        timer_start( time );

        magma_tally4_zstedx( range, n, vl, vu, il, iu, w, &rwork[inde],
                      &work[indwrk], n, &rwork[indrwk],
                      llrwk, iwork, liwork, dwork, info );

        timer_stop( time );
        timer_printf( "time zstedx = %6.2f\n", time );
        timer_start( time );

        magma_tally4_dmove_eig( range, n, w, &il, &iu, vl, vu, mout );

        magma_tally4_zsetmatrix( n, *mout, &work[indwrk + n * (il-1) ], n, dC, lddc );

        magma_tally4_zunmtr_gpu( Magma_tally4Left, uplo, Magma_tally4NoTrans, n, *mout, dA, ldda, &work[indtau],
                          dC, lddc, wA, ldwa, &iinfo );

        magma_tally4_zcopymatrix( n, *mout, dC, lddc, dA, ldda );

        timer_stop( time );
        timer_printf( "time zunmtr_gpu + copy = %6.2f\n", time );
    }

    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1) {
        if (*info == 0) {
            imax = n;
        } else {
            imax = *info - 1;
        }
        d__1 = 1. / sigma;
        blasf77_dscal( &imax, &d__1, w, &ione );
    }

    work[0]  = MAGMA_tally4_Z_MAKE( lwmin * one_eps, 0 );  // round up
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    magma_tally4_queue_destroy( stream );
    magma_tally4_free( dwork );

    return *info;
} /* magma_tally4_zheevdx_gpu */