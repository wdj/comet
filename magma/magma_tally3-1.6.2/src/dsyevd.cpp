/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Mark Gates

       @precisions normal d -> s

*/
#include "common_magma_tally3.h"
#include "magma_tally3_timer.h"

#define PRECISION_d
#define REAL

/**
    Purpose
    -------
    DSYEVD computes all eigenvalues and, optionally, eigenvectors of
    a real symmetric matrix A.  If eigenvectors are desired, it uses a
    divide and conquer algorithm.

    The divide and conquer algorithm makes very mild assumptions about
    floating point arithmetic. It will work on machines with a guard
    digit in add/subtract, or on those binary machines without guard
    digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    Cray-2. It could conceivably fail on hexadecimal or decimal machines
    without guard digits, but we know of none.

    Arguments
    ---------
    @param[in]
    jobz    magma_tally3_vec_t
      -     = Magma_tally3NoVec:  Compute eigenvalues only;
      -     = Magma_tally3Vec:    Compute eigenvalues and eigenvectors.

    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  Upper triangle of A is stored;
      -     = Magma_tally3Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       DOUBLE_PRECISION array, dimension (LDA, N)
            On entry, the symmetric matrix A.  If UPLO = Magma_tally3Upper, the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = Magma_tally3Lower,
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.
            On exit, if JOBZ = Magma_tally3Vec, then if INFO = 0, A contains the
            orthonormal eigenvectors of the matrix A.
            If JOBZ = Magma_tally3NoVec, then on exit the lower triangle (if UPLO=Magma_tally3Lower)
            or the upper triangle (if UPLO=Magma_tally3Upper) of A, including the
            diagonal, is destroyed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    w       DOUBLE PRECISION array, dimension (N)
            If INFO = 0, the eigenvalues in ascending order.

    @param[out]
    work    (workspace) DOUBLE_PRECISION array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The length of the array WORK.
            If N <= 1,                      LWORK >= 1.
            If JOBZ = Magma_tally3NoVec and N > 1, LWORK >= 2*N + N*NB.
            If JOBZ = Magma_tally3Vec   and N > 1, LWORK >= max( 2*N + N*NB, 1 + 6*N + 2*N**2 ).
            NB can be obtained through magma_tally3_get_dsytrd_nb(N).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal sizes of the WORK and IWORK
            arrays, returns these values as the first entries of the WORK
            and IWORK arrays, and no error message related to LWORK or
            LIWORK is issued by XERBLA.

    @param[out]
    iwork   (workspace) INTEGER array, dimension (MAX(1,LIWORK))
            On exit, if INFO = 0, IWORK[0] returns the optimal LIWORK.

    @param[in]
    liwork  INTEGER
            The dimension of the array IWORK.
            If N <= 1,                      LIWORK >= 1.
            If JOBZ = Magma_tally3NoVec and N > 1, LIWORK >= 1.
            If JOBZ = Magma_tally3Vec   and N > 1, LIWORK >= 3 + 5*N.
    \n
            If LIWORK = -1, then a workspace query is assumed; the
            routine only calculates the optimal sizes of the WORK and
            IWORK arrays, returns these values as the first entries of
            the WORK and IWORK arrays, and no error message related to
            LWORK or LIWORK is issued by XERBLA.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i and JOBZ = Magma_tally3NoVec, then the algorithm failed
                  to converge; i off-diagonal elements of an intermediate
                  tridiagonal form did not converge to zero;
                  if INFO = i and JOBZ = Magma_tally3Vec, then the algorithm failed
                  to compute an eigenvalue while working on the submatrix
                  lying in rows and columns INFO/(N+1) through
                  mod(INFO,N+1).

    Further Details
    ---------------
    Based on contributions by
       Jeff Rutter, Computer Science Division, University of California
       at Berkeley, USA

    Modified description of INFO. Sven, 16 Feb 05.

    @ingroup magma_tally3_dsyev_driver
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_dsyevd(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info)
{
    const char* uplo_ = lapack_uplo_const_tally3( uplo );
    const char* jobz_ = lapack_vec_const_tally3( jobz );
    magma_tally3_int_t ione = 1;
    magma_tally3_int_t izero = 0;
    double d_one = 1.;

    double d__1;

    double eps;
    magma_tally3_int_t inde;
    double anrm;
    double rmin, rmax;
    double sigma;
    magma_tally3_int_t iinfo, lwmin;
    magma_tally3_int_t lower;
    magma_tally3_int_t wantz;
    magma_tally3_int_t indwk2, llwrk2;
    magma_tally3_int_t iscale;
    double safmin;
    double bignum;
    magma_tally3_int_t indtau;
    magma_tally3_int_t indwrk, liwmin;
    magma_tally3_int_t llwork;
    double smlnum;
    magma_tally3_int_t lquery;

    double* dwork;

    wantz = (jobz == Magma_tally3Vec);
    lower = (uplo == Magma_tally3Lower);
    lquery = (lwork == -1 || liwork == -1);

    *info = 0;

    if (! (wantz || (jobz == Magma_tally3NoVec))) {
        *info = -1;
    } else if (! (lower || (uplo == Magma_tally3Upper))) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (lda < max(1,n)) {
        *info = -5;
    }

    magma_tally3_int_t nb = magma_tally3_get_dsytrd_nb( n );
    if ( n <= 1 ) {
        lwmin  = 1;
        liwmin = 1;
    }
    else if ( wantz ) {
        lwmin  = max( 2*n + n*nb, 1 + 6*n + 2*n*n );
        liwmin = 3 + 5*n;
    }
    else {
        lwmin  = 2*n + n*nb;
        liwmin = 1;
    }
    
    // multiply by 1+eps (in Double!) to ensure length gets rounded up,
    // if it cannot be exactly represented in floating point.
    real_Double_t one_eps = 1. + lapackf77_dlamch("Epsilon");
    work[0]  = lwmin * one_eps;
    iwork[0] = liwmin;

    if ((lwork < lwmin) && !lquery) {
        *info = -8;
    } else if ((liwork < liwmin) && ! lquery) {
        *info = -10;
    }

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (n == 0) {
        return *info;
    }

    if (n == 1) {
        w[0] = A[0];
        if (wantz) {
            A[0] = 1.;
        }
        return *info;
    }
    
    /* If matrix is very small, then just call LAPACK on CPU, no need for GPU */
    if (n <= 128) {
        lapackf77_dsyevd( jobz_, uplo_,
                          &n, A, &lda,
                          w, work, &lwork,
                          iwork, &liwork, info );
        return *info;
    }

    /* Get machine constants. */
    safmin = lapackf77_dlamch("Safe minimum");
    eps    = lapackf77_dlamch("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = magma_tally3_dsqrt( smlnum );
    rmax = magma_tally3_dsqrt( bignum );

    /* Scale matrix to allowable range, if necessary. */
    anrm = lapackf77_dlansy( "M", uplo_, &n, A, &lda, work );
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        lapackf77_dlascl( uplo_, &izero, &izero, &d_one, &sigma, &n, &n, A, &lda, info );
    }

    /* Call DSYTRD to reduce symmetric matrix to tridiagonal form. */
    // dsytrd work: e (n) + tau (n) + llwork (n*nb)  ==>  2n + n*nb
    // dstedx work: e (n) + tau (n) + z (n*n) + llwrk2 (1 + 4*n + n^2)  ==>  1 + 6n + 2n^2
    inde   = 0;
    indtau = inde   + n;
    indwrk = indtau + n;
    indwk2 = indwrk + n*n;
    llwork = lwork - indwrk;
    llwrk2 = lwork - indwk2;

    magma_tally3_timer_t time=0;
    timer_start( time );

    magma_tally3_dsytrd( uplo, n, A, lda, w, &work[inde],
                  &work[indtau], &work[indwrk], llwork, &iinfo );

    timer_stop( time );
    timer_printf( "time dsytrd = %6.2f\n", time );

    /* For eigenvalues only, call DSTERF.  For eigenvectors, first call
       DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
       tridiagonal matrix, then call DORMTR to multiply it to the Householder
       transformations represented as Householder vectors in A. */
    if (! wantz) {
        lapackf77_dsterf( &n, w, &work[inde], info );
    }
    else {
        timer_start( time );

        if (MAGMA_tally3_SUCCESS != magma_tally3_dmalloc( &dwork, 3*n*(n/2 + 1) )) {
            *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
            return *info;
        }

        // TTT Possible bug for n < 128
        magma_tally3_dstedx( Magma_tally3RangeAll, n, 0., 0., 0, 0, w, &work[inde],
                      &work[indwrk], n, &work[indwk2], llwrk2,
                      iwork, liwork, dwork, info );

        magma_tally3_free( dwork );

        timer_stop( time );
        timer_printf( "time dstedx = %6.2f\n", time );
        timer_start( time );

        magma_tally3_dormtr( Magma_tally3Left, uplo, Magma_tally3NoTrans, n, n, A, lda, &work[indtau],
                      &work[indwrk], n, &work[indwk2], llwrk2, &iinfo );

        lapackf77_dlacpy( "A", &n, &n, &work[indwrk], &n, A, &lda );

        timer_stop( time );
        timer_printf( "time dormtr + copy = %6.2f\n", time );
    }

    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1) {
        d__1 = 1. / sigma;
        blasf77_dscal( &n, &d__1, w, &ione );
    }

    work[0]  = lwmin * one_eps;  // round up
    iwork[0] = liwmin;

    return *info;
} /* magma_tally3_dsyevd */