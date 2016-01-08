/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov
       @author Azzam Haidar

       @precisions normal z -> c

*/
#include "common_magma_tally4.h"
#include "magma_tally4_timer.h"

#define PRECISION_z
#define COMPLEX

/**
    Purpose
    -------
    ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a
    complex Hermitian matrix A.  If eigenvectors are desired, it uses a
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
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    jobz    magma_tally4_vec_t
      -     = Magma_tally4NoVec:  Compute eigenvalues only;
      -     = Magma_tally4Vec:    Compute eigenvalues and eigenvectors.

    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of A is stored;
      -     = Magma_tally4Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       COMPLEX_16 array, dimension (LDA, N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally4Upper, the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = Magma_tally4Lower,
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.
            On exit, if JOBZ = Magma_tally4Vec, then if INFO = 0, A contains the
            orthonormal eigenvectors of the matrix A.
            If JOBZ = Magma_tally4NoVec, then on exit the lower triangle (if UPLO=Magma_tally4Lower)
            or the upper triangle (if UPLO=Magma_tally4Upper) of A, including the
            diagonal, is destroyed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

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
    rwork   (workspace) DOUBLE PRECISION array,
                                           dimension (LRWORK)
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
magma_tally4_zheevd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info)
{
    const char* uplo_ = lapack_uplo_const( uplo );
    const char* jobz_ = lapack_vec_const( jobz );
    magma_tally4_int_t ione = 1;
    magma_tally4_int_t izero = 0;
    double d_one = 1.;

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
    magma_tally4_int_t indwk2, llwrk2;
    magma_tally4_int_t iscale;
    double safmin;
    double bignum;
    magma_tally4_int_t indtau;
    magma_tally4_int_t indrwk, indwrk, liwmin;
    magma_tally4_int_t lrwmin, llwork;
    double smlnum;
    magma_tally4_int_t lquery;

    wantz = (jobz == Magma_tally4Vec);
    lower = (uplo == Magma_tally4Lower);
    lquery = (lwork == -1 || lrwork == -1 || liwork == -1);

    *info = 0;
    if (! (wantz || (jobz == Magma_tally4NoVec))) {
        *info = -1;
    } else if (! (lower || (uplo == Magma_tally4Upper))) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (lda < max(1,n)) {
        *info = -5;
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
    work[0]  = MAGMA_tally4_Z_MAKE( lwmin * one_eps, 0.);
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    if ((lwork < lwmin) && !lquery) {
        *info = -8;
    } else if ((lrwork < lrwmin) && ! lquery) {
        *info = -10;
    } else if ((liwork < liwmin) && ! lquery) {
        *info = -12;
    }

    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
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
        w[0] = MAGMA_tally4_Z_REAL(A[0]);
        if (wantz) {
            A[0] = MAGMA_tally4_Z_ONE;
        }
        return *info;
    }

    /* Check if matrix is very small then just call LAPACK on CPU, no need for GPU */
    if (n <= 128) {
        #ifdef ENABLE_DEBUG
        printf("--------------------------------------------------------------\n");
        printf("  warning matrix too small N=%d NB=%d, calling lapack on CPU  \n", (int) n, (int) nb);
        printf("--------------------------------------------------------------\n");
        #endif
        lapackf77_zheevd(jobz_, uplo_,
                         &n, A, &lda,
                         w, work, &lwork,
                         #if defined(PRECISION_z) || defined(PRECISION_c)
                         rwork, &lrwork,
                         #endif
                         iwork, &liwork, info);
        return *info;
    }

    /* Get machine constants. */
    safmin = lapackf77_dlamch("Safe minimum");
    eps = lapackf77_dlamch("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = magma_tally4_dsqrt(smlnum);
    rmax = magma_tally4_dsqrt(bignum);

    /* Scale matrix to allowable range, if necessary. */
    anrm = lapackf77_zlanhe("M", uplo_, &n, A, &lda, rwork);
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        lapackf77_zlascl(uplo_, &izero, &izero, &d_one, &sigma, &n, &n, A,
                         &lda, info);
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
    indwk2 = indwrk + n*n;
    llwork = lwork - indwrk;
    llwrk2 = lwork - indwk2;

    magma_tally4_timer_t time=0;
    timer_start( time );

    magma_tally4_zhetrd_mgpu(ngpu, 1, uplo, n, A, lda, w, &rwork[inde],
                      &work[indtau], &work[indwrk], llwork, &iinfo);

    timer_stop( time );
    timer_printf( "time zhetrd = %6.2f\n", time );

    /* For eigenvalues only, call DSTERF.  For eigenvectors, first call
       ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
       tridiagonal matrix, then call ZUNMTR to multiply it to the Householder
       transformations represented as Householder vectors in A. */
    if (! wantz) {
        lapackf77_dsterf(&n, w, &rwork[inde], info);
    }
    else {
        timer_start( time );

#ifdef USE_SINGLE_GPU
        if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc( &dwork, 3*n*(n/2 + 1) )) {
            *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
            return *info;
        }

        magma_tally4_zstedx(Magma_tally4RangeAll, n, 0, 0, 0, 0, w, &rwork[inde],
                     &work[indwrk], n, &rwork[indrwk],
                     llrwk, iwork, liwork, dwork, info);

        magma_tally4_free( dwork );
#else
        magma_tally4_zstedx_m(ngpu, Magma_tally4RangeAll, n, 0, 0, 0, 0, w, &rwork[inde],
                       &work[indwrk], n, &rwork[indrwk],
                       llrwk, iwork, liwork, info);
#endif

        timer_stop( time );
        timer_printf( "time zstedc = %6.2f\n", time );
        timer_start( time );

        magma_tally4_zunmtr_m(ngpu, Magma_tally4Left, uplo, Magma_tally4NoTrans, n, n, A, lda, &work[indtau],
                       &work[indwrk], n, &work[indwk2], llwrk2, &iinfo);

        lapackf77_zlacpy("A", &n, &n, &work[indwrk], &n, A, &lda);

        timer_stop( time );
        timer_printf( "time zunmtr + copy = %6.2f\n", time );
    }

    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1) {
        if (*info == 0) {
            imax = n;
        } else {
            imax = *info - 1;
        }
        d__1 = 1. / sigma;
        blasf77_dscal(&imax, &d__1, w, &ione);
    }

    work[0]  = MAGMA_tally4_Z_MAKE( lwmin * one_eps, 0.);  // round up
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    return *info;
} /* magma_tally4_zheevd_m */
