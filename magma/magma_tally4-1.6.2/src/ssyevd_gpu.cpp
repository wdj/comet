/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Raffaele Solca
       @author Mark Gates
       @author Azzam Haidar

       @generated from dsyevd_gpu.cpp normal d -> s, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_tally4.h"
#include "magma_tally4_timer.h"

#define REAL
#define FAST_SYMV

/**
    Purpose
    -------
    SSYEVD_GPU computes all eigenvalues and, optionally, eigenvectors of
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
    dA      REAL array on the GPU,
            dimension (LDDA, N).
            On entry, the symmetric matrix A.  If UPLO = Magma_tally4Upper, the
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
    ldda    INTEGER
            The leading dimension of the array DA.  LDDA >= max(1,N).

    @param[out]
    w       REAL array, dimension (N)
            If INFO = 0, the eigenvalues in ascending order.

    @param
    wA      (workspace) REAL array, dimension (LDWA, N)

    @param[in]
    ldwa    INTEGER
            The leading dimension of the array wA.  LDWA >= max(1,N).

    @param[out]
    work    (workspace) REAL array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The length of the array WORK.
            If N <= 1,                      LWORK >= 1.
            If JOBZ = Magma_tally4NoVec and N > 1, LWORK >= 2*N + N*NB.
            If JOBZ = Magma_tally4Vec   and N > 1, LWORK >= max( 2*N + N*NB, 1 + 6*N + 2*N**2 ).
            NB can be obtained through magma_tally4_get_ssytrd_nb(N).
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
            If N <= 1,                       LIWORK >= 1.
            If JOBZ = Magma_tally4NoVec and N > 1, LIWORK >= 1.
            If JOBZ  = Magma_tally4Vec   and N > 1, LIWORK >= 3 + 5*N.
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

    @ingroup magma_tally4_ssyev_driver
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_ssyevd_gpu(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *w,
    float *wA,  magma_tally4_int_t ldwa,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info)
{
    magma_tally4_int_t ione = 1;

    float d__1;

    float eps;
    magma_tally4_int_t inde;
    float anrm;
    float rmin, rmax;
    float sigma;
    magma_tally4_int_t iinfo, lwmin;
    magma_tally4_int_t lower;
    magma_tally4_int_t wantz;
    magma_tally4_int_t indwk2, llwrk2;
    magma_tally4_int_t iscale;
    float safmin;
    float bignum;
    magma_tally4_int_t indtau;
    magma_tally4_int_t indwrk, liwmin;
    magma_tally4_int_t llwork;
    float smlnum;
    magma_tally4_int_t lquery;

    magma_tally4Float_ptr dwork;
    magma_tally4_int_t lddc = ldda;

    wantz = (jobz == Magma_tally4Vec);
    lower = (uplo == Magma_tally4Lower);
    lquery = (lwork == -1 || liwork == -1);

    *info = 0;
    if (! (wantz || (jobz == Magma_tally4NoVec))) {
        *info = -1;
    } else if (! (lower || (uplo == Magma_tally4Upper))) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    }

    magma_tally4_int_t nb = magma_tally4_get_ssytrd_nb( n );
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
    real_Double_t one_eps = 1. + lapackf77_slamch("Epsilon");
    work[0]  = lwmin * one_eps;
    iwork[0] = liwmin;

    if ((lwork < lwmin) && !lquery) {
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

    /* If matrix is very small, then just call LAPACK on CPU, no need for GPU */
    if (n <= 128) {
        magma_tally4_int_t lda = n;
        float *A;
        magma_tally4_smalloc_cpu( &A, lda*n );
        magma_tally4_sgetmatrix( n, n, dA, ldda, A, lda );
        lapackf77_ssyevd( lapack_vec_const_tally4(jobz), lapack_uplo_const_tally4(uplo),
                          &n, A, &lda,
                          w, work, &lwork,
                          iwork, &liwork, info );
        magma_tally4_ssetmatrix( n, n, A, lda, dA, ldda );
        magma_tally4_free_cpu( A );
        return *info;
    }

    magma_tally4_queue_t stream;
    magma_tally4_queue_create( &stream );

    // ssytrd2_gpu requires ldda*ceildiv(n,64) + 2*ldda*nb
    // sormtr_gpu  requires lddc*n
    // slansy      requires n
    magma_tally4_int_t ldwork = max( ldda*ceildiv(n,64) + 2*ldda*nb, lddc*n );
    ldwork = max( ldwork, n );
    if ( wantz ) {
        // sstedx requires 3n^2/2
        ldwork = max( ldwork, 3*n*(n/2 + 1) );
    }
    if (MAGMA_tally4_SUCCESS != magma_tally4_smalloc( &dwork, ldwork )) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }

    /* Get machine constants. */
    safmin = lapackf77_slamch("Safe minimum");
    eps    = lapackf77_slamch("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = magma_tally4_ssqrt( smlnum );
    rmax = magma_tally4_ssqrt( bignum );

    /* Scale matrix to allowable range, if necessary. */
    anrm = magma_tally4blas_slansy( Magma_tally4MaxNorm, uplo, n, dA, ldda, dwork );
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
        magma_tally4blas_slascl( uplo, 0, 0, 1., sigma, n, n, dA, ldda, info );
    }

    /* Call SSYTRD to reduce symmetric matrix to tridiagonal form. */
    // ssytrd work: e (n) + tau (n) + llwork (n*nb)  ==>  2n + n*nb
    // sstedx work: e (n) + tau (n) + z (n*n) + llwrk2 (1 + 4*n + n^2)  ==>  1 + 6n + 2n^2
    inde   = 0;
    indtau = inde   + n;
    indwrk = indtau + n;
    indwk2 = indwrk + n*n;
    llwork = lwork - indwrk;
    llwrk2 = lwork - indwk2;

    magma_tally4_timer_t time=0;
    timer_start( time );

#ifdef FAST_SYMV
    magma_tally4_ssytrd2_gpu( uplo, n, dA, ldda, w, &work[inde],
                       &work[indtau], wA, ldwa, &work[indwrk], llwork,
                       dwork, ldwork, &iinfo );
#else
    magma_tally4_ssytrd_gpu(  uplo, n, dA, ldda, w, &work[inde],
                       &work[indtau], wA, ldwa, &work[indwrk], llwork,
                       &iinfo );
#endif

    timer_stop( time );
    #ifdef FAST_SYMV
    timer_printf( "time ssytrd2 = %6.2f\n", time );
    #else
    timer_printf( "time ssytrd = %6.2f\n", time );
    #endif

    /* For eigenvalues only, call SSTERF.  For eigenvectors, first call
       SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
       tridiagonal matrix, then call SORMTR to multiply it to the Householder
       transformations represented as Householder vectors in A. */
    if (! wantz) {
        lapackf77_ssterf( &n, w, &work[inde], info );
    }
    else {
        timer_start( time );

        magma_tally4_sstedx( Magma_tally4RangeAll, n, 0., 0., 0, 0, w, &work[inde],
                      &work[indwrk], n, &work[indwk2],
                      llwrk2, iwork, liwork, dwork, info );

        timer_stop( time );
        timer_printf( "time sstedx = %6.2f\n", time );
        timer_start( time );

        magma_tally4_ssetmatrix( n, n, &work[indwrk], n, dwork, lddc );

        magma_tally4_sormtr_gpu( Magma_tally4Left, uplo, Magma_tally4NoTrans, n, n, dA, ldda, &work[indtau],
                          dwork, lddc, wA, ldwa, &iinfo );

        magma_tally4_scopymatrix( n, n, dwork, lddc, dA, ldda );

        timer_stop( time );
        timer_printf( "time sormtr + copy = %6.2f\n", time );
    }

    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1) {
        d__1 = 1. / sigma;
        blasf77_sscal( &n, &d__1, w, &ione );
    }

    work[0]  = lwmin * one_eps;  // round up
    iwork[0] = liwmin;

    magma_tally4_queue_destroy( stream );
    magma_tally4_free( dwork );

    return *info;
} /* magma_tally4_ssyevd_gpu */