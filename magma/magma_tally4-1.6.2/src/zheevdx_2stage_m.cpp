/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Raffaele Solca
       @author Azzam Haidar

       @precisions normal z -> c

*/
#include "common_magma_tally4.h"
#include "magma_tally4_timer.h"
#include "magma_tally4_bulge.h"
#include "magma_tally4_zbulge.h"

#define PRECISION_z
#define COMPLEX

/**
    Purpose
    -------
    ZHEEVD_2STAGE computes all eigenvalues and, optionally, eigenvectors of a
    complex Hermitian matrix A. It uses a two-stage algorithm for the tridiagonalization.
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
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

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
    A       COMPLEX_16 array, dimension (LDA, N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally4Upper, the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = Magma_tally4Lower,
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.
            On exit, if JOBZ = Magma_tally4Vec, then if INFO = 0, the first m columns
            of A contains the required
            orthonormal eigenvectors of the matrix A.
            If JOBZ = Magma_tally4NoVec, then on exit the lower triangle (if UPLO=Magma_tally4Lower)
            or the upper triangle (if UPLO=Magma_tally4Upper) of A, including the
            diagonal, is destroyed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

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
    m       INTEGER
            The total number of eigenvalues found.  0 <= M <= N.
            If RANGE = Magma_tally4RangeAll, M = N, and if RANGE = Magma_tally4RangeI, M = IU-IL+1.

    @param[out]
    w       DOUBLE PRECISION array, dimension (N)
            If INFO = 0, the required m eigenvalues in ascending order.

    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The length of the array WORK.
            If N <= 1,                      LWORK >= 1.
            If JOBZ = Magma_tally4NoVec and N > 1, LWORK >= LQ2 + N + N*NB.
            If JOBZ = Magma_tally4Vec   and N > 1, LWORK >= LQ2 + 2*N + N**2.
            where LQ2 is the size needed to store the Q2 matrix
            and is returned by magma_tally4_bulge_get_lq2.
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
magma_tally4_zheevdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info)
{
    #define A( i_,j_) (A  + (i_) + (j_)*lda)
    #define A2(i_,j_) (A2 + (i_) + (j_)*lda2)
    
    const char* uplo_  = lapack_uplo_const_tally4( uplo  );
    const char* jobz_  = lapack_vec_const_tally4( jobz  );
    magma_tally4DoubleComplex c_one  = MAGMA_tally4_Z_ONE;
    double d_one = 1.;
    magma_tally4_int_t ione = 1;
    magma_tally4_int_t izero = 0;

    double d__1;

    double eps;
    double anrm;
    magma_tally4_int_t imax;
    double rmin, rmax;
    double sigma;
    //magma_tally4_int_t iinfo;
    magma_tally4_int_t lwmin, lrwmin, liwmin;
    magma_tally4_int_t lower;
    magma_tally4_int_t wantz;
    magma_tally4_int_t iscale;
    double safmin;
    double bignum;
    double smlnum;
    magma_tally4_int_t lquery;
    magma_tally4_int_t alleig, valeig, indeig;
    magma_tally4_int_t len;

    /* determine the number of threads */
    magma_tally4_int_t parallel_threads = magma_tally4_get_parallel_numthreads();

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
    } else if (lda < max(1,n)) {
        *info = -6;
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

    magma_tally4_int_t nb = magma_tally4_get_zbulge_nb(n, parallel_threads);
    magma_tally4_int_t Vblksiz = magma_tally4_zbulge_get_Vblksiz(n, nb, parallel_threads);

    magma_tally4_int_t ldt = Vblksiz;
    magma_tally4_int_t ldv = nb + Vblksiz;
    magma_tally4_int_t blkcnt = magma_tally4_bulge_get_blkcnt(n, nb, Vblksiz);
    magma_tally4_int_t lq2 = magma_tally4_zbulge_get_lq2(n, parallel_threads);

    if (wantz) {
        lwmin  = lq2 + 2*n + n*n;
        lrwmin = 1 + 5*n + 2*n*n;
        liwmin = 5*n + 3;
    } else {
        lwmin  = lq2 + n + n*nb;
        lrwmin = n;
        liwmin = 1;
    }

    // multiply by 1+eps (in Double!) to ensure length gets rounded up,
    // if it cannot be exactly represented in floating point.
    real_Double_t one_eps = 1. + lapackf77_dlamch("Epsilon");
    work[0]  = MAGMA_tally4_Z_MAKE( lwmin * one_eps, 0.);  // round up
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    if ((lwork < lwmin) && !lquery) {
        *info = -14;
    } else if ((lrwork < lrwmin) && ! lquery) {
        *info = -16;
    } else if ((liwork < liwmin) && ! lquery) {
        *info = -18;
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

    magma_tally4_device_t orig_dev;
    magma_tally4_getdevice( &orig_dev );
    
    timer_printf("using %d parallel_threads\n", (int) parallel_threads);

    /* Check if matrix is very small then just call LAPACK on CPU, no need for GPU */
    magma_tally4_int_t ntiles = n/nb;
    if ( ( ntiles < 2 ) || ( n <= 128 ) ) {
        #ifdef ENABLE_DEBUG
        printf("--------------------------------------------------------------\n");
        printf("  warning matrix too small N=%d NB=%d, calling lapack on CPU  \n", (int) n, (int) nb);
        printf("--------------------------------------------------------------\n");
        #endif
        lapackf77_zheevd(jobz_, uplo_, &n,
                         A, &lda, w,
                         work, &lwork,
                         #if defined(PRECISION_z) || defined(PRECISION_c)
                         rwork, &lrwork,
                         #endif
                         iwork, &liwork,
                         info);
        *m = n;
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

    magma_tally4_int_t indT2   = 0;
    magma_tally4_int_t indTAU2 = indT2  + blkcnt*ldt*Vblksiz;
    magma_tally4_int_t indV2   = indTAU2+ blkcnt*Vblksiz;
    magma_tally4_int_t indtau1 = indV2  + blkcnt*ldv*Vblksiz;
    magma_tally4_int_t indwrk  = indtau1+ n;
    magma_tally4_int_t indwk2  = indwrk + n*n;
    magma_tally4_int_t llwork = lwork - indwrk;
    magma_tally4_int_t llwrk2 = lwork - indwk2;
    magma_tally4_int_t inde = 0;
    magma_tally4_int_t indrwk = inde + n;
    magma_tally4_int_t llrwk = lrwork - indrwk;

    magma_tally4_timer_t time=0, time_total=0, time_alloc=0, time_dist=0, time_band=0;
    timer_start( time_total );

#ifdef HE2HB_SINGLEGPU
    magma_tally4DoubleComplex *dT1;
    if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc( &dT1, n*nb)) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }
    timer_start( time_band );
    magma_tally4_zhetrd_he2hb(uplo, n, nb, A, lda, &work[indtau1], &work[indwrk], llwork, dT1, info);
    timer_stop( time_band );
    timer_printf( "    1 GPU seq code time zhetrd_he2hb only = %7.4f\n", time_band );
    magma_tally4_free(dT1);
#else
    magma_tally4_int_t nstream = max(3,ngpu+2);
    magma_tally4_queue_t streams[Magma_tally4MaxGPUs][20];
    magma_tally4DoubleComplex *da[Magma_tally4MaxGPUs], *dT1[Magma_tally4MaxGPUs];
    magma_tally4_int_t ldda = ((n+31)/32)*32;

    magma_tally4_int_t ver = 0;
    magma_tally4_int_t distblk = max(256, 4*nb);

    #ifdef ENABLE_DEBUG
    printf("voici ngpu %d distblk %d NB %d nstream %d version %d \n ", ngpu, distblk, nb, nstream, ver);
    #endif

    timer_start( time_alloc );
    for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally4_int_t mlocal = ((n / distblk) / ngpu + 1) * distblk;
        magma_tally4_setdevice( dev );
        // TODO check malloc
        magma_tally4_zmalloc(&da[dev], ldda*mlocal );
        magma_tally4_zmalloc(&dT1[dev], (n*nb) );
        for( int i = 0; i < nstream; ++i ) {
            magma_tally4_queue_create( &streams[dev][i] );
        }
    }
    timer_stop( time_alloc );
    
    timer_start( time_dist );
    magma_tally4_zsetmatrix_1D_col_bcyclic( n, n, A, lda, da, ldda, ngpu, distblk );
    magma_tally4_setdevice(0);
    timer_stop( time_dist );

    timer_start( time_band );
    if (ver == 30) {
        magma_tally4_zhetrd_he2hb_mgpu_spec(uplo, n, nb, A, lda, &work[indtau1], &work[indwrk], llwork, da, ldda, dT1, nb, ngpu, distblk, streams, nstream, info);
    } else {
        magma_tally4_zhetrd_he2hb_mgpu(uplo, n, nb, A, lda, &work[indtau1], &work[indwrk], llwork, da, ldda, dT1, nb, ngpu, distblk, streams, nstream, info);
    }
    timer_stop( time_band );
    timer_printf("    time alloc %7.4f, ditribution %7.4f, zhetrd_he2hb only = %7.4f\n", time_alloc, time_dist, time_band );

    for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally4_setdevice( dev );
        magma_tally4_free( da[dev] );
        magma_tally4_free( dT1[dev] );
        for( int i = 0; i < nstream; ++i ) {
            magma_tally4_queue_destroy( streams[dev][i] );
        }
    }
#endif // not HE2HB_SINGLEGPU

    timer_stop( time_total );
    timer_printf( "    time zhetrd_he2hb_mgpu = %6.2f\n", time_total );
    timer_start( time_total );
    timer_start( time );

    /* copy the input matrix into WORK(INDWRK) with band storage */
    /* PAY ATTENTION THAT work[indwrk] should be able to be of size lda2*n which it should be checked in any future modification of lwork.*/
    magma_tally4_int_t lda2 = 2*nb; //nb+1+(nb-1);
    magma_tally4DoubleComplex* A2 = &work[indwrk];
    memset(A2, 0, n*lda2*sizeof(magma_tally4DoubleComplex));

    for (magma_tally4_int_t j = 0; j < n-nb; j++) {
        len = nb+1;
        blasf77_zcopy( &len, A(j,j), &ione, A2(0,j), &ione );
        memset(A(j,j), 0, (nb+1)*sizeof(magma_tally4DoubleComplex));
        *A(nb+j,j) = c_one;
    }
    for (magma_tally4_int_t j = 0; j < nb; j++) {
        len = nb-j;
        blasf77_zcopy( &len, A(j+n-nb,j+n-nb), &ione, A2(0,j+n-nb), &ione );
        memset(A(j+n-nb,j+n-nb), 0, (nb-j)*sizeof(magma_tally4DoubleComplex));
    }

    timer_stop( time );
    timer_printf( "    time zhetrd_convert = %6.2f\n", time );
    timer_start( time );

    magma_tally4_zhetrd_hb2st(uplo, n, nb, Vblksiz, A2, lda2, w, &rwork[inde], &work[indV2], ldv, &work[indTAU2], wantz, &work[indT2], ldt);

    timer_stop( time );
    timer_stop( time_total );
    timer_printf( "    time zhetrd_hb2st = %6.2f\n", time );
    timer_printf( "  time zhetrd = %6.2f\n", time_total );

    /* For eigenvalues only, call DSTERF.  For eigenvectors, first call
       ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
       tridiagonal matrix, then call ZUNMTR to multiply it to the Householder
       transformations represented as Householder vectors in A. */
    if (! wantz) {
        timer_start( time );

        lapackf77_dsterf(&n, w, &rwork[inde], info);
        magma_tally4_dmove_eig(range, n, w, &il, &iu, vl, vu, m);

        timer_stop( time );
        timer_printf( "  time dstedc = %6.2f\n", time );
    }
    else {
        timer_start( time_total );
        timer_start( time );

        magma_tally4_zstedx_m(ngpu, range, n, vl, vu, il, iu, w, &rwork[inde],
                       &work[indwrk], n, &rwork[indrwk],
                       llrwk, iwork, liwork, info);

        timer_stop( time );
        timer_printf( "  time zstedx_m = %6.2f\n", time );
        timer_start( time );

        magma_tally4_dmove_eig(range, n, w, &il, &iu, vl, vu, m);
/*
        magma_tally4DoubleComplex *dZ;
        magma_tally4_int_t lddz = n;

        if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc( &dZ, *m*lddz)) {
            *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
            return *info;
        }

        magma_tally4_zbulge_back(uplo, n, nb, *m, Vblksiz, &work[indwrk + n * (il-1)], n, dZ, lddz,
                          &work[indV2], ldv, &work[indTAU2], &work[indT2], ldt, info);

        magma_tally4_zgetmatrix( n, *m, dZ, lddz, &work[indwrk], n);

        magma_tally4_free(dZ);

*/

        magma_tally4_zbulge_back_m(ngpu, uplo, n, nb, *m, Vblksiz, &work[indwrk + n * (il-1)], n,
                            &work[indV2], ldv, &work[indTAU2], &work[indT2], ldt, info);

        timer_stop( time );
        timer_printf( "    time zbulge_back_m = %6.2f\n", time );
        timer_start( time );

        magma_tally4_zunmqr_m(ngpu, Magma_tally4Left, Magma_tally4NoTrans, n-nb, *m, n-nb, A+nb, lda, &work[indtau1],
                       &work[indwrk + n * (il-1) + nb], n, &work[indwk2], llwrk2, info);

        lapackf77_zlacpy("A", &n, m, &work[indwrk  + n * (il-1)], &n, A, &lda);

        timer_stop( time );
        timer_stop( time_total );
        timer_printf( "    time zunmqr_m + copy = %6.2f\n", time );
        timer_printf( "  time eigenvectors backtransf. = %6.2f\n", time_total );
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

    magma_tally4_setdevice( orig_dev );
    
    return *info;
} /* magma_tally4_zheevdx_2stage_m */
