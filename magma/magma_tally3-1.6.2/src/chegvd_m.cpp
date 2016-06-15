/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar
       @author Stan Tomov
       @author Mark Gates

       @generated from zhegvd_m.cpp normal z -> c, Fri Jan 30 19:00:18 2015

*/
#include "common_magma_tally3.h"
#include "magma_tally3_timer.h"

#define COMPLEX

/**
    Purpose
    -------
    CHEGVD computes all the eigenvalues, and optionally, the eigenvectors
    of a complex generalized Hermitian-definite eigenproblem, of the form
    A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
    B are assumed to be Hermitian and B is also positive definite.
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
    itype   INTEGER
            Specifies the problem type to be solved:
            = 1:  A*x = (lambda)*B*x
            = 2:  A*B*x = (lambda)*x
            = 3:  B*A*x = (lambda)*x

    @param[in]
    jobz    magma_tally3_vec_t
      -     = Magma_tally3NoVec:  Compute eigenvalues only;
      -     = Magma_tally3Vec:    Compute eigenvalues and eigenvectors.

    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  Upper triangles of A and B are stored;
      -     = Magma_tally3Lower:  Lower triangles of A and B are stored.

    @param[in]
    n       INTEGER
            The order of the matrices A and B.  N >= 0.

    @param[in,out]
    A       COMPLEX array, dimension (LDA, N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally3Upper, the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = Magma_tally3Lower,
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.
    \n
            On exit, if JOBZ = Magma_tally3Vec, then if INFO = 0, A contains the
            matrix Z of eigenvectors.  The eigenvectors are normalized
            as follows:
            if ITYPE = 1 or 2, Z**H*B*Z = I;
            if ITYPE = 3, Z**H*inv(B)*Z = I.
            If JOBZ = Magma_tally3NoVec, then on exit the upper triangle (if UPLO=Magma_tally3Upper)
            or the lower triangle (if UPLO=Magma_tally3Lower) of A, including the
            diagonal, is destroyed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in,out]
    B       COMPLEX array, dimension (LDB, N)
            On entry, the Hermitian matrix B.  If UPLO = Magma_tally3Upper, the
            leading N-by-N upper triangular part of B contains the
            upper triangular part of the matrix B.  If UPLO = Magma_tally3Lower,
            the leading N-by-N lower triangular part of B contains
            the lower triangular part of the matrix B.
    \n
            On exit, if INFO <= N, the part of B containing the matrix is
            overwritten by the triangular factor U or L from the Cholesky
            factorization B = U**H*U or B = L*L**H.

    @param[in]
    ldb     INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    w       REAL array, dimension (N)
            If INFO = 0, the eigenvalues in ascending order.

    @param[out]
    work    (workspace) COMPLEX array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The length of the array WORK.
            If N <= 1,                      LWORK >= 1.
            If JOBZ = Magma_tally3NoVec and N > 1, LWORK >= N + N*NB.
            If JOBZ = Magma_tally3Vec   and N > 1, LWORK >= max( N + N*NB, 2*N + N**2 ).
            NB can be obtained through magma_tally3_get_chetrd_nb(N).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal sizes of the WORK, RWORK and
            IWORK arrays, returns these values as the first entries of
            the WORK, RWORK and IWORK arrays, and no error message
            related to LWORK or LRWORK or LIWORK is issued by XERBLA.

    @param[out]
    rwork   (workspace) REAL array, dimension (MAX(1,LRWORK))
            On exit, if INFO = 0, RWORK[0] returns the optimal LRWORK.

    @param[in]
    lrwork  INTEGER
            The dimension of the array RWORK.
            If N <= 1,                      LRWORK >= 1.
            If JOBZ = Magma_tally3NoVec and N > 1, LRWORK >= N.
            If JOBZ = Magma_tally3Vec   and N > 1, LRWORK >= 1 + 5*N + 2*N**2.
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
            If JOBZ = Magma_tally3NoVec and N > 1, LIWORK >= 1.
            If JOBZ = Magma_tally3Vec   and N > 1, LIWORK >= 3 + 5*N.
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
      -     > 0:  CPOTRF or CHEEVD returned an error code:
               <= N:  if INFO = i and JOBZ = Magma_tally3NoVec, then the algorithm
                      failed to converge; i off-diagonal elements of an
                      intermediate tridiagonal form did not converge to
                      zero;
                      if INFO = i and JOBZ = Magma_tally3Vec, then the algorithm
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

    Modified so that no backsubstitution is performed if CHEEVD fails to
    converge (NEIG in old code could be greater than N causing out of
    bounds reference to A - reported by Ralf Meyer).  Also corrected the
    description of INFO and the test on ITYPE. Sven, 16 Feb 05.

    @ingroup magma_tally3_chegv_driver
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_chegvd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float *w, magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info)
{
    const char* uplo_ = lapack_uplo_const_tally3( uplo );
    const char* jobz_ = lapack_vec_const_tally3( jobz );

    magma_tally3FloatComplex c_one = MAGMA_tally3_C_ONE;

    magma_tally3_int_t lower;
    magma_tally3_trans_t trans;
    magma_tally3_int_t wantz;
    magma_tally3_int_t lquery;

    magma_tally3_int_t lwmin;
    magma_tally3_int_t liwmin;
    magma_tally3_int_t lrwmin;

    magma_tally3_queue_t stream;
    magma_tally3_queue_create( &stream );

    wantz = (jobz == Magma_tally3Vec);
    lower = (uplo == Magma_tally3Lower);
    lquery = (lwork == -1 || lrwork == -1 || liwork == -1);

    *info = 0;
    if (itype < 1 || itype > 3) {
        *info = -1;
    } else if (! (wantz || (jobz == Magma_tally3NoVec))) {
        *info = -2;
    } else if (! (lower || (uplo == Magma_tally3Upper))) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (lda < max(1,n)) {
        *info = -6;
    } else if (ldb < max(1,n)) {
        *info = -8;
    }

    magma_tally3_int_t nb = magma_tally3_get_chetrd_nb( n );
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
    real_Double_t one_eps = 1. + lapackf77_slamch("Epsilon");
    work[0]  = MAGMA_tally3_C_MAKE( lwmin * one_eps, 0.);  // round up
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    if (lwork < lwmin && ! lquery) {
        *info = -11;
    } else if (lrwork < lrwmin && ! lquery) {
        *info = -13;
    } else if (liwork < liwmin && ! lquery) {
        *info = -15;
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

    /* If matrix is very small, then just call LAPACK on CPU, no need for GPU */
    if (n <= 128) {
        lapackf77_chegvd( &itype, jobz_, uplo_,
                          &n, A, &lda, B, &ldb,
                          w, work, &lwork,
                          #ifdef COMPLEX
                          rwork, &lrwork,
                          #endif
                          iwork, &liwork, info);
        return *info;
    }

    magma_tally3_timer_t time=0;
    timer_start( time );

    magma_tally3_cpotrf_m( ngpu, uplo, n, B, ldb, info );
    if (*info != 0) {
        *info = n + *info;
        return *info;
    }

    timer_stop( time );
    timer_printf( "time cpotrf = %6.2f\n", time );
    timer_start( time );

    /*  Transform problem to standard eigenvalue problem and solve. */
    magma_tally3_chegst_m( ngpu, itype, uplo, n, A, lda, B, ldb, info );

    timer_stop( time );
    timer_printf( "time chegst = %6.2f\n", time );
    timer_start( time );

    magma_tally3_cheevd_m( ngpu, jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info );

    timer_stop( time );
    timer_printf( "time cheevd = %6.2f\n", time );

    if (wantz && *info == 0) {
        timer_start( time );

        /* Backtransform eigenvectors to the original problem. */
        if (itype == 1 || itype == 2) {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
               backtransform eigenvectors: x = inv(L)'*y or inv(U)*y */
            if (lower) {
                trans = Magma_tally3ConjTrans;
            } else {
                trans = Magma_tally3NoTrans;
            }
            magma_tally3_ctrsm_m( ngpu, Magma_tally3Left, uplo, trans, Magma_tally3NonUnit,
                           n, n, c_one, B, ldb, A, lda );
        }
        else if (itype == 3) {
            /* For B*A*x=(lambda)*x;
               backtransform eigenvectors: x = L*y or U'*y */
            if (lower) {
                trans = Magma_tally3NoTrans;
            } else {
                trans = Magma_tally3ConjTrans;
            }
            #ifdef ENABLE_DEBUG
            printf("--- the multi GPU version is falling back to 1 GPU to perform the last TRMM since there is no TRMM_mgpu --- \n");
            #endif
            magma_tally3FloatComplex *dA=NULL, *dB=NULL;
            magma_tally3_int_t ldda = roundup( n, 32 );
            magma_tally3_int_t lddb = ldda;
            
            if (MAGMA_tally3_SUCCESS != magma_tally3_cmalloc( &dA, n*ldda ) ||
                MAGMA_tally3_SUCCESS != magma_tally3_cmalloc( &dB, n*lddb )) {
                magma_tally3_free( dA );
                magma_tally3_free( dB );
                *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
                return *info;
            }
            magma_tally3_csetmatrix( n, n, B, ldb, dB, lddb );
            magma_tally3_csetmatrix( n, n, A, lda, dA, ldda );
            magma_tally3_ctrmm( Magma_tally3Left, uplo, trans, Magma_tally3NonUnit,
                         n, n, c_one, dB, lddb, dA, ldda );
            magma_tally3_cgetmatrix( n, n, dA, ldda, A, lda );
            
            magma_tally3_free( dA );
            magma_tally3_free( dB );
        }

        timer_stop( time );
        timer_printf( "time setmatrices trsm/mm + getmatrices = %6.2f\n", time );
    }

    work[0]  = MAGMA_tally3_C_MAKE( lwmin * one_eps, 0 );  // round up
    rwork[0] = lrwmin * one_eps;
    iwork[0] = liwmin;

    return *info;
} /* magma_tally3_chegvd_m */