/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
    
       @author Raffaele Solca
    
       @precisions normal z -> c

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    ZHEGVX computes selected eigenvalues, and optionally, eigenvectors
    of a complex generalized Hermitian-definite eigenproblem, of the form
    A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
    B are assumed to be Hermitian and B is also positive definite.
    Eigenvalues and eigenvectors can be selected by specifying either a
    range of values or a range of indices for the desired eigenvalues.
    
    Arguments
    ---------
    @param[in]
    itype   INTEGER
            Specifies the problem type to be solved:
            = 1:  A*x = (lambda)*B*x
            = 2:  A*B*x = (lambda)*x
            = 3:  B*A*x = (lambda)*x
    
    @param[in]
    jobz    magma_tally2_vec_t
      -     = Magma_tally2NoVec:  Compute eigenvalues only;
      -     = Magma_tally2Vec:    Compute eigenvalues and eigenvectors.
    
    @param[in]
    range   magma_tally2_range_t
      -     = Magma_tally2RangeAll: all eigenvalues will be found.
      -     = Magma_tally2RangeV:   all eigenvalues in the half-open interval (VL,VU]
                   will be found.
      -     = Magma_tally2RangeI:   the IL-th through IU-th eigenvalues will be found.
    
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangles of A and B are stored;
      -     = Magma_tally2Lower:  Lower triangles of A and B are stored.
    
    @param[in]
    n       INTEGER
            The order of the matrices A and B.  N >= 0.
    
    @param[in,out]
    A       COMPLEX_16 array, dimension (LDA, N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally2Upper, the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = Magma_tally2Lower,
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.
    \n
            On exit,  the lower triangle (if UPLO=Magma_tally2Lower) or the upper
            triangle (if UPLO=Magma_tally2Upper) of A, including the diagonal, is
            destroyed.
    
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).
    
    @param[in,out]
    B       COMPLEX_16 array, dimension (LDB, N)
            On entry, the Hermitian matrix B.  If UPLO = Magma_tally2Upper, the
            leading N-by-N upper triangular part of B contains the
            upper triangular part of the matrix B.  If UPLO = Magma_tally2Lower,
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
            If RANGE=Magma_tally2RangeV, the lower and upper bounds of the interval to
            be searched for eigenvalues. VL < VU.
            Not referenced if RANGE = Magma_tally2RangeAll or Magma_tally2RangeI.
    
    @param[in]
    il      INTEGER
    @param[in]
    iu      INTEGER
            If RANGE=Magma_tally2RangeI, the indices (in ascending order) of the
            smallest and largest eigenvalues to be returned.
            1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
            Not referenced if RANGE = Magma_tally2RangeAll or Magma_tally2RangeV.
    
    @param[in]
    abstol  DOUBLE PRECISION
            The absolute error tolerance for the eigenvalues.
            An approximate eigenvalue is accepted as converged
            when it is determined to lie in an interval [a,b]
            of width less than or equal to
    \n
                    ABSTOL + EPS * max( |a|,|b| ),
    \n
            where EPS is the machine precision.  If ABSTOL is less than
            or equal to zero, then  EPS*|T|  will be used in its place,
            where |T| is the 1-norm of the tridiagonal matrix obtained
            by reducing A to tridiagonal form.
    \n
            Eigenvalues will be computed most accurately when ABSTOL is
            set to twice the underflow threshold 2*DLAMCH('S'), not zero.
            If this routine returns with INFO > 0, indicating that some
            eigenvectors did not converge, try setting ABSTOL to
            2*DLAMCH('S').
    
    @param[out]
    m       INTEGER
            The total number of eigenvalues found.  0 <= M <= N.
            If RANGE = Magma_tally2RangeAll, M = N, and if RANGE = Magma_tally2RangeI, M = IU-IL+1.
    
    @param[out]
    w       DOUBLE PRECISION array, dimension (N)
            The first M elements contain the selected
            eigenvalues in ascending order.
    
    @param[out]
    Z       COMPLEX_16 array, dimension (LDZ, max(1,M))
            If JOBZ = Magma_tally2NoVec, then Z is not referenced.
            If JOBZ = Magma_tally2Vec, then if INFO = 0, the first M columns of Z
            contain the orthonormal eigenvectors of the matrix A
            corresponding to the selected eigenvalues, with the i-th
            column of Z holding the eigenvector associated with W(i).
            The eigenvectors are normalized as follows:
            if ITYPE = 1 or 2, Z**T*B*Z = I;
            if ITYPE = 3, Z**T*inv(B)*Z = I.
    \n
            If an eigenvector fails to converge, then that column of Z
            contains the latest approximation to the eigenvector, and the
            index of the eigenvector is returned in IFAIL.
            Note: the user must ensure that at least max(1,M) columns are
            supplied in the array Z; if RANGE = Magma_tally2RangeV, the exact value of M
            is not known in advance and an upper bound must be used.
    
    @param[in]
    ldz     INTEGER
            The leading dimension of the array Z.  LDZ >= 1, and if
            JOBZ = Magma_tally2Vec, LDZ >= max(1,N).
    
    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.
    
    @param[in]
    lwork   INTEGER
            The length of the array WORK.  LWORK >= max(1,2*N).
            For optimal efficiency, LWORK >= (NB+1)*N,
            where NB is the blocksize for ZHETRD returned by ILAENV.
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.
    
    @param
    rwork   (workspace) DOUBLE PRECISION array, dimension (7*N)
    
    @param
    iwork   (workspace) INTEGER array, dimension (5*N)
    
    @param[out]
    ifail   INTEGER array, dimension (N)
            If JOBZ = Magma_tally2Vec, then if INFO = 0, the first M elements of
            IFAIL are zero.  If INFO > 0, then IFAIL contains the
            indices of the eigenvectors that failed to converge.
            If JOBZ = Magma_tally2NoVec, then IFAIL is not referenced.
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  ZPOTRF or ZHEEVX returned an error code:
            <= N: if INFO = i, ZHEEVX failed to converge;
                  i eigenvectors failed to converge.  Their indices
                  are stored in array IFAIL.
            > N:  if INFO = N + i, for 1 <= i <= N, then the leading
                  minor of order i of B is not positive definite.
                  The factorization of B could not be completed and
                  no eigenvalues or eigenvectors were computed.
    
    Further Details
    ---------------
    Based on contributions by
       Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA

    @ingroup magma_tally2_zhegv_driver
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zhegvx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu, double abstol,
    magma_tally2_int_t *m, double *w,  magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork, double *rwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t *ifail,
    magma_tally2_int_t *info)
{
    magma_tally2DoubleComplex c_one = MAGMA_tally2_Z_ONE;
    
    magma_tally2DoubleComplex *dA;
    magma_tally2DoubleComplex *dB;
    magma_tally2DoubleComplex *dZ;
    magma_tally2_int_t ldda = n;
    magma_tally2_int_t lddb = n;
    magma_tally2_int_t lddz = n;
    
    magma_tally2_int_t lower;
    magma_tally2_trans_t trans;
    magma_tally2_int_t wantz;
    magma_tally2_int_t lquery;
    magma_tally2_int_t alleig, valeig, indeig;
    
    magma_tally2_int_t lwmin;
    
    magma_tally2_queue_t stream;
    magma_tally2_queue_create( &stream );
    
    wantz  = (jobz  == Magma_tally2Vec);
    lower  = (uplo  == Magma_tally2Lower);
    alleig = (range == Magma_tally2RangeAll);
    valeig = (range == Magma_tally2RangeV);
    indeig = (range == Magma_tally2RangeI);
    lquery = (lwork == -1);
    
    *info = 0;
    if (itype < 1 || itype > 3) {
        *info = -1;
    } else if (! (alleig || valeig || indeig)) {
        *info = -2;
    } else if (! (wantz || (jobz == Magma_tally2NoVec))) {
        *info = -3;
    } else if (! (lower || (uplo == Magma_tally2Upper))) {
        *info = -4;
    } else if (n < 0) {
        *info = -5;
    } else if (lda < max(1,n)) {
        *info = -7;
    } else if (ldb < max(1,n)) {
        *info = -9;
    } else if (ldz < 1 || (wantz && ldz < n)) {
        *info = -18;
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
    
    magma_tally2_int_t nb = magma_tally2_get_zhetrd_nb(n);
    
    lwmin = n * (nb + 1);
    
    work[0] = MAGMA_tally2_Z_MAKE( lwmin, 0 );
    
    
    if (lwork < lwmin && ! lquery) {
        *info = -20;
    }
    
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info));
        return *info;
    } else if (lquery) {
        return *info;
    }
    
    /* Quick return if possible */
    if (n == 0) {
        return *info;
    }
    
    if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dA, n*ldda ) ||
        MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dB, n*lddb ) ||
        MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dZ, n*lddz )) {
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    /*     Form a Cholesky factorization of B. */
    
    magma_tally2_zsetmatrix( n, n, B, ldb, dB, lddb );
    
    magma_tally2_zsetmatrix_async( n, n,
                            A,  lda,
                            dA, ldda, stream );
    
    magma_tally2_zpotrf_gpu(uplo, n, dB, lddb, info);
    if (*info != 0) {
        *info = n + *info;
        return *info;
    }
    
    magma_tally2_queue_sync( stream );
    
    magma_tally2_zgetmatrix_async( n, n,
                            dB, lddb,
                            B,  ldb, stream );
    
    /* Transform problem to standard eigenvalue problem and solve. */
    magma_tally2_zhegst_gpu(itype, uplo, n, dA, ldda, dB, lddb, info);
    magma_tally2_zheevx_gpu(jobz, range, uplo, n, dA, ldda, vl, vu, il, iu, abstol, m, w, dZ, lddz, A, lda, Z, ldz, work, lwork, rwork, iwork, ifail, info);
    
    if (wantz && *info == 0) {
        /* Backtransform eigenvectors to the original problem. */
        if (itype == 1 || itype == 2) {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
               backtransform eigenvectors: x = inv(L)'*y or inv(U)*y */
            if (lower) {
                trans = Magma_tally2ConjTrans;
            } else {
                trans = Magma_tally2NoTrans;
            }
            magma_tally2_ztrsm(Magma_tally2Left, uplo, trans, Magma_tally2NonUnit, n, *m, c_one, dB, lddb, dZ, lddz);
        }
        else if (itype == 3) {
            /* For B*A*x=(lambda)*x;
               backtransform eigenvectors: x = L*y or U'*y */
            if (lower) {
                trans = Magma_tally2NoTrans;
            } else {
                trans = Magma_tally2ConjTrans;
            }
            magma_tally2_ztrmm(Magma_tally2Left, uplo, trans, Magma_tally2NonUnit, n, *m, c_one, dB, lddb, dZ, lddz);
        }
        
        magma_tally2_zgetmatrix( n, *m, dZ, lddz, Z, ldz );
    }
    
    magma_tally2_queue_sync( stream );
    magma_tally2_queue_destroy( stream );
    
    magma_tally2_free( dA );
    magma_tally2_free( dB );
    magma_tally2_free( dZ );
    
    return *info;
} /* magma_tally2_zhegvx */