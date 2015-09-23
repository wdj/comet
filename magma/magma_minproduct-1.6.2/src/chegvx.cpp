/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
    
       @author Raffaele Solca
    
       @generated from zhegvx.cpp normal z -> c, Fri Jan 30 19:00:18 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    CHEGVX computes selected eigenvalues, and optionally, eigenvectors
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
    A       COMPLEX array, dimension (LDA, N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_minproductUpper, the
            leading N-by-N upper triangular part of A contains the
            upper triangular part of the matrix A.  If UPLO = Magma_minproductLower,
            the leading N-by-N lower triangular part of A contains
            the lower triangular part of the matrix A.
    \n
            On exit,  the lower triangle (if UPLO=Magma_minproductLower) or the upper
            triangle (if UPLO=Magma_minproductUpper) of A, including the diagonal, is
            destroyed.
    
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).
    
    @param[in,out]
    B       COMPLEX array, dimension (LDB, N)
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
    vl      REAL
    @param[in]
    vu      REAL
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
    
    @param[in]
    abstol  REAL
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
            set to twice the underflow threshold 2*SLAMCH('S'), not zero.
            If this routine returns with INFO > 0, indicating that some
            eigenvectors did not converge, try setting ABSTOL to
            2*SLAMCH('S').
    
    @param[out]
    m       INTEGER
            The total number of eigenvalues found.  0 <= M <= N.
            If RANGE = Magma_minproductRangeAll, M = N, and if RANGE = Magma_minproductRangeI, M = IU-IL+1.
    
    @param[out]
    w       REAL array, dimension (N)
            The first M elements contain the selected
            eigenvalues in ascending order.
    
    @param[out]
    Z       COMPLEX array, dimension (LDZ, max(1,M))
            If JOBZ = Magma_minproductNoVec, then Z is not referenced.
            If JOBZ = Magma_minproductVec, then if INFO = 0, the first M columns of Z
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
            supplied in the array Z; if RANGE = Magma_minproductRangeV, the exact value of M
            is not known in advance and an upper bound must be used.
    
    @param[in]
    ldz     INTEGER
            The leading dimension of the array Z.  LDZ >= 1, and if
            JOBZ = Magma_minproductVec, LDZ >= max(1,N).
    
    @param[out]
    work    (workspace) COMPLEX array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.
    
    @param[in]
    lwork   INTEGER
            The length of the array WORK.  LWORK >= max(1,2*N).
            For optimal efficiency, LWORK >= (NB+1)*N,
            where NB is the blocksize for CHETRD returned by ILAENV.
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.
    
    @param
    rwork   (workspace) REAL array, dimension (7*N)
    
    @param
    iwork   (workspace) INTEGER array, dimension (5*N)
    
    @param[out]
    ifail   INTEGER array, dimension (N)
            If JOBZ = Magma_minproductVec, then if INFO = 0, the first M elements of
            IFAIL are zero.  If INFO > 0, then IFAIL contains the
            indices of the eigenvectors that failed to converge.
            If JOBZ = Magma_minproductNoVec, then IFAIL is not referenced.
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  CPOTRF or CHEEVX returned an error code:
            <= N: if INFO = i, CHEEVX failed to converge;
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

    @ingroup magma_minproduct_chegv_driver
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_chegvx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu, float abstol,
    magma_minproduct_int_t *m, float *w,  magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork, float *rwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info)
{
    magma_minproductFloatComplex c_one = MAGMA_minproduct_C_ONE;
    
    magma_minproductFloatComplex *dA;
    magma_minproductFloatComplex *dB;
    magma_minproductFloatComplex *dZ;
    magma_minproduct_int_t ldda = n;
    magma_minproduct_int_t lddb = n;
    magma_minproduct_int_t lddz = n;
    
    magma_minproduct_int_t lower;
    magma_minproduct_trans_t trans;
    magma_minproduct_int_t wantz;
    magma_minproduct_int_t lquery;
    magma_minproduct_int_t alleig, valeig, indeig;
    
    magma_minproduct_int_t lwmin;
    
    magma_minproduct_queue_t stream;
    magma_minproduct_queue_create( &stream );
    
    wantz  = (jobz  == Magma_minproductVec);
    lower  = (uplo  == Magma_minproductLower);
    alleig = (range == Magma_minproductRangeAll);
    valeig = (range == Magma_minproductRangeV);
    indeig = (range == Magma_minproductRangeI);
    lquery = (lwork == -1);
    
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
    
    magma_minproduct_int_t nb = magma_minproduct_get_chetrd_nb(n);
    
    lwmin = n * (nb + 1);
    
    work[0] = MAGMA_minproduct_C_MAKE( lwmin, 0 );
    
    
    if (lwork < lwmin && ! lquery) {
        *info = -20;
    }
    
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info));
        return *info;
    } else if (lquery) {
        return *info;
    }
    
    /* Quick return if possible */
    if (n == 0) {
        return *info;
    }
    
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dA, n*ldda ) ||
        MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dB, n*lddb ) ||
        MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dZ, n*lddz )) {
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    /*     Form a Cholesky factorization of B. */
    
    magma_minproduct_csetmatrix( n, n, B, ldb, dB, lddb );
    
    magma_minproduct_csetmatrix_async( n, n,
                            A,  lda,
                            dA, ldda, stream );
    
    magma_minproduct_cpotrf_gpu(uplo, n, dB, lddb, info);
    if (*info != 0) {
        *info = n + *info;
        return *info;
    }
    
    magma_minproduct_queue_sync( stream );
    
    magma_minproduct_cgetmatrix_async( n, n,
                            dB, lddb,
                            B,  ldb, stream );
    
    /* Transform problem to standard eigenvalue problem and solve. */
    magma_minproduct_chegst_gpu(itype, uplo, n, dA, ldda, dB, lddb, info);
    magma_minproduct_cheevx_gpu(jobz, range, uplo, n, dA, ldda, vl, vu, il, iu, abstol, m, w, dZ, lddz, A, lda, Z, ldz, work, lwork, rwork, iwork, ifail, info);
    
    if (wantz && *info == 0) {
        /* Backtransform eigenvectors to the original problem. */
        if (itype == 1 || itype == 2) {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
               backtransform eigenvectors: x = inv(L)'*y or inv(U)*y */
            if (lower) {
                trans = Magma_minproductConjTrans;
            } else {
                trans = Magma_minproductNoTrans;
            }
            magma_minproduct_ctrsm(Magma_minproductLeft, uplo, trans, Magma_minproductNonUnit, n, *m, c_one, dB, lddb, dZ, lddz);
        }
        else if (itype == 3) {
            /* For B*A*x=(lambda)*x;
               backtransform eigenvectors: x = L*y or U'*y */
            if (lower) {
                trans = Magma_minproductNoTrans;
            } else {
                trans = Magma_minproductConjTrans;
            }
            magma_minproduct_ctrmm(Magma_minproductLeft, uplo, trans, Magma_minproductNonUnit, n, *m, c_one, dB, lddb, dZ, lddz);
        }
        
        magma_minproduct_cgetmatrix( n, *m, dZ, lddz, Z, ldz );
    }
    
    magma_minproduct_queue_sync( stream );
    magma_minproduct_queue_destroy( stream );
    
    magma_minproduct_free( dA );
    magma_minproduct_free( dB );
    magma_minproduct_free( dZ );
    
    return *info;
} /* magma_minproduct_chegvx */
