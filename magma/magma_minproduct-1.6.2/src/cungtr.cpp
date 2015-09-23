/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zungtr.cpp normal z -> c, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    CUNGTR generates a complex unitary matrix Q which is defined as the
    product of n-1 elementary reflectors of order N, as returned by
    CHETRD:

    if UPLO = Magma_minproductUpper, Q = H(n-1) . . . H(2) H(1),

    if UPLO = Magma_minproductLower, Q = H(1) H(2) . . . H(n-1).

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper: Upper triangle of A contains elementary reflectors
                          from CHETRD;
      -     = Magma_minproductLower: Lower triangle of A contains elementary reflectors
                          from CHETRD.

    @param[in]
    n       INTEGER
            The order of the matrix Q. N >= 0.

    @param[in,out]
    A       COMPLEX array, dimension (LDA,N)
            On entry, the vectors which define the elementary reflectors,
            as returned by CHETRD.
            On exit, the N-by-N unitary matrix Q.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A. LDA >= N.

    @param[in]
    tau     COMPLEX array, dimension (N-1)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by CHETRD.

    @param[out]
    work    (workspace) COMPLEX array, dimension (LWORK)
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK. LWORK >= N-1.
            For optimum performance LWORK >= N*NB, where NB is
            the optimal blocksize.
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    @param[in]
    dT      COMPLEX array on the GPU device.
            DT contains the T matrices used in blocking the elementary
            reflectors H(i) as returned by magma_minproduct_chetrd.

    @param[in]
    nb      INTEGER
            This is the block size used in CHETRD, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_cheev_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_cungtr(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproductFloatComplex *dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info)
{
#define A(i,j) (A + (j)*lda+ (i))

    magma_minproduct_int_t i__1;
    magma_minproduct_int_t i, j;
    magma_minproduct_int_t iinfo;
    magma_minproduct_int_t upper, lwkopt, lquery;

    *info = 0;
    lquery = (lwork == -1);
    upper = (uplo == Magma_minproductUpper);
    if (! upper && uplo != Magma_minproductLower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,n)) {
        *info = -4;
    } else /* if (complicated condition) */ {
        /* Computing MAX */
        if (lwork < max(1, n-1) && ! lquery) {
            *info = -7;
        }
    }

    lwkopt = max(1, n) * nb;
    if (*info == 0) {
        work[0] = MAGMA_minproduct_C_MAKE( lwkopt, 0 );
    }

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info));
        return *info;
    } else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (n == 0) {
        work[0] = MAGMA_minproduct_C_ONE;
        return *info;
    }

    if (upper) {
        /*  Q was determined by a call to CHETRD with UPLO = 'U'
            Shift the vectors which define the elementary reflectors one
            column to the left, and set the last row and column of Q to
            those of the unit matrix                                    */
        for (j = 0; j < n-1; ++j) {
            for (i = 0; i < j-1; ++i)
                *A(i, j) = *A(i, j + 1);

            *A(n-1, j) = MAGMA_minproduct_C_ZERO;
        }
        for (i = 0; i < n-1; ++i) {
            *A(i, n-1) = MAGMA_minproduct_C_ZERO;
        }
        *A(n-1, n-1) = MAGMA_minproduct_C_ONE;
        
        /* Generate Q(1:n-1,1:n-1) */
        i__1 = n - 1;
        lapackf77_cungql(&i__1, &i__1, &i__1, A(0,0), &lda, tau, work,
                         &lwork, &iinfo);
    } else {
        
        /*  Q was determined by a call to CHETRD with UPLO = 'L'.
            Shift the vectors which define the elementary reflectors one
            column to the right, and set the first row and column of Q to
            those of the unit matrix                                      */
        for (j = n-1; j > 0; --j) {
            *A(0, j) = MAGMA_minproduct_C_ZERO;
            for (i = j; i < n-1; ++i)
                *A(i, j) = *A(i, j - 1);
        }

        *A(0, 0) = MAGMA_minproduct_C_ONE;
        for (i = 1; i < n-1; ++i)
            *A(i, 0) = MAGMA_minproduct_C_ZERO;
        
        if (n > 1) {
            /* Generate Q(2:n,2:n) */
            magma_minproduct_cungqr(n-1, n-1, n-1, A(1, 1), lda, tau, dT, nb, &iinfo);
        }
    }
    
    work[0] = MAGMA_minproduct_C_MAKE( lwkopt, 0 );

    return *info;
} /* magma_minproduct_cungtr */

#undef A
