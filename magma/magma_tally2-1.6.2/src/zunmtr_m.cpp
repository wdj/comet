/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Raffaele Solca

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    ZUNMTR overwrites the general complex M-by-N matrix C with

                                SIDE = Magma_tally2Left    SIDE = Magma_tally2Right
    TRANS = Magma_tally2NoTrans:       Q * C               C * Q
    TRANS = Magma_tally2_ConjTrans:    Q**H * C            C * Q**H

    where Q is a complex unitary matrix of order nq, with nq = m if
    SIDE = Magma_tally2Left and nq = n if SIDE = Magma_tally2Right. Q is defined as the product of
    nq-1 elementary reflectors, as returned by SSYTRD:

    if UPLO = Magma_tally2Upper, Q = H(nq-1) . . . H(2) H(1);

    if UPLO = Magma_tally2Lower, Q = H(1) H(2) . . . H(nq-1).

    Arguments
    ---------
    @param[in]
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    side    magma_tally2_side_t
      -     = Magma_tally2Left:      apply Q or Q**H from the Left;
      -     = Magma_tally2Right:     apply Q or Q**H from the Right.

    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper: Upper triangle of A contains elementary reflectors
                   from SSYTRD;
      -     = Magma_tally2Lower: Lower triangle of A contains elementary reflectors
                   from SSYTRD.

    @param[in]
    trans   magma_tally2_trans_t
      -     = Magma_tally2NoTrans:    No transpose, apply Q;
      -     = Magma_tally2_ConjTrans: Conjugate transpose, apply Q**H.

    @param[in]
    m       INTEGER
            The number of rows of the matrix C. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C. N >= 0.

    @param[in]
    A       COMPLEX_16 array, dimension
                                 (LDA,M) if SIDE = Magma_tally2Left
                                 (LDA,N) if SIDE = Magma_tally2Right
            The vectors which define the elementary reflectors, as
            returned by SSYTRD.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.
            LDA >= max(1,M) if SIDE = Magma_tally2Left; LDA >= max(1,N) if SIDE = Magma_tally2Right.

    @param[in]
    tau     COMPLEX_16 array, dimension
                                 (M-1) if SIDE = Magma_tally2Left
                                 (N-1) if SIDE = Magma_tally2Right
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SSYTRD.

    @param[in,out]
    C       COMPLEX_16 array, dimension (LDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

    @param[in]
    ldc     INTEGER
            The leading dimension of the array C. LDC >= max(1,M).

    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.
            If SIDE = Magma_tally2Left,  LWORK >= max(1,N);
            if SIDE = Magma_tally2Right, LWORK >= max(1,M).
            For optimum performance LWORK >= N*NB if SIDE = Magma_tally2Left, and
            LWORK >= M*NB if SIDE = Magma_tally2Right, where NB is the optimal
            blocksize.
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_zheev_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zunmtr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *C,    magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info)
{
    #define A(i_,j_) (A + (i_) + (j_)*lda)
    #define C(i_,j_) (C + (i_) + (j_)*ldc)
    
    magma_tally2DoubleComplex c_one = MAGMA_tally2_Z_ONE;

    magma_tally2_int_t  i__2;
    magma_tally2_int_t i1, i2, nb, mi, ni, nq, nw;
    int left, upper, lquery;
    magma_tally2_int_t iinfo;
    magma_tally2_int_t lwkopt;

    *info = 0;
    left   = (side == Magma_tally2Left);
    upper  = (uplo == Magma_tally2Upper);
    lquery = (lwork == -1);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if (! left && side != Magma_tally2Right) {
        *info = -1;
    } else if (! upper && uplo != Magma_tally2Lower) {
        *info = -2;
    } else if (trans != Magma_tally2NoTrans &&
               trans != Magma_tally2_ConjTrans) {
        *info = -3;
    } else if (m < 0) {
        *info = -4;
    } else if (n < 0) {
        *info = -5;
    } else if (lda < max(1,nq)) {
        *info = -7;
    } else if (ldc < max(1,m)) {
        *info = -10;
    } else if (lwork < max(1,nw) && ! lquery) {
        *info = -12;
    }

    nb = 32;
    lwkopt = max(1,nw) * nb;
    if (*info == 0) {
        work[0] = MAGMA_tally2_Z_MAKE( lwkopt, 0 );
    }
    
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || nq == 1) {
        work[0] = c_one;
        return *info;
    }

    if (left) {
        mi = m - 1;
        ni = n;
    } else {
        mi = m;
        ni = n - 1;
    }

    if (upper) {
        /* Q was determined by a call to SSYTRD with UPLO = 'U' */
        i__2 = nq - 1;
        // TODO: upper case is not yet implemented for multiple GPUs -- see above
        // for now use one GPU
        //lapackf77_zunmql(side_, trans_, &mi, &ni, &i__2, A(0,1), &lda,
        //                 tau, C, &ldc, work, &lwork, &iinfo);
        //magma_tally2_zunmql_m(ngpu, side, trans, mi, ni, i__2, A(0,1), lda, tau,
        //               C, ldc, work, lwork, &iinfo);
        magma_tally2_zunmql(side, trans, mi, ni, i__2, A(0,1), lda, tau,
                       C, ldc, work, lwork, &iinfo); 
    }
    else {
        /* Q was determined by a call to SSYTRD with UPLO = 'L' */
        if (left) {
            i1 = 1;
            i2 = 0;
        } else {
            i1 = 0;
            i2 = 1;
        }
        i__2 = nq - 1;
        magma_tally2_zunmqr_m(ngpu, side, trans, mi, ni, i__2, A(1,0), lda, tau,
                       C(i1,i2), ldc, work, lwork, &iinfo);
    }

    work[0] = MAGMA_tally2_Z_MAKE( lwkopt, 0 );

    return *info;
} /* magma_tally2_zunmtr */
