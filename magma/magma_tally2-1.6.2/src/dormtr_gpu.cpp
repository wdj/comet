/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov

       @generated from zunmtr_gpu.cpp normal z -> d, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    DORMTR overwrites the general real M-by-N matrix C with

                                SIDE = Magma_tally2Left    SIDE = Magma_tally2Right
    TRANS = Magma_tally2NoTrans:       Q * C               C * Q
    TRANS = Magma_tally2Trans:    Q**H * C            C * Q**H

    where Q is a real unitary matrix of order nq, with nq = m if
    SIDE = Magma_tally2Left and nq = n if SIDE = Magma_tally2Right. Q is defined as the product of
    nq-1 elementary reflectors, as returned by DSYTRD:

    if UPLO = Magma_tally2Upper, Q = H(nq-1) . . . H(2) H(1);

    if UPLO = Magma_tally2Lower, Q = H(1) H(2) . . . H(nq-1).

    Arguments
    ---------
    @param[in]
    side    magma_tally2_side_t
      -     = Magma_tally2Left:      apply Q or Q**H from the Left;
      -     = Magma_tally2Right:     apply Q or Q**H from the Right.

    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper: Upper triangle of A contains elementary reflectors
                   from DSYTRD;
      -     = Magma_tally2Lower: Lower triangle of A contains elementary reflectors
                   from DSYTRD.

    @param[in]
    trans   magma_tally2_trans_t
      -     = Magma_tally2NoTrans:    No transpose, apply Q;
      -     = Magma_tally2Trans: Conjugate transpose, apply Q**H.

    @param[in]
    m       INTEGER
            The number of rows of the matrix C. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C. N >= 0.

    @param[in]
    dA      DOUBLE_PRECISION array, dimension
                                 (LDDA,M) if SIDE = Magma_tally2Left
                                 (LDDA,N) if SIDE = Magma_tally2Right
            The vectors which define the elementary reflectors, as
            returned by DSYTRD_GPU. On output the diagonal, the subdiagonal and the
            upper part (UPLO=Magma_tally2Lower) or lower part (UPLO=Magma_tally2Upper) are destroyed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            LDDA >= max(1,M) if SIDE = Magma_tally2Left; LDDA >= max(1,N) if SIDE = Magma_tally2Right.

    @param[in]
    tau     DOUBLE_PRECISION array, dimension
                                 (M-1) if SIDE = Magma_tally2Left
                                 (N-1) if SIDE = Magma_tally2Right
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by DSYTRD.

    @param[in,out]
    dC      DOUBLE_PRECISION array, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by (Q*C) or (Q**H * C) or (C * Q**H) or (C*Q).

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDDC >= max(1,M).

    @param[in]
    wA      (workspace) DOUBLE_PRECISION array, dimension
                                 (LDWA,M) if SIDE = Magma_tally2Left
                                 (LDWA,N) if SIDE = Magma_tally2Right
            The vectors which define the elementary reflectors, as
            returned by DSYTRD_GPU.

    @param[in]
    ldwa    INTEGER
            The leading dimension of the array wA.
            LDWA >= max(1,M) if SIDE = Magma_tally2Left; LDWA >= max(1,N) if SIDE = Magma_tally2Right.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_dsyev_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_dormtr_gpu(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA,    magma_tally2_int_t ldda,
    double   *tau,
    magma_tally2Double_ptr dC,    magma_tally2_int_t lddc,
    double    *wA,    magma_tally2_int_t ldwa,
    magma_tally2_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    #define wA(i_,j_) (wA + (i_) + (j_)*ldwa)
    
    magma_tally2_int_t i1, i2, mi, ni, nq;
    int left, upper;
    magma_tally2_int_t iinfo;

    *info = 0;
    left   = (side == Magma_tally2Left);
    upper  = (uplo == Magma_tally2Upper);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        //nw = n;
    } else {
        nq = n;
        //nw = m;
    }
    if (! left && side != Magma_tally2Right) {
        *info = -1;
    } else if (! upper && uplo != Magma_tally2Lower) {
        *info = -2;
    } else if (trans != Magma_tally2NoTrans &&
               trans != Magma_tally2Trans) {
        *info = -3;
    } else if (m < 0) {
        *info = -4;
    } else if (n < 0) {
        *info = -5;
    } else if (ldda < max(1,nq)) {
        *info = -7;
    } else if (lddc < max(1,m)) {
        *info = -10;
    } else if (ldwa < max(1,nq)) {
        *info = -12;
    }

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || nq == 1) {
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
        magma_tally2_dormql2_gpu(side, trans, mi, ni, nq-1, dA(0,1), ldda, tau,
                          dC, lddc, wA(0,1), ldwa, &iinfo);
    }
    else {
        /* Q was determined by a call to DSYTRD with UPLO = 'L' */
        if (left) {
            i1 = 1;
            i2 = 0;
        } else {
            i1 = 0;
            i2 = 1;
        }
        magma_tally2_dormqr2_gpu(side, trans, mi, ni, nq-1, dA(1,0), ldda, tau,
                          dC(i1,i2), lddc, wA(1,0), ldwa, &iinfo);
    }

    return *info;
} /* magma_tally2_dormtr */
