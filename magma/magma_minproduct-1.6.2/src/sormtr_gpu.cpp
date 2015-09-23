/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov

       @generated from zunmtr_gpu.cpp normal z -> s, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    SORMTR overwrites the general real M-by-N matrix C with

                                SIDE = Magma_minproductLeft    SIDE = Magma_minproductRight
    TRANS = Magma_minproductNoTrans:       Q * C               C * Q
    TRANS = Magma_minproductTrans:    Q**H * C            C * Q**H

    where Q is a real unitary matrix of order nq, with nq = m if
    SIDE = Magma_minproductLeft and nq = n if SIDE = Magma_minproductRight. Q is defined as the product of
    nq-1 elementary reflectors, as returned by SSYTRD:

    if UPLO = Magma_minproductUpper, Q = H(nq-1) . . . H(2) H(1);

    if UPLO = Magma_minproductLower, Q = H(1) H(2) . . . H(nq-1).

    Arguments
    ---------
    @param[in]
    side    magma_minproduct_side_t
      -     = Magma_minproductLeft:      apply Q or Q**H from the Left;
      -     = Magma_minproductRight:     apply Q or Q**H from the Right.

    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper: Upper triangle of A contains elementary reflectors
                   from SSYTRD;
      -     = Magma_minproductLower: Lower triangle of A contains elementary reflectors
                   from SSYTRD.

    @param[in]
    trans   magma_minproduct_trans_t
      -     = Magma_minproductNoTrans:    No transpose, apply Q;
      -     = Magma_minproductTrans: Conjugate transpose, apply Q**H.

    @param[in]
    m       INTEGER
            The number of rows of the matrix C. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C. N >= 0.

    @param[in]
    dA      REAL array, dimension
                                 (LDDA,M) if SIDE = Magma_minproductLeft
                                 (LDDA,N) if SIDE = Magma_minproductRight
            The vectors which define the elementary reflectors, as
            returned by SSYTRD_GPU. On output the diagonal, the subdiagonal and the
            upper part (UPLO=Magma_minproductLower) or lower part (UPLO=Magma_minproductUpper) are destroyed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            LDDA >= max(1,M) if SIDE = Magma_minproductLeft; LDDA >= max(1,N) if SIDE = Magma_minproductRight.

    @param[in]
    tau     REAL array, dimension
                                 (M-1) if SIDE = Magma_minproductLeft
                                 (N-1) if SIDE = Magma_minproductRight
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SSYTRD.

    @param[in,out]
    dC      REAL array, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by (Q*C) or (Q**H * C) or (C * Q**H) or (C*Q).

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDDC >= max(1,M).

    @param[in]
    wA      (workspace) REAL array, dimension
                                 (LDWA,M) if SIDE = Magma_minproductLeft
                                 (LDWA,N) if SIDE = Magma_minproductRight
            The vectors which define the elementary reflectors, as
            returned by SSYTRD_GPU.

    @param[in]
    ldwa    INTEGER
            The leading dimension of the array wA.
            LDWA >= max(1,M) if SIDE = Magma_minproductLeft; LDWA >= max(1,N) if SIDE = Magma_minproductRight.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_ssyev_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_sormtr_gpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA,    magma_minproduct_int_t ldda,
    float   *tau,
    magma_minproductFloat_ptr dC,    magma_minproduct_int_t lddc,
    float    *wA,    magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    #define wA(i_,j_) (wA + (i_) + (j_)*ldwa)
    
    magma_minproduct_int_t i1, i2, mi, ni, nq;
    int left, upper;
    magma_minproduct_int_t iinfo;

    *info = 0;
    left   = (side == Magma_minproductLeft);
    upper  = (uplo == Magma_minproductUpper);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        //nw = n;
    } else {
        nq = n;
        //nw = m;
    }
    if (! left && side != Magma_minproductRight) {
        *info = -1;
    } else if (! upper && uplo != Magma_minproductLower) {
        *info = -2;
    } else if (trans != Magma_minproductNoTrans &&
               trans != Magma_minproductTrans) {
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
        magma_minproduct_xerbla( __func__, -(*info) );
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
        magma_minproduct_sormql2_gpu(side, trans, mi, ni, nq-1, dA(0,1), ldda, tau,
                          dC, lddc, wA(0,1), ldwa, &iinfo);
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
        magma_minproduct_sormqr2_gpu(side, trans, mi, ni, nq-1, dA(1,0), ldda, tau,
                          dC(i1,i2), lddc, wA(1,0), ldwa, &iinfo);
    }

    return *info;
} /* magma_minproduct_sormtr */
