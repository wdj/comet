/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZUNMTR overwrites the general complex M-by-N matrix C with

                                SIDE = Magma_tally4Left    SIDE = Magma_tally4Right
    TRANS = Magma_tally4NoTrans:       Q * C               C * Q
    TRANS = Magma_tally4_ConjTrans:    Q**H * C            C * Q**H

    where Q is a complex unitary matrix of order nq, with nq = m if
    SIDE = Magma_tally4Left and nq = n if SIDE = Magma_tally4Right. Q is defined as the product of
    nq-1 elementary reflectors, as returned by ZHETRD:

    if UPLO = Magma_tally4Upper, Q = H(nq-1) . . . H(2) H(1);

    if UPLO = Magma_tally4Lower, Q = H(1) H(2) . . . H(nq-1).

    Arguments
    ---------
    @param[in]
    side    magma_tally4_side_t
      -     = Magma_tally4Left:      apply Q or Q**H from the Left;
      -     = Magma_tally4Right:     apply Q or Q**H from the Right.

    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper: Upper triangle of A contains elementary reflectors
                   from ZHETRD;
      -     = Magma_tally4Lower: Lower triangle of A contains elementary reflectors
                   from ZHETRD.

    @param[in]
    trans   magma_tally4_trans_t
      -     = Magma_tally4NoTrans:    No transpose, apply Q;
      -     = Magma_tally4_ConjTrans: Conjugate transpose, apply Q**H.

    @param[in]
    m       INTEGER
            The number of rows of the matrix C. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C. N >= 0.

    @param[in]
    dA      COMPLEX_16 array, dimension
                                 (LDDA,M) if SIDE = Magma_tally4Left
                                 (LDDA,N) if SIDE = Magma_tally4Right
            The vectors which define the elementary reflectors, as
            returned by ZHETRD_GPU. On output the diagonal, the subdiagonal and the
            upper part (UPLO=Magma_tally4Lower) or lower part (UPLO=Magma_tally4Upper) are destroyed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            LDDA >= max(1,M) if SIDE = Magma_tally4Left; LDDA >= max(1,N) if SIDE = Magma_tally4Right.

    @param[in]
    tau     COMPLEX_16 array, dimension
                                 (M-1) if SIDE = Magma_tally4Left
                                 (N-1) if SIDE = Magma_tally4Right
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by ZHETRD.

    @param[in,out]
    dC      COMPLEX_16 array, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by (Q*C) or (Q**H * C) or (C * Q**H) or (C*Q).

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDDC >= max(1,M).

    @param[in]
    wA      (workspace) COMPLEX_16 array, dimension
                                 (LDWA,M) if SIDE = Magma_tally4Left
                                 (LDWA,N) if SIDE = Magma_tally4Right
            The vectors which define the elementary reflectors, as
            returned by ZHETRD_GPU.

    @param[in]
    ldwa    INTEGER
            The leading dimension of the array wA.
            LDWA >= max(1,M) if SIDE = Magma_tally4Left; LDWA >= max(1,N) if SIDE = Magma_tally4Right.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_zheev_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zunmtr_gpu(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex   *tau,
    magma_tally4DoubleComplex_ptr dC,    magma_tally4_int_t lddc,
    magma_tally4DoubleComplex    *wA,    magma_tally4_int_t ldwa,
    magma_tally4_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    #define wA(i_,j_) (wA + (i_) + (j_)*ldwa)
    
    magma_tally4_int_t i1, i2, mi, ni, nq;
    int left, upper;
    magma_tally4_int_t iinfo;

    *info = 0;
    left   = (side == Magma_tally4Left);
    upper  = (uplo == Magma_tally4Upper);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        //nw = n;
    } else {
        nq = n;
        //nw = m;
    }
    if (! left && side != Magma_tally4Right) {
        *info = -1;
    } else if (! upper && uplo != Magma_tally4Lower) {
        *info = -2;
    } else if (trans != Magma_tally4NoTrans &&
               trans != Magma_tally4_ConjTrans) {
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
        magma_tally4_xerbla( __func__, -(*info) );
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
        magma_tally4_zunmql2_gpu(side, trans, mi, ni, nq-1, dA(0,1), ldda, tau,
                          dC, lddc, wA(0,1), ldwa, &iinfo);
    }
    else {
        /* Q was determined by a call to ZHETRD with UPLO = 'L' */
        if (left) {
            i1 = 1;
            i2 = 0;
        } else {
            i1 = 0;
            i2 = 1;
        }
        magma_tally4_zunmqr2_gpu(side, trans, mi, ni, nq-1, dA(1,0), ldda, tau,
                          dC(i1,i2), lddc, wA(1,0), ldwa, &iinfo);
    }

    return *info;
} /* magma_tally4_zunmtr */