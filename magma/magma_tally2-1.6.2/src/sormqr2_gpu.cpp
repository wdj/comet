/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov

       @generated from zunmqr2_gpu.cpp normal z -> s, Fri Jan 30 19:00:16 2015

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    SORMQR overwrites the general real M-by-N matrix C with

    @verbatim
                               SIDE = Magma_tally2Left    SIDE = Magma_tally2Right
    TRANS = Magma_tally2NoTrans:      Q * C               C * Q
    TRANS = Magma_tally2Trans:   Q**H * C            C * Q**H
    @endverbatim

    where Q is a real unitary matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by SGEQRF. Q is of order M if SIDE = Magma_tally2Left and of order N
    if SIDE = Magma_tally2Right.

    Arguments
    ---------
    @param[in]
    side    magma_tally2_side_t
      -     = Magma_tally2Left:      apply Q or Q**H from the Left;
      -     = Magma_tally2Right:     apply Q or Q**H from the Right.

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
    k       INTEGER
            The number of elementary reflectors whose product defines
            the matrix Q.
            If SIDE = Magma_tally2Left,  M >= K >= 0;
            if SIDE = Magma_tally2Right, N >= K >= 0.

    @param[in]
    dA      REAL array, dimension (LDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            SGEQRF in the first k columns of its array argument A.
            The diagonal and the upper part
            are destroyed, the reflectors are not modified.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            LDDA >= max(1,M) if SIDE = Magma_tally2Left; LDDA >= max(1,N) if SIDE = Magma_tally2Right.

    @param[in]
    tau     REAL array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SGEQRF.

    @param[in,out]
    dC      REAL array, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by (Q*C) or (Q**H * C) or (C * Q**H) or (C*Q).

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDDC >= max(1,M).

    @param[in]
    wA      (workspace) REAL array, dimension
                                 (LDWA,M) if SIDE = Magma_tally2Left
                                 (LDWA,N) if SIDE = Magma_tally2Right
            The vectors which define the elementary reflectors, as
            returned by SSYTRD_GPU.

    @param[in]
    ldwa    INTEGER
            The leading dimension of the array wA.
            LDWA >= max(1,M) if SIDE = Magma_tally2Left; LDWA >= max(1,N) if SIDE = Magma_tally2Right.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_sgeqrf_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_sormqr2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float    *tau,
    magma_tally2Float_ptr dC, magma_tally2_int_t lddc,
    float    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    #define wA(i_,j_) (wA + (i_) + (j_)*ldwa)
    
    /* Allocate work space on the GPU */
    magma_tally2Float_ptr dwork;

    float c_zero = MAGMA_tally2_S_ZERO;
    float c_one  = MAGMA_tally2_S_ONE;
    
    magma_tally2_int_t i, i__4, lddwork;
    float T[2*4160]        /* was [65][64] */;
    magma_tally2_int_t i1, i2, step, ib, ic, jc, nb, mi, ni, nq;
    int left, notran;

    wA -= 1 + ldwa;
    dC -= 1 + lddc;
    --tau;

    *info = 0;
    left   = (side == Magma_tally2Left);
    notran = (trans == Magma_tally2NoTrans);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        //nw = n;
        magma_tally2_smalloc( &dwork, (n + 64)*64 );  // TODO after checking args, else memory leak!
    } else {
        nq = n;
        //nw = m;
        magma_tally2_smalloc( &dwork, (m + 64)*64 );  // TODO after checking args, else memory leak!
    }
    if (! left && side != Magma_tally2Right) {
        *info = -1;
    } else if (! notran && trans != Magma_tally2Trans) {
        *info = -2;
    } else if (m < 0) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (k < 0 || k > nq) {
        *info = -5;
    } else if (ldda < max(1,nq)) {
        *info = -7;
    } else if (lddc < max(1,m)) {
        *info = -10;
    } else if (ldwa < max(1,nq)) {
        *info = -12;
    }

    // size of the block
    nb = 64;

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
        return *info;
    }

    /* Use hybrid CPU-GPU code */
    if ( ( left && (! notran) ) ||  ( (! left) && notran ) ) {
        i1 = 1;
        i2 = k;
        step = nb;
    } else {
        i1 = ((k - 1)/nb)*nb + 1;
        i2 = 1;
        step = -nb;
    }

    // silence "uninitialized" warnings
    mi = 0;
    ni = 0;
    
    if (left) {
        ni = n;
        jc = 1;
    } else {
        mi = m;
        ic = 1;
    }

    // set nb-1 super-diagonals to 0, and diagonal to 1.
    // This way we can copy V directly to the GPU,
    // with the upper triangle parts already set to identity.
    magma_tally2blas_slaset_band( Magma_tally2Upper, k, k, nb, c_zero, c_one, dA, ldda );

    // for i=i1 to i2 by step
    for (i = i1; (step < 0 ? i >= i2 : i <= i2); i += step) {
        ib = min(nb, k - i + 1);

        /* Form the triangular factor of the block reflector
           H = H(i) H(i+1) . . . H(i+ib-1) */
        i__4 = nq - i + 1;
        lapackf77_slarft("Forward", "Columnwise", &i__4, &ib,
                         wA(i,i), &ldwa, &tau[i], T, &ib);

        if (left) {
            /* H or H' is applied to C(i:m,1:n) */
            mi = m - i + 1;
            ic = i;
        }
        else {
            /* H or H' is applied to C(1:m,i:n) */
            ni = n - i + 1;
            jc = i;
        }

        if (left)
            lddwork = ni;
        else
            lddwork = mi;

        /* Apply H or H'; First copy T to the GPU */
        magma_tally2_ssetmatrix( ib, ib, T, ib, dwork, ib );
        magma_tally2_slarfb_gpu( side, trans, Magma_tally2Forward, Magma_tally2Columnwise,
                          mi, ni, ib,
                          dA(i-1,i-1), ldda, dwork, ib,  // dA using 0-based indices here
                          dC(ic,jc), lddc,
                          dwork + ib*ib, lddwork);
    }

    magma_tally2_free( dwork );

    return *info;
} /* magma_tally2_sormqr */
