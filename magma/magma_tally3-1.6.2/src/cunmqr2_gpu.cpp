/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov

       @generated from zunmqr2_gpu.cpp normal z -> c, Fri Jan 30 19:00:16 2015

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    CUNMQR overwrites the general complex M-by-N matrix C with

    @verbatim
                               SIDE = Magma_tally3Left    SIDE = Magma_tally3Right
    TRANS = Magma_tally3NoTrans:      Q * C               C * Q
    TRANS = Magma_tally3_ConjTrans:   Q**H * C            C * Q**H
    @endverbatim

    where Q is a complex unitary matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by CGEQRF. Q is of order M if SIDE = Magma_tally3Left and of order N
    if SIDE = Magma_tally3Right.

    Arguments
    ---------
    @param[in]
    side    magma_tally3_side_t
      -     = Magma_tally3Left:      apply Q or Q**H from the Left;
      -     = Magma_tally3Right:     apply Q or Q**H from the Right.

    @param[in]
    trans   magma_tally3_trans_t
      -     = Magma_tally3NoTrans:    No transpose, apply Q;
      -     = Magma_tally3_ConjTrans: Conjugate transpose, apply Q**H.

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
            If SIDE = Magma_tally3Left,  M >= K >= 0;
            if SIDE = Magma_tally3Right, N >= K >= 0.

    @param[in]
    dA      COMPLEX array, dimension (LDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            CGEQRF in the first k columns of its array argument A.
            The diagonal and the upper part
            are destroyed, the reflectors are not modified.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            LDDA >= max(1,M) if SIDE = Magma_tally3Left; LDDA >= max(1,N) if SIDE = Magma_tally3Right.

    @param[in]
    tau     COMPLEX array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by CGEQRF.

    @param[in,out]
    dC      COMPLEX array, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by (Q*C) or (Q**H * C) or (C * Q**H) or (C*Q).

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDDC >= max(1,M).

    @param[in]
    wA      (workspace) COMPLEX array, dimension
                                 (LDWA,M) if SIDE = Magma_tally3Left
                                 (LDWA,N) if SIDE = Magma_tally3Right
            The vectors which define the elementary reflectors, as
            returned by CHETRD_GPU.

    @param[in]
    ldwa    INTEGER
            The leading dimension of the array wA.
            LDWA >= max(1,M) if SIDE = Magma_tally3Left; LDWA >= max(1,N) if SIDE = Magma_tally3Right.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_cgeqrf_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_cunmqr2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex    *tau,
    magma_tally3FloatComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3FloatComplex    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    #define wA(i_,j_) (wA + (i_) + (j_)*ldwa)
    
    /* Allocate work space on the GPU */
    magma_tally3FloatComplex_ptr dwork;

    magma_tally3FloatComplex c_zero = MAGMA_tally3_C_ZERO;
    magma_tally3FloatComplex c_one  = MAGMA_tally3_C_ONE;
    
    magma_tally3_int_t i, i__4, lddwork;
    magma_tally3FloatComplex T[2*4160]        /* was [65][64] */;
    magma_tally3_int_t i1, i2, step, ib, ic, jc, nb, mi, ni, nq;
    int left, notran;

    wA -= 1 + ldwa;
    dC -= 1 + lddc;
    --tau;

    *info = 0;
    left   = (side == Magma_tally3Left);
    notran = (trans == Magma_tally3NoTrans);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        //nw = n;
        magma_tally3_cmalloc( &dwork, (n + 64)*64 );  // TODO after checking args, else memory leak!
    } else {
        nq = n;
        //nw = m;
        magma_tally3_cmalloc( &dwork, (m + 64)*64 );  // TODO after checking args, else memory leak!
    }
    if (! left && side != Magma_tally3Right) {
        *info = -1;
    } else if (! notran && trans != Magma_tally3_ConjTrans) {
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
        magma_tally3_xerbla( __func__, -(*info) );
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
    magma_tally3blas_claset_band( Magma_tally3Upper, k, k, nb, c_zero, c_one, dA, ldda );

    // for i=i1 to i2 by step
    for (i = i1; (step < 0 ? i >= i2 : i <= i2); i += step) {
        ib = min(nb, k - i + 1);

        /* Form the triangular factor of the block reflector
           H = H(i) H(i+1) . . . H(i+ib-1) */
        i__4 = nq - i + 1;
        lapackf77_clarft("Forward", "Columnwise", &i__4, &ib,
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
        magma_tally3_csetmatrix( ib, ib, T, ib, dwork, ib );
        magma_tally3_clarfb_gpu( side, trans, Magma_tally3Forward, Magma_tally3Columnwise,
                          mi, ni, ib,
                          dA(i-1,i-1), ldda, dwork, ib,  // dA using 0-based indices here
                          dC(ic,jc), lddc,
                          dwork + ib*ib, lddwork);
    }

    magma_tally3_free( dwork );

    return *info;
} /* magma_tally3_cunmqr */
