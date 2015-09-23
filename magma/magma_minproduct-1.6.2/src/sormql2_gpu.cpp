/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca

       @generated from zunmql2_gpu.cpp normal z -> s, Fri Jan 30 19:00:15 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    SORMQL overwrites the general real M-by-N matrix C with

    @verbatim
                               SIDE = Magma_minproductLeft   SIDE = Magma_minproductRight
    TRANS = Magma_minproductNoTrans:      Q * C              C * Q
    TRANS = Magma_minproductTrans:   Q**H * C           C * Q**H
    @endverbatim

    where Q is a real unitary matrix defined as the product of k
    elementary reflectors

          Q = H(k) . . . H(2) H(1)

    as returned by SGEQLF. Q is of order M if SIDE = Magma_minproductLeft and of order N
    if SIDE = Magma_minproductRight.

    Arguments
    ---------
    @param[in]
    side    magma_minproduct_side_t
      -     = Magma_minproductLeft:      apply Q or Q**H from the Left;
      -     = Magma_minproductRight:     apply Q or Q**H from the Right.

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
    k       INTEGER
            The number of elementary reflectors whose product defines
            the matrix Q.
            If SIDE = Magma_minproductLeft,  M >= K >= 0;
            if SIDE = Magma_minproductRight, N >= K >= 0.

    @param[in]
    dA      REAL array, dimension (LDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            SGEQLF in the last k columns of its array argument A.
            The diagonal and the lower part
            are destroyed, the reflectors are not modified.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            LDDA >= max(1,M) if SIDE = Magma_minproductLeft;
            LDDA >= max(1,N) if SIDE = Magma_minproductRight.

    @param[in]
    tau     REAL array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SGEQLF.

    @param[in,out]
    dC      REAL array, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

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

    @ingroup magma_minproduct_sgeqlf_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_sormql2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float    *tau,
    magma_minproductFloat_ptr dC, magma_minproduct_int_t lddc,
    float    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    #define wA(i_,j_) (wA + (i_) + (j_)*ldwa)
    
    /* Allocate work space on the GPU */
    magma_minproductFloat_ptr dwork;
    magma_minproduct_smalloc( &dwork, 2*(m + 64)*64 );

    float c_zero = MAGMA_minproduct_S_ZERO;
    float c_one  = MAGMA_minproduct_S_ONE;
    
    magma_minproduct_int_t i, i__4;
    float T[2*4160]        /* was [65][64] */;
    magma_minproduct_int_t i1, i2, step, ib, nb, mi, ni, nq, nw;
    magma_minproduct_int_t ldwork;
    int left, notran;

    wA -= 1 + ldwa;
    dC -= 1 + lddc;
    --tau;

    *info  = 0;
    left   = (side == Magma_minproductLeft);
    notran = (trans == Magma_minproductNoTrans);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = max(1,n);
    } else {
        nq = n;
        nw = max(1,m);
    }
    if (! left && side != Magma_minproductRight) {
        *info = -1;
    } else if (! notran && trans != Magma_minproductTrans) {
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
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0) {
        return *info;
    }

    ldwork = nw;
        
    /* Use hybrid CPU-GPU code */
    if ((left && notran) || (! left && ! notran)) {
        i1 = 1;
        i2 = k;
        step = nb;
    } else {
        i1 = (k - 1) / nb * nb + 1;
        i2 = 1;
        step = -nb;
    }
    
    // silence "uninitialized" warnings
    mi = 0;
    ni = 0;
    
    if (left) {
        ni = n;
    } else {
        mi = m;
    }
    
    // set nb-1 sub-diagonals to 0, and diagonal to 1.
    // This way we can copy V directly to the GPU,
    // already with the lower triangle parts already set to identity.
    magma_minproductblas_slaset_band( Magma_minproductLower, k, k, nb, c_zero, c_one, dA, ldda );
    
    for (i = i1; (step < 0 ? i >= i2 : i <= i2); i += step) {
        ib = min(nb, k - i + 1);
        
        /* Form the triangular factor of the block reflector
           H = H(i+ib-1) . . . H(i+1) H(i) */
        i__4 = nq - k + i + ib - 1;
        lapackf77_slarft("Backward", "Columnwise", &i__4, &ib,
                         wA(1,i), &ldwa, &tau[i], T, &ib);
    
        if (left) {
            /* H or H' is applied to C(1:m-k+i+ib-1,1:n) */
            mi = m - k + i + ib - 1;
        }
        else {
            /* H or H' is applied to C(1:m,1:n-k+i+ib-1) */
            ni = n - k + i + ib - 1;
        }
        
        /* Apply H or H'; First copy T to the GPU */
        magma_minproduct_ssetmatrix( ib, ib, T, ib, dwork+i__4*ib, ib );
        magma_minproduct_slarfb_gpu(side, trans, Magma_minproductBackward, Magma_minproductColumnwise,
                         mi, ni, ib,
                         dA(0,i-1), ldda, dwork+i__4*ib, ib,  // dA using 0-based indices here
                         dC(1,1), lddc,
                         dwork+i__4*ib + ib*ib, ldwork);
    }

    magma_minproduct_free( dwork );

    return *info;
} /* magma_minproduct_sormql */
