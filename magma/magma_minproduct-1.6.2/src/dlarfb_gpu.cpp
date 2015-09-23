/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Mark Gates
       @generated from zlarfb_gpu.cpp normal z -> d, Fri Jan 30 19:00:15 2015
*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    DLARFB applies a real block reflector H or its transpose H^H to a
    DOUBLE_PRECISION m by n matrix C, from the left.

    Arguments
    ---------
    @param[in]
    side    magma_minproduct_side_t
      -     = Magma_minproductLeft:      apply H or H^H from the Left
      -     = Magma_minproductRight:     apply H or H^H from the Right

    @param[in]
    trans   magma_minproduct_trans_t
      -     = Magma_minproductNoTrans:    apply H   (No transpose)
      -     = Magma_minproductTrans: apply H^H (Conjugate transpose)

    @param[in]
    direct  magma_minproduct_direct_t
            Indicates how H is formed from a product of elementary
            reflectors
      -     = Magma_minproductForward:  H = H(1) H(2) . . . H(k) (Forward)
      -     = Magma_minproductBackward: H = H(k) . . . H(2) H(1) (Backward)

    @param[in]
    storev  magma_minproduct_storev_t
            Indicates how the vectors which define the elementary
            reflectors are stored:
      -     = Magma_minproductColumnwise: Columnwise
      -     = Magma_minproductRowwise:    Rowwise

    @param[in]
    m       INTEGER
            The number of rows of the matrix C.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C.

    @param[in]
    k       INTEGER
            The order of the matrix T (= the number of elementary
            reflectors whose product defines the block reflector).

    @param[in]
    dV      DOUBLE_PRECISION array on the GPU, dimension
                (LDDV,K) if STOREV = Magma_minproductColumnwise
                (LDDV,M) if STOREV = Magma_minproductRowwise and SIDE = Magma_minproductLeft
                (LDDV,N) if STOREV = Magma_minproductRowwise and SIDE = Magma_minproductRight
            The matrix V. See further details.

    @param[in]
    lddv    INTEGER
            The leading dimension of the array V.
            If STOREV = Magma_minproductColumnwise and SIDE = Magma_minproductLeft, LDDV >= max(1,M);
            if STOREV = Magma_minproductColumnwise and SIDE = Magma_minproductRight, LDDV >= max(1,N);
            if STOREV = Magma_minproductRowwise, LDDV >= K.

    @param[in]
    dT      DOUBLE_PRECISION array on the GPU, dimension (LDDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    @param[in]
    lddt    INTEGER
            The leading dimension of the array T. LDDT >= K.

    @param[in,out]
    dC      DOUBLE_PRECISION array on the GPU, dimension (LDDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C, or H^H*C, or C*H, or C*H^H.

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDA >= max(1,M).

    @param
    dwork   (workspace) DOUBLE_PRECISION array, dimension (LDWORK,K)

    @param[in]
    ldwork  INTEGER
            The leading dimension of the array WORK.
            If SIDE = Magma_minproductLeft,  LDWORK >= max(1,N);
            if SIDE = Magma_minproductRight, LDWORK >= max(1,M);

    Further Details
    ---------------
    The shape of the matrix V and the storage of the vectors which define
    the H(i) is best illustrated by the following example with n = 5 and
    k = 3.
    All elements including 0's and 1's are stored, unlike LAPACK.

        DIRECT = Magma_minproductForward and         DIRECT = Magma_minproductForward and
        STOREV = Magma_minproductColumnwise:         STOREV = Magma_minproductRowwise:

                 V = (  1  0  0 )                 V = (  1 v1 v1 v1 v1 )
                     ( v1  1  0 )                     (  0  1 v2 v2 v2 )
                     ( v1 v2  1 )                     (  0  0  1 v3 v3 )
                     ( v1 v2 v3 )
                     ( v1 v2 v3 )

        DIRECT = Magma_minproductBackward and        DIRECT = Magma_minproductBackward and
        STOREV = Magma_minproductColumnwise:         STOREV = Magma_minproductRowwise:

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1  0  0 )
                     ( v1 v2 v3 )                     ( v2 v2 v2  1  0 )
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
                     (  0  1 v3 )
                     (  0  0  1 )

    @ingroup magma_minproduct_daux3
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_dlarfb_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_const_ptr dV,    magma_minproduct_int_t lddv,
    magma_minproductDouble_const_ptr dT,    magma_minproduct_int_t lddt,
    magma_minproductDouble_ptr dC,          magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dwork,       magma_minproduct_int_t ldwork )
{
    #define dV(i_,j_)  (dV    + (i_) + (j_)*lddv)
    #define dT(i_,j_)  (dT    + (i_) + (j_)*lddt)
    #define dC(i_,j_)  (dC    + (i_) + (j_)*lddc)
    #define dwork(i_)  (dwork + (i_))
    
    double c_zero    = MAGMA_minproduct_D_ZERO;
    double c_one     = MAGMA_minproduct_D_ONE;
    double c_neg_one = MAGMA_minproduct_D_NEG_ONE;

    /* Check input arguments */
    magma_minproduct_int_t info = 0;
    if (m < 0) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (k < 0) {
        info = -7;
    } else if ( ((storev == Magma_minproductColumnwise) && (side == Magma_minproductLeft) && lddv < max(1,m)) ||
                ((storev == Magma_minproductColumnwise) && (side == Magma_minproductRight) && lddv < max(1,n)) ||
                ((storev == Magma_minproductRowwise) && lddv < k) ) {
        info = -9;
    } else if (lddt < k) {
        info = -11;
    } else if (lddc < max(1,m)) {
        info = -13;
    } else if ( ((side == Magma_minproductLeft) && ldwork < max(1,n)) ||
                ((side == Magma_minproductRight) && ldwork < max(1,m)) ) {
        info = -15;
    }
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return info;
    }
    
    /* Function Body */
    if (m <= 0 || n <= 0) {
        return info;
    }

    // opposite of trans
    magma_minproduct_trans_t transt;
    if (trans == Magma_minproductNoTrans)
        transt = Magma_minproductTrans;
    else
        transt = Magma_minproductNoTrans;
    
    // whether T is upper or lower triangular
    magma_minproduct_uplo_t uplo;
    if (direct == Magma_minproductForward)
        uplo = Magma_minproductUpper;
    else
        uplo = Magma_minproductLower;
    
    // whether V is stored transposed or not
    magma_minproduct_trans_t notransV, transV;
    if (storev == Magma_minproductColumnwise) {
        notransV = Magma_minproductNoTrans;
        transV   = Magma_minproductTrans;
    }
    else {
        notransV = Magma_minproductTrans;
        transV   = Magma_minproductNoTrans;
    }

    if ( side == Magma_minproductLeft ) {
        // Form H C or H^H C
        // Comments assume H C. When forming H^H C, T gets transposed via transt.
        
        // W = C^H V
        magma_minproduct_dgemm( Magma_minproductTrans, notransV,
                     n, k, m,
                     c_one,  dC(0,0),  lddc,
                             dV(0,0),  lddv,
                     c_zero, dwork(0), ldwork);

        // W = W T^H = C^H V T^H
        magma_minproduct_dtrmm( Magma_minproductRight, uplo, transt, Magma_minproductNonUnit,
                     n, k,
                     c_one, dT(0,0),  lddt,
                            dwork(0), ldwork);

        // C = C - V W^H = C - V T V^H C = (I - V T V^H) C = H C
        magma_minproduct_dgemm( notransV, Magma_minproductTrans,
                     m, n, k,
                     c_neg_one, dV(0,0),  lddv,
                                dwork(0), ldwork,
                     c_one,     dC(0,0),  lddc);
    }
    else {
        // Form C H or C H^H
        // Comments assume C H. When forming C H^H, T gets transposed via trans.
        
        // W = C V
        magma_minproduct_dgemm( Magma_minproductNoTrans, notransV,
                     m, k, n,
                     c_one,  dC(0,0),  lddc,
                             dV(0,0),  lddv,
                     c_zero, dwork(0), ldwork);

        // W = W T = C V T
        magma_minproduct_dtrmm( Magma_minproductRight, uplo, trans, Magma_minproductNonUnit,
                     m, k,
                     c_one, dT(0,0),  lddt,
                            dwork(0), ldwork);

        // C = C - W V^H = C - C V T V^H = C (I - V T V^H) = C H
        magma_minproduct_dgemm( Magma_minproductNoTrans, transV,
                     m, n, k,
                     c_neg_one, dwork(0), ldwork,
                                dV(0,0),  lddv,
                     c_one,     dC(0,0),  lddc);
    }

    return info;
} /* magma_minproduct_dlarfb */
