/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Mark Gates
       @generated from zlarfb_gpu.cpp normal z -> d, Fri Jan 30 19:00:15 2015
*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    DLARFB applies a real block reflector H or its transpose H^H to a
    DOUBLE_PRECISION m by n matrix C, from the left.

    Arguments
    ---------
    @param[in]
    side    magma_tally3_side_t
      -     = Magma_tally3Left:      apply H or H^H from the Left
      -     = Magma_tally3Right:     apply H or H^H from the Right

    @param[in]
    trans   magma_tally3_trans_t
      -     = Magma_tally3NoTrans:    apply H   (No transpose)
      -     = Magma_tally3Trans: apply H^H (Conjugate transpose)

    @param[in]
    direct  magma_tally3_direct_t
            Indicates how H is formed from a product of elementary
            reflectors
      -     = Magma_tally3Forward:  H = H(1) H(2) . . . H(k) (Forward)
      -     = Magma_tally3Backward: H = H(k) . . . H(2) H(1) (Backward)

    @param[in]
    storev  magma_tally3_storev_t
            Indicates how the vectors which define the elementary
            reflectors are stored:
      -     = Magma_tally3Columnwise: Columnwise
      -     = Magma_tally3Rowwise:    Rowwise

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
                (LDDV,K) if STOREV = Magma_tally3Columnwise
                (LDDV,M) if STOREV = Magma_tally3Rowwise and SIDE = Magma_tally3Left
                (LDDV,N) if STOREV = Magma_tally3Rowwise and SIDE = Magma_tally3Right
            The matrix V. See further details.

    @param[in]
    lddv    INTEGER
            The leading dimension of the array V.
            If STOREV = Magma_tally3Columnwise and SIDE = Magma_tally3Left, LDDV >= max(1,M);
            if STOREV = Magma_tally3Columnwise and SIDE = Magma_tally3Right, LDDV >= max(1,N);
            if STOREV = Magma_tally3Rowwise, LDDV >= K.

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
            If SIDE = Magma_tally3Left,  LDWORK >= max(1,N);
            if SIDE = Magma_tally3Right, LDWORK >= max(1,M);

    Further Details
    ---------------
    The shape of the matrix V and the storage of the vectors which define
    the H(i) is best illustrated by the following example with n = 5 and
    k = 3.
    All elements including 0's and 1's are stored, unlike LAPACK.

        DIRECT = Magma_tally3Forward and         DIRECT = Magma_tally3Forward and
        STOREV = Magma_tally3Columnwise:         STOREV = Magma_tally3Rowwise:

                 V = (  1  0  0 )                 V = (  1 v1 v1 v1 v1 )
                     ( v1  1  0 )                     (  0  1 v2 v2 v2 )
                     ( v1 v2  1 )                     (  0  0  1 v3 v3 )
                     ( v1 v2 v3 )
                     ( v1 v2 v3 )

        DIRECT = Magma_tally3Backward and        DIRECT = Magma_tally3Backward and
        STOREV = Magma_tally3Columnwise:         STOREV = Magma_tally3Rowwise:

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1  0  0 )
                     ( v1 v2 v3 )                     ( v2 v2 v2  1  0 )
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
                     (  0  1 v3 )
                     (  0  0  1 )

    @ingroup magma_tally3_daux3
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_dlarfb_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_const_ptr dV,    magma_tally3_int_t lddv,
    magma_tally3Double_const_ptr dT,    magma_tally3_int_t lddt,
    magma_tally3Double_ptr dC,          magma_tally3_int_t lddc,
    magma_tally3Double_ptr dwork,       magma_tally3_int_t ldwork )
{
    #define dV(i_,j_)  (dV    + (i_) + (j_)*lddv)
    #define dT(i_,j_)  (dT    + (i_) + (j_)*lddt)
    #define dC(i_,j_)  (dC    + (i_) + (j_)*lddc)
    #define dwork(i_)  (dwork + (i_))
    
    double c_zero    = MAGMA_tally3_D_ZERO;
    double c_one     = MAGMA_tally3_D_ONE;
    double c_neg_one = MAGMA_tally3_D_NEG_ONE;

    /* Check input arguments */
    magma_tally3_int_t info = 0;
    if (m < 0) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (k < 0) {
        info = -7;
    } else if ( ((storev == Magma_tally3Columnwise) && (side == Magma_tally3Left) && lddv < max(1,m)) ||
                ((storev == Magma_tally3Columnwise) && (side == Magma_tally3Right) && lddv < max(1,n)) ||
                ((storev == Magma_tally3Rowwise) && lddv < k) ) {
        info = -9;
    } else if (lddt < k) {
        info = -11;
    } else if (lddc < max(1,m)) {
        info = -13;
    } else if ( ((side == Magma_tally3Left) && ldwork < max(1,n)) ||
                ((side == Magma_tally3Right) && ldwork < max(1,m)) ) {
        info = -15;
    }
    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return info;
    }
    
    /* Function Body */
    if (m <= 0 || n <= 0) {
        return info;
    }

    // opposite of trans
    magma_tally3_trans_t transt;
    if (trans == Magma_tally3NoTrans)
        transt = Magma_tally3Trans;
    else
        transt = Magma_tally3NoTrans;
    
    // whether T is upper or lower triangular
    magma_tally3_uplo_t uplo;
    if (direct == Magma_tally3Forward)
        uplo = Magma_tally3Upper;
    else
        uplo = Magma_tally3Lower;
    
    // whether V is stored transposed or not
    magma_tally3_trans_t notransV, transV;
    if (storev == Magma_tally3Columnwise) {
        notransV = Magma_tally3NoTrans;
        transV   = Magma_tally3Trans;
    }
    else {
        notransV = Magma_tally3Trans;
        transV   = Magma_tally3NoTrans;
    }

    if ( side == Magma_tally3Left ) {
        // Form H C or H^H C
        // Comments assume H C. When forming H^H C, T gets transposed via transt.
        
        // W = C^H V
        magma_tally3_dgemm( Magma_tally3Trans, notransV,
                     n, k, m,
                     c_one,  dC(0,0),  lddc,
                             dV(0,0),  lddv,
                     c_zero, dwork(0), ldwork);

        // W = W T^H = C^H V T^H
        magma_tally3_dtrmm( Magma_tally3Right, uplo, transt, Magma_tally3NonUnit,
                     n, k,
                     c_one, dT(0,0),  lddt,
                            dwork(0), ldwork);

        // C = C - V W^H = C - V T V^H C = (I - V T V^H) C = H C
        magma_tally3_dgemm( notransV, Magma_tally3Trans,
                     m, n, k,
                     c_neg_one, dV(0,0),  lddv,
                                dwork(0), ldwork,
                     c_one,     dC(0,0),  lddc);
    }
    else {
        // Form C H or C H^H
        // Comments assume C H. When forming C H^H, T gets transposed via trans.
        
        // W = C V
        magma_tally3_dgemm( Magma_tally3NoTrans, notransV,
                     m, k, n,
                     c_one,  dC(0,0),  lddc,
                             dV(0,0),  lddv,
                     c_zero, dwork(0), ldwork);

        // W = W T = C V T
        magma_tally3_dtrmm( Magma_tally3Right, uplo, trans, Magma_tally3NonUnit,
                     m, k,
                     c_one, dT(0,0),  lddt,
                            dwork(0), ldwork);

        // C = C - W V^H = C - C V T V^H = C (I - V T V^H) = C H
        magma_tally3_dgemm( Magma_tally3NoTrans, transV,
                     m, n, k,
                     c_neg_one, dwork(0), ldwork,
                                dV(0,0),  lddv,
                     c_one,     dC(0,0),  lddc);
    }

    return info;
} /* magma_tally3_dlarfb */
