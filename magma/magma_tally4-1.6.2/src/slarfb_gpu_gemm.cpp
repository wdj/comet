/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @author Azzam Haidar
       @generated from zlarfb_gpu_gemm.cpp normal z -> s, Fri Jan 30 19:00:16 2015
*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    SLARFB applies a real block reflector H or its transpose H^H to a
    REAL m by n matrix C, from the left.
    
    __Note that this function assumes__ that the upper part of dV is 0
    because it is referenced. Same for upper/lower part of dT.

    Arguments
    ---------
    @param[in]
    side    magma_tally4_side_t
      -     = Magma_tally4Left:      apply H or H^H from the Left
      -     = Magma_tally4Right:     apply H or H^H from the Right

    @param[in]
    trans   magma_tally4_trans_t
      -     = Magma_tally4NoTrans:    apply H   (No transpose)
      -     = Magma_tally4Trans: apply H^H (Conjugate transpose)

    @param[in]
    direct  magma_tally4_direct_t
            Indicates how H is formed from a product of elementary
            reflectors
      -     = Magma_tally4Forward:  H = H(1) H(2) . . . H(k) (Forward)
      -     = Magma_tally4Backward: H = H(k) . . . H(2) H(1) (Backward)

    @param[in]
    storev  magma_tally4_storev_t
            Indicates how the vectors which define the elementary
            reflectors are stored:
      -     = Magma_tally4Columnwise: Columnwise
      -     = Magma_tally4Rowwise:    Rowwise

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
    dV      REAL array on the GPU, dimension
                (LDDV,K) if STOREV = Magma_tally4Columnwise
                (LDDV,M) if STOREV = Magma_tally4Rowwise and SIDE = Magma_tally4Left
                (LDDV,N) if STOREV = Magma_tally4Rowwise and SIDE = Magma_tally4Right
            The matrix V. See further details.

    @param[in]
    lddv    INTEGER
            The leading dimension of the array V.
            If STOREV = Magma_tally4Columnwise and SIDE = Magma_tally4Left, LDDV >= max(1,M);
            if STOREV = Magma_tally4Columnwise and SIDE = Magma_tally4Right, LDDV >= max(1,N);
            if STOREV = Magma_tally4Rowwise, LDDV >= K.

    @param[in]
    dT      REAL array on the GPU, dimension (LDDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    @param[in]
    lddt    INTEGER
            The leading dimension of the array T. LDDT >= K.

    @param[in,out]
    dC      REAL array on the GPU, dimension (LDDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C, or H^H*C, or C*H, or C*H^H.

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDA >= max(1,M).

    @param
    dwork   (workspace) REAL array, dimension (LDWORK,K)

    @param[in]
    ldwork  INTEGER
            The leading dimension of the array WORK.
            If SIDE = Magma_tally4Left,  LDWORK >= max(1,N);
            if SIDE = Magma_tally4Right, LDWORK >= max(1,M);

    @param
    dworkvt (workspace) REAL array, dimension (LDWORKT,K)

    @param[in]
    ldworkvt INTEGER
            The leading dimension of the array WORKVT.
            LDWORKVT >= max(1,min(M,N));

    Further Details
    ---------------
    The shape of the matrix V and the storage of the vectors which define
    the H(i) is best illustrated by the following example with n = 5 and
    k = 3.
    All elements including 0's and 1's are stored, unlike LAPACK.

        DIRECT = Magma_tally4Forward and         DIRECT = Magma_tally4Forward and
        STOREV = Magma_tally4Columnwise:         STOREV = Magma_tally4Rowwise:

                 V = (  1  0  0 )                 V = (  1 v1 v1 v1 v1 )
                     ( v1  1  0 )                     (  0  1 v2 v2 v2 )
                     ( v1 v2  1 )                     (  0  0  1 v3 v3 )
                     ( v1 v2 v3 )
                     ( v1 v2 v3 )

        DIRECT = Magma_tally4Backward and        DIRECT = Magma_tally4Backward and
        STOREV = Magma_tally4Columnwise:         STOREV = Magma_tally4Rowwise:

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1  0  0 )
                     ( v1 v2 v3 )                     ( v2 v2 v2  1  0 )
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
                     (  0  1 v3 )
                     (  0  0  1 )

    @ingroup magma_tally4_saux3
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_slarfb_gpu_gemm(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_const_ptr dV,    magma_tally4_int_t lddv,
    magma_tally4Float_const_ptr dT,    magma_tally4_int_t lddt,
    magma_tally4Float_ptr dC,          magma_tally4_int_t lddc,
    magma_tally4Float_ptr dwork,       magma_tally4_int_t ldwork,
    magma_tally4Float_ptr dworkvt,     magma_tally4_int_t ldworkvt)
{
    float c_zero    = MAGMA_tally4_S_ZERO;
    float c_one     = MAGMA_tally4_S_ONE;
    float c_neg_one = MAGMA_tally4_S_NEG_ONE;

    magma_tally4_int_t info = 0;
    
    /* Function Body */
    if (m <= 0 || n <= 0) {
        return info;
    }
    
    // internal variable
    magma_tally4_int_t ldwvt = (m > n ?  k : m);
    magma_tally4_int_t ldw;
    if ( side == Magma_tally4Left ) {
        ldw = k;
    } else {
        ldw = m;
    }
    
    // opposite of trans
    magma_tally4_trans_t transt;
    if (trans == Magma_tally4NoTrans)
        transt = Magma_tally4Trans;
    else
        transt = Magma_tally4NoTrans;
    
    MAGMA_tally4_UNUSED( transt );  // TODO: is this a bug that it isn't used?
    
    // whether V is stored transposed or not
    magma_tally4_trans_t notransV, transV;
    if (storev == Magma_tally4Columnwise) {
        notransV = Magma_tally4NoTrans;
        transV   = Magma_tally4Trans;
    }
    else {
        notransV = Magma_tally4Trans;
        transV   = Magma_tally4NoTrans;
    }

    if ( side == Magma_tally4Left ) {
        // Form H C or H^H C
        // Comments assume H C.
        // When forming H^H C, T gets transposed via transt for m >= n or by trans for m < n.
        
        // W = V^H C
        magma_tally4_sgemm( Magma_tally4Trans, notransV,
                     k, n, m,
                     c_one,  dV,    lddv,
                             dC,    lddc,
                     c_zero, dwork, ldw);

        if (m <= n) {
            // W2 = V T
            magma_tally4_sgemm( notransV, trans,
                         m, k, k,
                         c_one,  dV, lddv,
                                 dT, lddt,
                         c_zero, dworkvt, ldwvt);
            // C = C - W2 W = C - V T V^H C = (I - V T V^H) C = H C
            magma_tally4_sgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                         m, n, k,
                         c_neg_one, dworkvt,  ldwvt,
                                    dwork,    ldw,
                         c_one,     dC,       lddc);
        } else {
            // W2 = T W  = T  V^H C
            magma_tally4_sgemm( trans, Magma_tally4NoTrans,
                         k, n, k,
                         c_one,  dT, lddt,
                                 dwork, ldw,
                         c_zero, dworkvt, ldwvt);
            // C = C - V W2 = C - V T V^H C = (I - V T V^H) C = H C
            magma_tally4_sgemm( notransV, Magma_tally4NoTrans,
                         m, n, k,
                         c_neg_one, dV,  lddv,
                                    dworkvt,  ldwvt,
                         c_one,     dC,       lddc);
        }
    }
    else {
        // Form C H or C H^H
        // Comments assume C H.
        // When forming C H^H, T gets transposed via trans.
        
        // W = C V
        magma_tally4_sgemm( Magma_tally4NoTrans, notransV,
                     m, k, n,
                     c_one,  dC,    lddc,
                             dV,    lddv,
                     c_zero, dwork, ldw);
        if (m <= n) {
            // W2 = W T = C V T
            magma_tally4_sgemm( Magma_tally4NoTrans, trans,
                         m, k, k,
                         c_one,  dwork, ldw,
                                 dT, lddt,
                         c_zero, dworkvt, ldwvt);
            // C = C - W2 V^H = C - C V T V^H = C (I - V T V^H) = C H
            magma_tally4_sgemm( Magma_tally4NoTrans, transV,
                         m, n, k,
                         c_neg_one, dworkvt, ldwvt,
                                    dV,    lddv,
                         c_one,     dC,    lddc);
        } else {
            // W2 = T V^H
            magma_tally4_sgemm( trans, transV,
                         k, n, k,
                         c_one,  dT, lddt,
                                 dV, lddv,
                         c_zero, dworkvt, ldwvt);
            // C = C - W W2 = C - C V T V^H = C (I - V T V^H) = C H
            magma_tally4_sgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                         m, n, k,
                         c_neg_one, dwork,   ldw,
                                    dworkvt, ldwvt,
                         c_one,     dC,      lddc);
        }
    }

    return info;
} /* magma_tally4_slarfb */
