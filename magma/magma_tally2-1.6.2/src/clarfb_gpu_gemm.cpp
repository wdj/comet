/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @author Azzam Haidar
       @generated from zlarfb_gpu_gemm.cpp normal z -> c, Fri Jan 30 19:00:15 2015
*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    CLARFB applies a complex block reflector H or its transpose H^H to a
    COMPLEX m by n matrix C, from the left.
    
    __Note that this function assumes__ that the upper part of dV is 0
    because it is referenced. Same for upper/lower part of dT.

    Arguments
    ---------
    @param[in]
    side    magma_tally2_side_t
      -     = Magma_tally2Left:      apply H or H^H from the Left
      -     = Magma_tally2Right:     apply H or H^H from the Right

    @param[in]
    trans   magma_tally2_trans_t
      -     = Magma_tally2NoTrans:    apply H   (No transpose)
      -     = Magma_tally2_ConjTrans: apply H^H (Conjugate transpose)

    @param[in]
    direct  magma_tally2_direct_t
            Indicates how H is formed from a product of elementary
            reflectors
      -     = Magma_tally2Forward:  H = H(1) H(2) . . . H(k) (Forward)
      -     = Magma_tally2Backward: H = H(k) . . . H(2) H(1) (Backward)

    @param[in]
    storev  magma_tally2_storev_t
            Indicates how the vectors which define the elementary
            reflectors are stored:
      -     = Magma_tally2Columnwise: Columnwise
      -     = Magma_tally2Rowwise:    Rowwise

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
    dV      COMPLEX array on the GPU, dimension
                (LDDV,K) if STOREV = Magma_tally2Columnwise
                (LDDV,M) if STOREV = Magma_tally2Rowwise and SIDE = Magma_tally2Left
                (LDDV,N) if STOREV = Magma_tally2Rowwise and SIDE = Magma_tally2Right
            The matrix V. See further details.

    @param[in]
    lddv    INTEGER
            The leading dimension of the array V.
            If STOREV = Magma_tally2Columnwise and SIDE = Magma_tally2Left, LDDV >= max(1,M);
            if STOREV = Magma_tally2Columnwise and SIDE = Magma_tally2Right, LDDV >= max(1,N);
            if STOREV = Magma_tally2Rowwise, LDDV >= K.

    @param[in]
    dT      COMPLEX array on the GPU, dimension (LDDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    @param[in]
    lddt    INTEGER
            The leading dimension of the array T. LDDT >= K.

    @param[in,out]
    dC      COMPLEX array on the GPU, dimension (LDDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C, or H^H*C, or C*H, or C*H^H.

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDA >= max(1,M).

    @param
    dwork   (workspace) COMPLEX array, dimension (LDWORK,K)

    @param[in]
    ldwork  INTEGER
            The leading dimension of the array WORK.
            If SIDE = Magma_tally2Left,  LDWORK >= max(1,N);
            if SIDE = Magma_tally2Right, LDWORK >= max(1,M);

    @param
    dworkvt (workspace) COMPLEX array, dimension (LDWORKT,K)

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

        DIRECT = Magma_tally2Forward and         DIRECT = Magma_tally2Forward and
        STOREV = Magma_tally2Columnwise:         STOREV = Magma_tally2Rowwise:

                 V = (  1  0  0 )                 V = (  1 v1 v1 v1 v1 )
                     ( v1  1  0 )                     (  0  1 v2 v2 v2 )
                     ( v1 v2  1 )                     (  0  0  1 v3 v3 )
                     ( v1 v2 v3 )
                     ( v1 v2 v3 )

        DIRECT = Magma_tally2Backward and        DIRECT = Magma_tally2Backward and
        STOREV = Magma_tally2Columnwise:         STOREV = Magma_tally2Rowwise:

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1  0  0 )
                     ( v1 v2 v3 )                     ( v2 v2 v2  1  0 )
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
                     (  0  1 v3 )
                     (  0  0  1 )

    @ingroup magma_tally2_caux3
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_clarfb_gpu_gemm(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_const_ptr dV,    magma_tally2_int_t lddv,
    magma_tally2FloatComplex_const_ptr dT,    magma_tally2_int_t lddt,
    magma_tally2FloatComplex_ptr dC,          magma_tally2_int_t lddc,
    magma_tally2FloatComplex_ptr dwork,       magma_tally2_int_t ldwork,
    magma_tally2FloatComplex_ptr dworkvt,     magma_tally2_int_t ldworkvt)
{
    magma_tally2FloatComplex c_zero    = MAGMA_tally2_C_ZERO;
    magma_tally2FloatComplex c_one     = MAGMA_tally2_C_ONE;
    magma_tally2FloatComplex c_neg_one = MAGMA_tally2_C_NEG_ONE;

    magma_tally2_int_t info = 0;
    
    /* Function Body */
    if (m <= 0 || n <= 0) {
        return info;
    }
    
    // internal variable
    magma_tally2_int_t ldwvt = (m > n ?  k : m);
    magma_tally2_int_t ldw;
    if ( side == Magma_tally2Left ) {
        ldw = k;
    } else {
        ldw = m;
    }
    
    // opposite of trans
    magma_tally2_trans_t transt;
    if (trans == Magma_tally2NoTrans)
        transt = Magma_tally2_ConjTrans;
    else
        transt = Magma_tally2NoTrans;
    
    MAGMA_tally2_UNUSED( transt );  // TODO: is this a bug that it isn't used?
    
    // whether V is stored transposed or not
    magma_tally2_trans_t notransV, transV;
    if (storev == Magma_tally2Columnwise) {
        notransV = Magma_tally2NoTrans;
        transV   = Magma_tally2_ConjTrans;
    }
    else {
        notransV = Magma_tally2_ConjTrans;
        transV   = Magma_tally2NoTrans;
    }

    if ( side == Magma_tally2Left ) {
        // Form H C or H^H C
        // Comments assume H C.
        // When forming H^H C, T gets transposed via transt for m >= n or by trans for m < n.
        
        // W = V^H C
        magma_tally2_cgemm( Magma_tally2_ConjTrans, notransV,
                     k, n, m,
                     c_one,  dV,    lddv,
                             dC,    lddc,
                     c_zero, dwork, ldw);

        if (m <= n) {
            // W2 = V T
            magma_tally2_cgemm( notransV, trans,
                         m, k, k,
                         c_one,  dV, lddv,
                                 dT, lddt,
                         c_zero, dworkvt, ldwvt);
            // C = C - W2 W = C - V T V^H C = (I - V T V^H) C = H C
            magma_tally2_cgemm( Magma_tally2NoTrans, Magma_tally2NoTrans,
                         m, n, k,
                         c_neg_one, dworkvt,  ldwvt,
                                    dwork,    ldw,
                         c_one,     dC,       lddc);
        } else {
            // W2 = T W  = T  V^H C
            magma_tally2_cgemm( trans, Magma_tally2NoTrans,
                         k, n, k,
                         c_one,  dT, lddt,
                                 dwork, ldw,
                         c_zero, dworkvt, ldwvt);
            // C = C - V W2 = C - V T V^H C = (I - V T V^H) C = H C
            magma_tally2_cgemm( notransV, Magma_tally2NoTrans,
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
        magma_tally2_cgemm( Magma_tally2NoTrans, notransV,
                     m, k, n,
                     c_one,  dC,    lddc,
                             dV,    lddv,
                     c_zero, dwork, ldw);
        if (m <= n) {
            // W2 = W T = C V T
            magma_tally2_cgemm( Magma_tally2NoTrans, trans,
                         m, k, k,
                         c_one,  dwork, ldw,
                                 dT, lddt,
                         c_zero, dworkvt, ldwvt);
            // C = C - W2 V^H = C - C V T V^H = C (I - V T V^H) = C H
            magma_tally2_cgemm( Magma_tally2NoTrans, transV,
                         m, n, k,
                         c_neg_one, dworkvt, ldwvt,
                                    dV,    lddv,
                         c_one,     dC,    lddc);
        } else {
            // W2 = T V^H
            magma_tally2_cgemm( trans, transV,
                         k, n, k,
                         c_one,  dT, lddt,
                                 dV, lddv,
                         c_zero, dworkvt, ldwvt);
            // C = C - W W2 = C - C V T V^H = C (I - V T V^H) = C H
            magma_tally2_cgemm( Magma_tally2NoTrans, Magma_tally2NoTrans,
                         m, n, k,
                         c_neg_one, dwork,   ldw,
                                    dworkvt, ldwvt,
                         c_one,     dC,      lddc);
        }
    }

    return info;
} /* magma_tally2_clarfb */
