/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @author Azzam Haidar
       @precisions normal z -> s d c
*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    ZLARFB applies a complex block reflector H or its transpose H^H to a
    COMPLEX_16 m by n matrix C, from the left.
    
    __Note that this function assumes__ that the upper part of dV is 0
    because it is referenced. Same for upper/lower part of dT.

    Arguments
    ---------
    @param[in]
    side    magma_tally3_side_t
      -     = Magma_tally3Left:      apply H or H^H from the Left
      -     = Magma_tally3Right:     apply H or H^H from the Right

    @param[in]
    trans   magma_tally3_trans_t
      -     = Magma_tally3NoTrans:    apply H   (No transpose)
      -     = Magma_tally3_ConjTrans: apply H^H (Conjugate transpose)

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
    dV      COMPLEX_16 array on the GPU, dimension
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
    dT      COMPLEX_16 array on the GPU, dimension (LDDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    @param[in]
    lddt    INTEGER
            The leading dimension of the array T. LDDT >= K.

    @param[in,out]
    dC      COMPLEX_16 array on the GPU, dimension (LDDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C, or H^H*C, or C*H, or C*H^H.

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDA >= max(1,M).

    @param
    dwork   (workspace) COMPLEX_16 array, dimension (LDWORK,K)

    @param[in]
    ldwork  INTEGER
            The leading dimension of the array WORK.
            If SIDE = Magma_tally3Left,  LDWORK >= max(1,N);
            if SIDE = Magma_tally3Right, LDWORK >= max(1,M);

    @param
    dworkvt (workspace) COMPLEX_16 array, dimension (LDWORKT,K)

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

    @ingroup magma_tally3_zaux3
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zlarfb_gpu_gemm(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_const_ptr dV,    magma_tally3_int_t lddv,
    magma_tally3DoubleComplex_const_ptr dT,    magma_tally3_int_t lddt,
    magma_tally3DoubleComplex_ptr dC,          magma_tally3_int_t lddc,
    magma_tally3DoubleComplex_ptr dwork,       magma_tally3_int_t ldwork,
    magma_tally3DoubleComplex_ptr dworkvt,     magma_tally3_int_t ldworkvt)
{
    magma_tally3DoubleComplex c_zero    = MAGMA_tally3_Z_ZERO;
    magma_tally3DoubleComplex c_one     = MAGMA_tally3_Z_ONE;
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;

    magma_tally3_int_t info = 0;
    
    /* Function Body */
    if (m <= 0 || n <= 0) {
        return info;
    }
    
    // internal variable
    magma_tally3_int_t ldwvt = (m > n ?  k : m);
    magma_tally3_int_t ldw;
    if ( side == Magma_tally3Left ) {
        ldw = k;
    } else {
        ldw = m;
    }
    
    // opposite of trans
    magma_tally3_trans_t transt;
    if (trans == Magma_tally3NoTrans)
        transt = Magma_tally3_ConjTrans;
    else
        transt = Magma_tally3NoTrans;
    
    MAGMA_tally3_UNUSED( transt );  // TODO: is this a bug that it isn't used?
    
    // whether V is stored transposed or not
    magma_tally3_trans_t notransV, transV;
    if (storev == Magma_tally3Columnwise) {
        notransV = Magma_tally3NoTrans;
        transV   = Magma_tally3_ConjTrans;
    }
    else {
        notransV = Magma_tally3_ConjTrans;
        transV   = Magma_tally3NoTrans;
    }

    if ( side == Magma_tally3Left ) {
        // Form H C or H^H C
        // Comments assume H C.
        // When forming H^H C, T gets transposed via transt for m >= n or by trans for m < n.
        
        // W = V^H C
        magma_tally3_zgemm( Magma_tally3_ConjTrans, notransV,
                     k, n, m,
                     c_one,  dV,    lddv,
                             dC,    lddc,
                     c_zero, dwork, ldw);

        if (m <= n) {
            // W2 = V T
            magma_tally3_zgemm( notransV, trans,
                         m, k, k,
                         c_one,  dV, lddv,
                                 dT, lddt,
                         c_zero, dworkvt, ldwvt);
            // C = C - W2 W = C - V T V^H C = (I - V T V^H) C = H C
            magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3NoTrans,
                         m, n, k,
                         c_neg_one, dworkvt,  ldwvt,
                                    dwork,    ldw,
                         c_one,     dC,       lddc);
        } else {
            // W2 = T W  = T  V^H C
            magma_tally3_zgemm( trans, Magma_tally3NoTrans,
                         k, n, k,
                         c_one,  dT, lddt,
                                 dwork, ldw,
                         c_zero, dworkvt, ldwvt);
            // C = C - V W2 = C - V T V^H C = (I - V T V^H) C = H C
            magma_tally3_zgemm( notransV, Magma_tally3NoTrans,
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
        magma_tally3_zgemm( Magma_tally3NoTrans, notransV,
                     m, k, n,
                     c_one,  dC,    lddc,
                             dV,    lddv,
                     c_zero, dwork, ldw);
        if (m <= n) {
            // W2 = W T = C V T
            magma_tally3_zgemm( Magma_tally3NoTrans, trans,
                         m, k, k,
                         c_one,  dwork, ldw,
                                 dT, lddt,
                         c_zero, dworkvt, ldwvt);
            // C = C - W2 V^H = C - C V T V^H = C (I - V T V^H) = C H
            magma_tally3_zgemm( Magma_tally3NoTrans, transV,
                         m, n, k,
                         c_neg_one, dworkvt, ldwvt,
                                    dV,    lddv,
                         c_one,     dC,    lddc);
        } else {
            // W2 = T V^H
            magma_tally3_zgemm( trans, transV,
                         k, n, k,
                         c_one,  dT, lddt,
                                 dV, lddv,
                         c_zero, dworkvt, ldwvt);
            // C = C - W W2 = C - C V T V^H = C (I - V T V^H) = C H
            magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3NoTrans,
                         m, n, k,
                         c_neg_one, dwork,   ldw,
                                    dworkvt, ldwvt,
                         c_one,     dC,      lddc);
        }
    }

    return info;
} /* magma_tally3_zlarfb */
