/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @author Azzam Haidar
       @author Tingxing Dong
       @generated from zlarfb_gemm_batched.cpp normal z -> d, Fri Jan 30 19:00:19 2015
*/


#include "../testing/testings.h"  // muse bed included in order to use cublas_trans_const_tally4
#include "common_magma_tally4.h"

//#define USE_CUBLAS

/**
    Purpose
    -------
    DLARFB applies a real block reflector H or its transpose H^H to a
    DOUBLE_PRECISION m by n matrix C, from the left.
    
    __Note that this function assumes__ that the upper part of dV_array is 0
    because it is referenced. Same for upper/lower part of dT_array.

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
    dV_array      DOUBLE_PRECISION array on the GPU, dimension
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
    dT_array      DOUBLE_PRECISION array on the GPU, dimension (LDDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    @param[in]
    lddt    INTEGER
            The leading dimension of the array T. LDDT >= K.

    @param[in,out]
    dC_array      DOUBLE_PRECISION array on the GPU, dimension (LDDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C, or H^H*C, or C*H, or C*H^H.

    @param[in]
    lddc    INTEGER
            The leading dimension of the array C. LDA >= max(1,M).

    @param
    dwork_array   (workspace) DOUBLE_PRECISION array, dimension (LDWORK,K)

    @param[in]
    ldwork  INTEGER
            The leading dimension of the array WORK.
            If SIDE = Magma_tally4Left,  LDWORK >= max(1,N);
            if SIDE = Magma_tally4Right, LDWORK >= max(1,M);

    @param
    dworkvt_array (workspace) DOUBLE_PRECISION array, dimension (LDWORKT,K)

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

    @ingroup magma_tally4_daux3
    ********************************************************************/
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally4_int_t
magma_tally4_dlarfb_gemm_batched_magem(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_const_ptr dV_array[],    magma_tally4_int_t lddv,
    magma_tally4Double_const_ptr dT_array[],    magma_tally4_int_t lddt,
    magma_tally4Double_ptr dC_array[],          magma_tally4_int_t lddc,
    magma_tally4Double_ptr dwork_array[],       magma_tally4_int_t ldwork,
    magma_tally4Double_ptr dworkvt_array[],     magma_tally4_int_t ldworkvt,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue);
extern "C" magma_tally4_int_t
magma_tally4_dlarfb_gemm_batched_cugem(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_const_ptr dV_array[],    magma_tally4_int_t lddv,
    magma_tally4Double_const_ptr dT_array[],    magma_tally4_int_t lddt,
    magma_tally4Double_ptr dC_array[],          magma_tally4_int_t lddc,
    magma_tally4Double_ptr dwork_array[],       magma_tally4_int_t ldwork,
    magma_tally4Double_ptr dworkvt_array[],     magma_tally4_int_t ldworkvt,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally4_int_t
magma_tally4_dlarfb_gemm_batched(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_const_ptr dV_array[],    magma_tally4_int_t lddv,
    magma_tally4Double_const_ptr dT_array[],    magma_tally4_int_t lddt,
    magma_tally4Double_ptr dC_array[],          magma_tally4_int_t lddc,
    magma_tally4Double_ptr dwork_array[],       magma_tally4_int_t ldwork,
    magma_tally4Double_ptr dworkvt_array[],     magma_tally4_int_t ldworkvt,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue)
{

    if(m >= 32 && n >= 32 && k >= 32)
    {
        magma_tally4_dlarfb_gemm_batched_magem(
            side, trans, direct, storev,
            m, n,  k,
            (magma_tally4Double_const_ptr *) dV_array, lddv,
            (magma_tally4Double_const_ptr *) dT_array, lddt,
            dC_array, lddc,
            dwork_array, ldwork,
            dworkvt_array, ldworkvt,
            batchCount, queue);
    }
    else{

        magma_tally4_dlarfb_gemm_batched_cugem(
            side, trans, direct, storev,
            m, n,  k,
            (magma_tally4Double_const_ptr *) dV_array, lddv,
            (magma_tally4Double_const_ptr *) dT_array, lddt,
            dC_array, lddc,
            dwork_array, ldwork,
            dworkvt_array, ldworkvt,
            batchCount, myhandle, queue);
   }
    return MAGMA_tally4_SUCCESS;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally4_int_t
magma_tally4_dlarfb_gemm_batched_magem(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_const_ptr dV_array[],    magma_tally4_int_t lddv,
    magma_tally4Double_const_ptr dT_array[],    magma_tally4_int_t lddt,
    magma_tally4Double_ptr dC_array[],          magma_tally4_int_t lddc,
    magma_tally4Double_ptr dwork_array[],       magma_tally4_int_t ldwork,
    magma_tally4Double_ptr dworkvt_array[],     magma_tally4_int_t ldworkvt,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue)
{
    double c_zero    = MAGMA_tally4_D_ZERO;
    double c_one     = MAGMA_tally4_D_ONE;
    double c_neg_one = MAGMA_tally4_D_NEG_ONE;

    /* Function Body */
    magma_tally4_int_t info = 0;
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
    
    // whether T is upper or lower triangular
    magma_tally4_uplo_t uplo;
    if (direct == Magma_tally4Forward)
        uplo = Magma_tally4Upper;
    else
        uplo = Magma_tally4Lower;
    
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
        
        // W = V' C                              
        magma_tally4blas_dgemm_batched( Magma_tally4Trans,notransV, /*NontransLeft*/
                     k, n, m,
                     c_one,  dV_array,    lddv,
                             dC_array,    lddc,
                     c_zero, dwork_array, ldw, batchCount, queue);

        if (m <= n) {
            // W2 = V T
            magma_tally4blas_dgemm_batched( notransV, trans, /* (NoTrans), trans(ConjTrans),*/
                         m, k, k,
                         c_one,  dV_array, lddv,
                                 dT_array, lddt,
                         c_zero, dworkvt_array, ldwvt, batchCount, queue);


            // C = C - W2 W = C - V T V' C = (I - V T V') C = H C
            magma_tally4blas_dgemm_batched( Magma_tally4NoTrans, Magma_tally4NoTrans,
                         m, n, k,
                         c_neg_one, dworkvt_array,  ldwvt,
                                    dwork_array,    ldw,
                         c_one,     dC_array,       lddc, batchCount, queue);
        } else 
        {
            // W2 = T W  = T  V' C
            magma_tally4blas_dgemm_batched( trans, Magma_tally4NoTrans,
                         k, n, k,
                         c_one,  dT_array, lddt,
                                 dwork_array, ldw,
                         c_zero, dworkvt_array, ldwvt, batchCount, queue);

            // C = C - V W2 = C - V T V' C = (I - V T V') C = H C
            magma_tally4blas_dgemm_batched( notransV, Magma_tally4NoTrans,
                         m, n, k,
                         c_neg_one, dV_array,  lddv,
                                    dworkvt_array,  ldwvt,
                         c_one,     dC_array,       lddc, batchCount, queue);
        }
    }
    else {
        // Form C H or C H^H
        // Comments assume C H.
        // When forming C H^H, T gets transposed via trans.
        
        // W = C V
        magma_tally4blas_dgemm_batched( Magma_tally4NoTrans, notransV,
                     m, k, n,
                     c_one,  dC_array,    lddc,
                             dV_array,    lddv,
                     c_zero, dwork_array, ldw, batchCount, queue);
        if (m <= n) {
            // W2 = W T = C V T
            magma_tally4blas_dgemm_batched( Magma_tally4NoTrans, trans,
                         m, k, k,
                         c_one,  dwork_array, ldw,
                                 dT_array, lddt,
                         c_zero, dworkvt_array, ldwvt, batchCount, queue);

            // C = C - W2 V' = C - C V T V' = C (I - V T V') = C H
            magma_tally4blas_dgemm_batched( Magma_tally4NoTrans, transV,
                         m, n, k,
                         c_neg_one, dworkvt_array, ldwvt,
                                    dV_array,    lddv,
                         c_one,     dC_array,    lddc, batchCount, queue);
        } else {
            // W2 = T V'
            magma_tally4blas_dgemm_batched( trans, transV,
                         k, n, k,
                         c_one,  dT_array, lddt,
                                 dV_array, lddv,
                         c_zero, dworkvt_array, ldwvt, batchCount, queue);
            // C = C - W W2 = C - C V T V' = C (I - V T V') = C H
            magma_tally4blas_dgemm_batched( Magma_tally4NoTrans, Magma_tally4NoTrans,
                         m, n, k,
                         c_neg_one, dwork_array,   ldw,
                                    dworkvt_array, ldwvt,
                         c_one,     dC_array,      lddc, batchCount, queue);
        }
    }
    return MAGMA_tally4_SUCCESS;
} /* magma_tally4_dlarfb */
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally4_int_t
magma_tally4_dlarfb_gemm_batched_cugem(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_const_ptr dV_array[],    magma_tally4_int_t lddv,
    magma_tally4Double_const_ptr dT_array[],    magma_tally4_int_t lddt,
    magma_tally4Double_ptr dC_array[],          magma_tally4_int_t lddc,
    magma_tally4Double_ptr dwork_array[],       magma_tally4_int_t ldwork,
    magma_tally4Double_ptr dworkvt_array[],     magma_tally4_int_t ldworkvt,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue)
{
    double c_zero    = MAGMA_tally4_D_ZERO;
    double c_one     = MAGMA_tally4_D_ONE;
    double c_neg_one = MAGMA_tally4_D_NEG_ONE;

    /* Function Body */
    magma_tally4_int_t info = 0;
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
    
    // whether T is upper or lower triangular
    magma_tally4_uplo_t uplo;
    if (direct == Magma_tally4Forward)
        uplo = Magma_tally4Upper;
    else
        uplo = Magma_tally4Lower;
    
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
        
        // W = V' C
        cublasDgemmBatched(myhandle, cublas_trans_const_tally4(Magma_tally4Trans), cublas_trans_const_tally4(notransV),
                     k, n, m,
                     &c_one,  (const double**)dV_array,    lddv,
                              (const double**)dC_array,    lddc,
                     &c_zero, dwork_array, ldw, batchCount);
        if (m <= n) {
            // W2 = V T
            cublasDgemmBatched(myhandle, cublas_trans_const_tally4(notransV), cublas_trans_const_tally4(trans),
                         m, k, k,
                         &c_one,  (const double**)dV_array, lddv,
                                  (const double**)dT_array, lddt,
                         &c_zero, dworkvt_array, ldwvt, batchCount);
            // C = C - W2 W = C - V T V' C = (I - V T V') C = H C
            cublasDgemmBatched(myhandle, cublas_trans_const_tally4(Magma_tally4NoTrans), cublas_trans_const_tally4(Magma_tally4NoTrans),
                         m, n, k,
                         &c_neg_one, (const double**)dworkvt_array,  ldwvt,
                                     (const double**)dwork_array,    ldw,
                         &c_one,     dC_array,       lddc, batchCount);
        } else {
            // W2 = T W  = T  V' C
            cublasDgemmBatched(myhandle, cublas_trans_const_tally4(trans), cublas_trans_const_tally4(Magma_tally4NoTrans),
                         k, n, k,
                         &c_one,  (const double**)dT_array, lddt,
                                  (const double**)dwork_array, ldw,
                         &c_zero, dworkvt_array, ldwvt, batchCount);
            // C = C - V W2 = C - V T V' C = (I - V T V') C = H C
            cublasDgemmBatched(myhandle, cublas_trans_const_tally4(notransV), cublas_trans_const_tally4(Magma_tally4NoTrans),
                         m, n, k,
                         &c_neg_one, (const double**)dV_array,  lddv,
                                     (const double**)dworkvt_array,  ldwvt,
                         &c_one,     dC_array,       lddc, batchCount);
        }
    }
    else {
        // Form C H or C H^H
        // Comments assume C H.
        // When forming C H^H, T gets transposed via trans.
        
        // W = C V
        cublasDgemmBatched(myhandle, cublas_trans_const_tally4(Magma_tally4NoTrans), cublas_trans_const_tally4(notransV),
                     m, k, n,
                     &c_one,  (const double**)dC_array,    lddc,
                              (const double**)dV_array,    lddv,
                     &c_zero, dwork_array, ldw, batchCount);
        if (m <= n) {
            // W2 = W T = C V T
           cublasDgemmBatched(myhandle, cublas_trans_const_tally4(Magma_tally4NoTrans), cublas_trans_const_tally4(trans),
                         m, k, k,
                         &c_one,  (const double**)dwork_array, ldw,
                                  (const double**)dT_array, lddt,
                         &c_zero, dworkvt_array, ldwvt, batchCount);
            // C = C - W2 V' = C - C V T V' = C (I - V T V') = C H
            cublasDgemmBatched(myhandle, cublas_trans_const_tally4(Magma_tally4NoTrans), cublas_trans_const_tally4(transV),
                         m, n, k,
                         &c_neg_one, (const double**)dworkvt_array, ldwvt,
                                     (const double**)dV_array,    lddv,
                         &c_one,     dC_array,    lddc, batchCount);
        } else {
            // W2 = T V'
            cublasDgemmBatched(myhandle, cublas_trans_const_tally4(trans), cublas_trans_const_tally4(transV),
                         k, n, k,
                         &c_one,  (const double**)dT_array, lddt,
                                  (const double**)dV_array, lddv,
                         &c_zero, dworkvt_array, ldwvt, batchCount);
            // C = C - W W2 = C - C V T V' = C (I - V T V') = C H
           cublasDgemmBatched(myhandle, cublas_trans_const_tally4(Magma_tally4NoTrans), cublas_trans_const_tally4(Magma_tally4NoTrans),
                         m, n, k,
                         &c_neg_one, (const double**)dwork_array,   ldw,
                                     (const double**)dworkvt_array, ldwvt,
                         &c_one,     dC_array,      lddc, batchCount);
        }
    }
    return MAGMA_tally4_SUCCESS;
} /* magma_tally4_dlarfb */

