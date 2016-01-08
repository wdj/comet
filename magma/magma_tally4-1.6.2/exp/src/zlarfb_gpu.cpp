/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

// === Define what BLAS to use ============================================
#define PRECISION_z
#if (defined(PRECISION_s) || defined(PRECISION_d))
//  #define cublasZgemm magma_tally4blas_zgemm
#endif
// === End defining what BLAS to use =======================================

extern "C" magma_tally4_int_t
magma_tally4_zlarfb_gpu( char side, char trans, char direct, char storev,
          magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
          cuDoubleComplex *dV,    magma_tally4_int_t ldv, 
          cuDoubleComplex *dT,    magma_tally4_int_t ldt,
          cuDoubleComplex *dC,    magma_tally4_int_t ldc, 
          cuDoubleComplex *dwork, magma_tally4_int_t ldwork)
{
/*  -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Univ. of California Berkeley
       @date January 2015

    Purpose
    =======
    ZLARFB applies a complex block reflector H or its transpose H' to a
    COMPLEX_16 m by n matrix C, from the left.

    Arguments
    =========
    SIDE    (input) CHARACTER
            = 'L': apply H or H' from the Left
            = 'R': apply H or H' from the Right (Not implemented)

    TRANS   (input) CHARACTER
            = 'N': apply H  (No transpose)      (Not implemented)
            = 'C': apply H' (Conjugate transpose)

    DIRECT  (input) CHARACTER
            Indicates how H is formed from a product of elementary
            reflectors
            = 'F': H = H(1) H(2) . . . H(k) (Forward)
            = 'B': H = H(k) . . . H(2) H(1) (Backward)

    STOREV  (input) CHARACTER
            Indicates how the vectors which define the elementary
            reflectors are stored:
            = 'C': Columnwise
            = 'R': Rowwise

    M       (input) INTEGER
            The number of rows of the matrix C.

    N       (input) INTEGER
            The number of columns of the matrix C.

    K       (input) INTEGER
            The order of the matrix T (= the number of elementary
            reflectors whose product defines the block reflector).

    DV      (input) COMPLEX_16 array, dimension (LDV,K)
            The matrix V. See further details.

    LDV     (input) INTEGER
            The leading dimension of the array V. LDV >= max(1,M);

    DT      (input) COMPLEX_16 array, dimension (LDT,K)
            The triangular k by k matrix T in the representation of the
            block reflector.

    LDT     (input) INTEGER
            The leading dimension of the array T. LDT >= K.

    DC      (input/output) COMPLEX_16 array, dimension (LDC,N)
            On entry, the m by n matrix C.
            On exit, C is overwritten by H*C.

    LDC     (input) INTEGER
            The leading dimension of the array C. LDA >= max(1,M).

    WORK    (workspace) COMPLEX_16 array, dimension (LDWORK,K)

    LDWORK  (input) INTEGER
            The leading dimension of the array WORK. LDWORK >= max(1,N);
    ===================================================================      */

    cuDoubleComplex c_zero    = MAGMA_tally4_Z_ZERO;
    cuDoubleComplex c_one     = MAGMA_tally4_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;

    /* Function Body */
    if (m <= 0 || n <= 0) {
        return MAGMA_tally4_SUCCESS;
    }

    char transt;
    if (trans == 'N' || trans == 'n')
      transt = Magma_tally4ConjTrans;
    else
      transt = Magma_tally4NoTrans;

    if ( ( side  == 'r' || side  == 'R') ) {
        fprintf(stderr, "The case (side == right) is not implemented\n");
        magma_tally4_xerbla( __func__, 1 );
        return MAGMA_tally4_ERR_ILLEGAL_VALUE;
    }

    if ( storev == 'c' || storev == 'C') {
        /*
          if (n==1 && m%32==0){
          // This is used when we have to apply H on only one vector
          magma_tally4blas_zgemvt(m, k, 1., dv_ref(0,0), ldv, dc_ref(0, 0), dwork);
          printf("m= %d, n = %d, ldwork = %d\n", m, k, ldwork);
          }
          else
        */
        cublasZgemm( Magma_tally4ConjTrans, Magma_tally4NoTrans,
                     n, k, m,
                     c_one,  dC,    ldc,
                             dV,    ldv,
                     c_zero, dwork, ldwork);

        if (direct == 'F' || direct =='f')
            cublasZtrmm( Magma_tally4Right, Magma_tally4Upper, transt, Magma_tally4NonUnit,
                         n, k, 
                         c_one, dT,    ldt, 
                                dwork, ldwork);
        else
            cublasZtrmm( Magma_tally4Right, Magma_tally4Lower, transt, Magma_tally4NonUnit,
                         n, k, 
                         c_one, dT,    ldt, 
                                dwork, ldwork);

        cublasZgemm( Magma_tally4NoTrans, Magma_tally4ConjTrans, 
                     m, n, k, 
                     c_neg_one, dV,    ldv,
                                dwork, ldwork, 
                     c_one,     dC,    ldc);
    }
    else {
        cublasZgemm( Magma_tally4NoTrans, Magma_tally4ConjTrans, 
                     m, k, n, 
                     c_one,  dC,    ldc,
                             dV,    ldv, 
                     c_zero, dwork, ldwork);

        cublasZtrmm( Magma_tally4Right, Magma_tally4Upper, transt, Magma_tally4NonUnit,
                     m, k, 
                     c_one, dT,    ldt, 
                            dwork, ldwork);

        cublasZgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, 
                     m, n, k, 
                     c_neg_one, dwork, ldwork,
                                dV,    ldv,
                     c_one,     dC,    ldc);
    }
    return MAGMA_tally4_SUCCESS;
} /* magma_tally4_zlarfb */
