/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from dgemm_tesla.cu normal d -> s, Fri Jan 30 19:00:10 2015
*/
#include "common_magma_minproduct.h"
#include "commonblas_s.h"

/**
    Purpose
    -------
    
    SGEMM performs one of the matrix-matrix operations
    
        C = alpha*op( A )*op( B ) + beta*C,
    
    where op( X ) is one of
    
        op( X ) = X   or   op( X ) = X**T,
    
    alpha and beta are scalars, and A, B and C are matrices, with op( A )
    an m by k matrix,  op( B ) a k by n matrix and C an m by n matrix.
    
    Parameters
    ----------
    
    @param[in]
    transA  magma_minproduct_trans_t.
            On entry, transA specifies the form of op( A ) to be used in
            the matrix multiplication as follows:
      -     = Magma_minproductNoTrans:    op( A ) = A.
      -     = Magma_minproductTrans:      op( A ) = A**T.
      -     = Magma_minproductConjTrans:  op( A ) = A**T.
    
    @param[in]
    transB  magma_minproduct_trans_t.
            On entry, transB specifies the form of op( B ) to be used in
            the matrix multiplication as follows:
      -     = Magma_minproductNoTrans:    op( B ) = B.
      -     = Magma_minproductTrans:      op( B ) = B**T.
      -     = Magma_minproductConjTrans:  op( B ) = B**T.
    
    @param[in]
    m       INTEGER.
            On entry,  M  specifies  the number  of rows  of the  matrix
            op( A )  and of the  matrix  C.  M  must  be at least  zero.
    
    @param[in]
    n       INTEGER.
            On entry,  N  specifies the number  of columns of the matrix
            op( B ) and the number of columns of the matrix C. N must be
            at least zero.
    
    @param[in]
    k       INTEGER.
            On entry,  K  specifies  the number of columns of the matrix
            op( A ) and the number of rows of the matrix op( B ). K must
            be at least  zero.
    
    @param[in]
    alpha   REAL.
            On entry, ALPHA specifies the scalar alpha.
    
    @param[in]
    A       REAL array of DIMENSION ( LDA, ka ), where ka is
            k  when  transA = Magma_minproductNoTrans,  and is  m  otherwise.
            Before entry with  transA = Magma_minproductNoTrans,  the leading  m by k
            part of the array  A  must contain the matrix  A,  otherwise
            the leading  k by m  part of the array  A  must contain  the
            matrix A.
    
    @param[in]
    lda     INTEGER.
            On entry, LDA specifies the first dimension of A as declared
            in the calling (sub) program. When  transA = Magma_minproductNoTrans then
            LDA must be at least  max( 1, m ), otherwise  LDA must be at
            least  max( 1, k ).
    
    @param[in]
    B       REAL array of DIMENSION ( LDB, kb ), where kb is
            n  when  transB = Magma_minproductNoTrans,  and is  k  otherwise.
            Before entry with  transB = Magma_minproductNoTrans,  the leading  k by n
            part of the array  B  must contain the matrix  B,  otherwise
            the leading  n by k  part of the array  B  must contain  the
            matrix B.
    
    @param[in]
    ldb     INTEGER.
            On entry, LDB specifies the first dimension of B as declared
            in the calling (sub) program. When  transB = Magma_minproductNoTrans then
            LDB must be at least  max( 1, k ), otherwise  LDB must be at
            least  max( 1, n ).
    
    @param[in]
    beta    REAL.
            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
            supplied as zero then C need not be set on input.
    
    @param[in,out]
    C       REAL array of DIMENSION ( LDC, n ).
            Before entry, the leading  m by n  part of the array  C must
            contain the matrix  C,  except when  beta  is zero, in which
            case C need not be set on entry.
            On exit, the array  C  is overwritten by the  m by n  matrix
            ( alpha*op( A )*op( B ) + beta*C ).
    
    @param[in]
    ldc     INTEGER.
            On entry, LDC specifies the first dimension of C as declared
            in  the  calling  (sub)  program.   LDC  must  be  at  least
            max( 1, m ).

    @ingroup magma_minproduct_sblas3
    ********************************************************************/
extern "C" void
magma_minproductblas_sgemm_tesla(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    const float *A, magma_minproduct_int_t lda,
    const float *B, magma_minproduct_int_t ldb,
    float beta,
    float *C, magma_minproduct_int_t ldc )
{
    magma_minproduct_int_t info = 0;
    if      ( transA != Magma_minproductNoTrans && transA != Magma_minproductTrans && transA != Magma_minproductConjTrans )
        info = -1;
    else if ( transB != Magma_minproductNoTrans && transB != Magma_minproductTrans && transB != Magma_minproductConjTrans )
        info = -2;
    else if ( m < 0 )
        info = -3;
    else if ( n < 0 )
        info = -4;
    else if ( k < 0 )
        info = -5;
    else if ( transA == Magma_minproductNoTrans ? lda < m : lda < k )
        info = -8;
    else if ( transB == Magma_minproductNoTrans ? ldb < k : lda < n )
        info = -10;
    else if ( ldc < m )
        info = -13;
    
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    if ( m == 0 || n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0) ) {
        return;
    }
    if ( alpha == 0.0 ) {
        if ( beta == 0.0 ) {
            magma_minproductblas_sgemm_ab_0(
                C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
            return;
        }
        else {
            magma_minproductblas_sgemm_a_0(
                C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
            return;
        }
    }
    
    if ( transA == Magma_minproductNoTrans ) {
        if ( transB == Magma_minproductNoTrans ) {
            /*=======================================================================
              ===================C = alpha * A * B + beta * C =======================
              =======================================================================*/
            if ( m > 512 && n > 512 ) {
                if ( m % 64 == 0 && n % 16 == 0 && k % 16 == 0 )
                    magma_minproductblas_sgemm_N_N_64_16_16_16_4_special(
                        C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
                else
                    magma_minproductblas_sgemm_N_N_64_16_16_16_4(
                        C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
            }
            else {
                magma_minproduct_sgemm(transA, transB,
                    m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
            }
        }
        else {
            /*=======================================================================
              ===================C = alpha * A * B^T + beta * C======================
              =======================================================================*/
            if ( m > 512 && n > 512 ) {
                //if ( m % 64 == 0 && n % 16 == 0 && k % 4 == 0 )
                //    magma_minproductblas_sgemm_N_T_64_16_4_16_4(
                //        C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
                //else
                    magma_minproductblas_sgemm_N_T_64_16_4_16_4(
                        C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
            }
            else {
                magma_minproduct_sgemm(transA, transB,
                    m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
            }
        }
    }
    else {
        if ( transB == Magma_minproductNoTrans ) {
            /*=======================================================================
              ===================C = alpha * A^T * B + beta * C======================
              =======================================================================*/
            if ( m > 512 && n > 512 ) {
                //if ( m % 32 == 0 && n % 32 == 0 && k % 8 == 0 )
                //    magma_minproductblas_sgemm_T_N_32_32_8_8_8(
                //        C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
                //else
                    magma_minproductblas_sgemm_T_N_32_32_8_8_8(
                        C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
            }
            else {
                magma_minproduct_sgemm(transA, transB,
                    m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
            }
        }
        else {
            /*=======================================================================
              ===================C = alpha * A^T * B^T + beta * C====================
              =======================================================================*/
            if ( m > 512 && n > 512 ) {
                if ( m % 64 == 0 && n % 16 == 0 && k % 16 == 0 )
                    magma_minproductblas_sgemm_T_T_64_16_16_16_4_special(
                        C, B, A, n, m, k, ldb, lda, ldc, alpha, beta );
                else
                    magma_minproductblas_sgemm_T_T_64_16_16_16_4(
                        C, B, A, n, m, k, ldb, lda, ldc, alpha, beta );
            }
            else {
                magma_minproduct_sgemm(transA, transB,
                    m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
            }
        }
    }
}
