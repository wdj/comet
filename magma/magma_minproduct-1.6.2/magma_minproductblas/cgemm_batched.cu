/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgemm_batched.cu normal z -> c, Fri Jan 30 19:00:10 2015

       @author Jakub Kurzak
       @author Stan Tomov
       @author Mark Gates
       @author Azzam Haidar

       [zcds]gemm_fermi.cu          defines the CPU driver.
       [zcds]gemm_fermi_kernels.h   defines the block sizes for each precision.
       gemm_stencil_defs.h          defines types and functions for precision-independent code.
       
       These files are included multiple times, once for each transpose version.
       gemm_stencil.cuh             defines the GPU kernel (device function).
       gemm_kernel_batched.cuh              defines the GPU kernel (global function).
       
       The batched version uses gemm_kernel_batched.cuh instead of gemm_kernel.cuh.
*/
#include "common_magma_minproduct.h"
#include "commonblas_c.h"

#define PRECISION_c

///////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    CGEMM performs one of the matrix-matrix operations
    
        C = alpha*op( A )*op( B ) + beta*C,
    
    where op( X ) is one of
    
        op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
    
    alpha and beta are scalars, and A, B and C are matrices, with
    op( A ) an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
    
    Parameters
    ----------
    @param[in]
    transA  CHARACTER*1.
            On entry, transA specifies the form of op( A ) to be used in
            the matrix multiplication as follows:
      -     = 'N':  op( A ) = A.
      -     = 'T':  op( A ) = A**T.
      -     = 'C':  op( A ) = A**H.
    
    @param[in]
    transB  CHARACTER*1.
            On entry, transB specifies the form of op( B ) to be used in
            the matrix multiplication as follows:
      -     = 'N':  op( B ) = B.
      -     = 'T':  op( B ) = B**T.
      -     = 'C':  op( B ) = B**H.
    
    @param[in]
    m       INTEGER.
            On entry,  M  specifies  the number  of rows  of the  matrix
            op( dA )  and of the  matrix dC.  M  must  be at least  zero.
    
    @param[in]
    n       INTEGER.
            On entry,  N  specifies the number  of columns of the matrix
            op( dB ) and the number of columns of the matrix dC. N must be
            at least zero.
    
    @param[in]
    k       INTEGER.
            On entry,  K  specifies  the number of columns of the matrix
            op( dA ) and the number of rows of the matrix op( dB ). K must
            be at least  zero.
    
    @param[in]
    alpha   COMPLEX
            On entry, ALPHA specifies the scalar alpha.
    
    @param[in]
    dA      COMPLEX array of DIMENSION ( LDA, ka ), where ka is
            k  when  transA = Magma_minproductNoTrans,  and is  m  otherwise.
            Before entry with  transA = Magma_minproductNoTrans,  the leading  m by k
            part of the array dA must contain the matrix dA, otherwise
            the leading  k by m  part of the array dA must contain  the
            matrix dA.
    
    @param[in]
    ldda    INTEGER.
            On entry, LDA specifies the first dimension of A as declared
            in the calling (sub) program. When  transA = Magma_minproductNoTrans then
            LDA must be at least  max( 1, m ), otherwise  LDA must be at
            least  max( 1, k ).
    
    @param[in]
    dB      COMPLEX array of DIMENSION ( LDB, kb ), where kb is
            n  when  transB = Magma_minproductNoTrans,  and is  k  otherwise.
            Before entry with  transB = Magma_minproductNoTrans,  the leading  k by n
            part of the array dB must contain the matrix dB, otherwise
            the leading  n by k  part of the array dB must contain  the
            matrix dB.
    
    @param[in]
    lddb    INTEGER.
            On entry, LDB specifies the first dimension of dB as declared
            in the calling (sub) program. When  transB = Magma_minproductNoTrans then
            LDB must be at least  max( 1, k ), otherwise  LDB must be at
            least  max( 1, n ).
    
    @param[in]
    beta    COMPLEX.
            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
            supplied as zero then dC need not be set on input.
    
    @param[in,out]
    dC      COMPLEX array of DIMENSION ( LDC, n ).
            Before entry, the leading  m by n  part of the array  dC must
            contain the matrix  dC,  except when  beta  is zero, in which
            case dC need not be set on entry.
            On exit, the array  dC  is overwritten by the  m by n  matrix
            ( alpha*op( dA )*op( dB ) + beta*dC ).
    
    @param[in]
    lddc    INTEGER.
            On entry, LDC specifies the first dimension of dC as declared
            in  the  calling  (sub)  program.   LDC  must  be  at  least
            max( 1, m ).

    @ingroup magma_minproduct_cblas3
    ********************************************************************/
extern "C" void
magma_minproductblas_cgemm_batched(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex const * const * dA_array, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex const * const * dB_array, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue )
{

    //if( (k < 32) || (k == 32 && m < 256) || (k == 32 && n < 256) || (k == 32 && transA==Magma_minproductTrans) || (k == 32 && transA==Magma_minproductConjTrans) || (m <= 32 && n <= 32) ) {
    if( (k <= 32) || (m <= 32 && n <= 32) ) {
        magma_minproductblas_cgemm_batched_k32(
                  transA, transB, m, n, k,
                  alpha, dA_array, ldda,
                         dB_array, lddb,
                  beta,  dC_array, lddc,
                  batchCount, queue );
    }
    else{
        magma_minproductblas_cgemm_batched_lg(
                  transA, transB, m, n, k,
                  alpha, dA_array, ldda,
                         dB_array, lddb,
                  beta,  dC_array, lddc,
                  batchCount, queue );
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
