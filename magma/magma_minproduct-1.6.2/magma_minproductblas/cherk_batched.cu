/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zherk_batched.cu normal z -> c, Fri Jan 30 19:00:10 2015

       @author Jakub Kurzak
       @author Stan Tomov
       @author Mark Gates
       @author Azzam Haidar

       [zcds]gemm_fermi.cu          defines the CPU driver.
       [zcds]gemm_fermi_kernels.h   defines the block sizes for each precision.
       gemm_stencil_defs.h          defines types and functions for precision-independent code.
       
       These files are included multiple times, once for each transpose version.
       herk_stencil.cuh             defines the GPU kernel (device function).
       herk_kernel_batched.cuh              defines the GPU kernel (global function).
       
       The batched version uses herk_kernel_batched.cuh instead of herk_kernel.cuh.
*/
#include "common_magma_minproduct.h"
#include "commonblas_c.h"

#define PRECISION_c

///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    CHERK performs one of the hermitian rank k operations

    C := alpha*A*A**H + beta*C,

    or

    C := alpha*A**H*A + beta*C,

    where alpha and beta are real scalars, C is an n by n hermitian
    matrix and A is an n by k matrix in the first case and a k by n
    matrix in the second case.
    
    Parameters
    ----------

    @param[in]
    uplo    CHARACTER*1.
           On entry, uplo specifies whether the upper or lower
           triangular part of the array C is to be referenced as
           follows:

           uplo = 'U' or 'u' Only the upper triangular part of C
           is to be referenced.

           uplo = 'L' or 'l' Only the lower triangular part of C
           is to be referenced.
    
    @param[in]
    trans   CHARACTER*1.
            On entry, trans specifies the operation to be performed as
            follows:

            trans = 'N' or 'n' C := alpha*A*A**H + beta*C.

            trans = 'C' or 'c' C := alpha*A**H*A + beta*C.

    @param[in]
    n       INTEGER.
            On entry,  specifies the order of the matrix C. N must be
            at least zero.
    
    @param[in]
    k       INTEGER.
            On entry with trans = 'N' or 'n', k specifies the number
            of columns of the matrix A, and on entry with
            trans = 'C' or 'c', k specifies the number of rows of the
            matrix A. K must be at least zero.

    @param[in]
    alpha   REAL
            On entry, ALPHA specifies the scalar alpha.
    
    @param[in]
    dA      COMPLEX array of DIMENSION ( ldda, ka ), where ka is
            k  when  trans = Magma_minproductNoTrans,  and is  n  otherwise.
            Before entry with  trans = Magma_minproductNoTrans,  the leading  m by k
            part of the array dA must contain the matrix dA, otherwise
            the leading  k by m  part of the array dA must contain  the
            matrix dA.
    
    @param[in]
    ldda    INTEGER.
            On entry, ldda specifies the first dimension of A as declared
            in the calling (sub) program. When  trans = Magma_minproductNoTrans then
            ldda must be at least  max( 1, n ), otherwise  ldda must be at
            least  max( 1, k ).
    
    @param[in]
    beta    REAL.
            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
            supplied as zero then dC need not be set on input.
    
    @param[in,out]
    dC      COMPLEX array of DIMENSION ( lddc, n ).
            Before entry with uplo = 'U' or 'u', the leading n by n
            upper triangular part of the array C must contain the upper
            triangular part of the hermitian matrix and the strictly
            lower triangular part of C is not referenced. On exit, the
            upper triangular part of the array C is overwritten by the
            upper triangular part of the updated matrix.
            Before entry with uplo = 'L' or 'l', the leading n by n
            lower triangular part of the array C must contain the lower
            triangular part of the hermitian matrix and the strictly
            upper triangular part of C is not referenced. On exit, the
            lower triangular part of the array C is overwritten by the
            lower triangular part of the updated matrix.
            Note that the imaginary parts of the diagonal elements need
            not be set, they are assumed to be zero, and on exit they
            are set to zero.

    @param[in]
    lddc    INTEGER.
            On entry, lddc specifies the first dimension of dC as declared
            in  the  calling  (sub)  program.   lddc  must  be  at  least
            max( 1, m ).
    @ingroup magma_minproduct_cblas3
    ********************************************************************/
extern "C" void
magma_minproductblas_cherk_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloatComplex const * const * dA_array, magma_minproduct_int_t ldda,
    float beta,
    magma_minproductFloatComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue )
{

    if( k <= 32 ) {
        magma_minproductblas_cherk_batched_k32(
                  uplo, trans, n, k,
                  alpha, dA_array, ldda,
                  beta,  dC_array, lddc,
                  batchCount, queue );
    }
    else{
        magma_minproductblas_cherk_batched_lg(
                  uplo, trans, n, k,
                  alpha, dA_array, ldda,
                  beta,  dC_array, lddc,
                  batchCount, queue );
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
