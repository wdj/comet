/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zherk_fermi_batched_k32.cu normal z -> c, Fri Jan 30 19:00:10 2015

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
#include "common_magma_tally2.h"
#include "commonblas_c.h"

#define PRECISION_c

///////////////////////////////////////////////////////////////////////////////////////////////////

#include "cgemm_fermi_kernels_batched_k32.h"

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
            k  when  trans = Magma_tally2NoTrans,  and is  n  otherwise.
            Before entry with  trans = Magma_tally2NoTrans,  the leading  m by k
            part of the array dA must contain the matrix dA, otherwise
            the leading  k by m  part of the array dA must contain  the
            matrix dA.
    
    @param[in]
    ldda    INTEGER.
            On entry, ldda specifies the first dimension of A as declared
            in the calling (sub) program. When  trans = Magma_tally2NoTrans then
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
    @ingroup magma_tally2_cblas3
    ********************************************************************/
extern "C" void
magma_tally2blas_cherk_batched_k32(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2FloatComplex const * const * dA_array, magma_tally2_int_t ldda,
    float beta,
    magma_tally2FloatComplex **dC_array, magma_tally2_int_t lddc, magma_tally2_int_t batchCount, magma_tally2_queue_t queue )
{
    magma_tally2FloatComplex cbeta  = MAGMA_tally2_C_MAKE( beta, 0. );
    magma_tally2FloatComplex calpha = MAGMA_tally2_C_MAKE( alpha, 0. );

    magma_tally2_int_t info = 0;
    if      ( uplo != Magma_tally2Upper && uplo != Magma_tally2Lower )
        info = -1;
    else if ( trans != Magma_tally2NoTrans && trans != Magma_tally2Trans && trans != Magma_tally2ConjTrans )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( k < 0 )
        info = -4;
    else if ( trans == Magma_tally2NoTrans ? ldda < n : ldda < k )
        info = -7;
    else if ( lddc < n )
        info = -10;

    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally2_int_t arch = magma_tally2_getdevice_arch();
    if ( arch < 200  ) {
        printf("not supported \n"); // TODO call cublas
        return;
    }
    
    // --------------------
    // CUDA ARCH 2.x (Fermi) version
    if ( n <= 0 || k <= 0 )
        return;
    
    size_t offsetA = 0;
    int TransA = 0, TransB = 0, uploA = 0;

    if      ( uplo == Magma_tally2Lower )
        uploA = 1;
    else if ( uplo == Magma_tally2Upper )
        uploA = 2;

    if      ( trans == Magma_tally2NoTrans )
        #if defined(PRECISION_z) || defined(PRECISION_c)     
        TransB = 2;
        #else
        TransB = 1;
        #endif
    else if ( trans == Magma_tally2Trans || trans == Magma_tally2ConjTrans)
        #if defined(PRECISION_z) || defined(PRECISION_c)     
        TransA = 2;
        #else
        TransA = 1;
        #endif


    #ifdef TEXTURE_1D
        size_t sizeA = (size_t) ldda * (size_t) (!TransA ? k : n);

        size_t CUBLAS_MAX_1DBUF_SIZE = ((1 << 27) - 512);
        if ( sizeA >= CUBLAS_MAX_1DBUF_SIZE )
        {
            printf("not supported \n"); // TODO call cublas
            return;
        }
        // Set textures parameters
        tex_ref_A.normalized = false;
        tex_ref_A.filterMode = cudaFilterModePoint;
        tex_ref_A.addressMode[0] = cudaAddressModeClamp;

        // Bind A and B to texture references
        cudaError_t err;
        err = cudaBindTexture(&offsetA, tex_ref_A, dA_array[0], sizeA*sizeof(magma_tally2FloatComplex));
        if ( err != cudaSuccess ) {
            fprintf( stderr, "cannot bind A to texture: %s (%d)\n", cudaGetErrorString(err), err );
            return;
        }
    #endif

    // Set up grids
    dim3 dimBlock(DIM_X, DIM_Y);

    offsetA = offsetA/sizeof(magma_tally2FloatComplex);
 
    if ( TransA == 0 && TransB == 1 ) {
        dim3 dimGrid( (n - 1)/BLK_M_nt + 1,
                      (n - 1)/BLK_N_nt + 1 ,
                      batchCount );
        magma_tally2blas_c_herk_kernel_fermi_nt_batched<<< dimGrid, dimBlock, 0, queue >>>(
            uploA, n, k, dA_array, ldda, dA_array, ldda, dC_array, lddc, calpha, cbeta,
            (int)offsetA, (int)offsetA );
    }
    else if ( TransA == 0 && TransB == 2 ) {
        dim3 dimGrid( (n - 1)/BLK_M_nc + 1,
                      (n - 1)/BLK_N_nc + 1 ,
                      batchCount );
         magma_tally2blas_c_herk_kernel_fermi_nc_batched<<< dimGrid, dimBlock, 0, queue >>>(
            uploA, n, k, dA_array, ldda, dA_array, ldda, dC_array, lddc, calpha, cbeta,
            (int)offsetA, (int)offsetA );
    }
    else if ( TransA == 1 && TransB == 0 ) {
        dim3 dimGrid( (n - 1)/BLK_M_tn + 1,
                      (n - 1)/BLK_N_tn + 1 ,
                      batchCount );
         magma_tally2blas_c_herk_kernel_fermi_tn_batched<<< dimGrid, dimBlock, 0, queue >>>(
            uploA, n, k, dA_array, ldda, dA_array, ldda, dC_array, lddc, calpha, cbeta,
            (int)offsetA, (int)offsetA );
    }
    else if ( TransA == 2 && TransB == 0 ) {
        dim3 dimGrid( (n - 1)/BLK_M_cn + 1,
                      (n - 1)/BLK_N_cn + 1 ,
                      batchCount );
         magma_tally2blas_c_herk_kernel_fermi_cn_batched<<< dimGrid, dimBlock, 0, queue >>>(
            uploA, n, k, dA_array, ldda, dA_array, ldda, dC_array, lddc, calpha, cbeta,
            (int)offsetA, (int)offsetA );
    }

    #ifdef TEXTURE_1D
        cudaUnbindTexture( tex_ref_A );
    #endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////
