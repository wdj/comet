/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mark Gates
       @generated from zset_get.cpp normal z -> c, Fri Jan 30 19:00:20 2015
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma_minproduct.h"
#include "error.h"

#ifdef HAVE_CUBLAS

// ========================================
// copying vectors
extern "C" void
magma_minproduct_csetvector_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex const* hx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr    dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, sizeof(magma_minproductFloatComplex),
        hx_src, incx,
        dy_dst, incy );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_minproduct_csetvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex const* hx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr    dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, sizeof(magma_minproductFloatComplex),
        hx_src, incx,
        dy_dst, incy, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_minproduct_cgetvector_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex*       hy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, sizeof(magma_minproductFloatComplex),
        dx_src, incx,
        hy_dst, incy );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_minproduct_cgetvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex*       hy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, sizeof(magma_minproductFloatComplex),
        dx_src, incx,
        hy_dst, incy, queue );
    check_xerror( status, func, file, line );
}

// --------------------
// TODO compare performance with cublasCcopy BLAS function.
// But this implementation can handle any element size, not just [sdcz] precisions.
extern "C" void
magma_minproduct_ccopyvector_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr    dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line )
{
    if ( incx == 1 && incy == 1 ) {
        cudaError_t status;
        status = cudaMemcpy(
            dy_dst,
            dx_src,
            n*sizeof(magma_minproductFloatComplex), cudaMemcpyDeviceToDevice );
        check_xerror( status, func, file, line );
    }
    else {
        magma_minproduct_ccopymatrix_internal(
            1, n, dx_src, incx, dy_dst, incy, func, file, line );
    }
}

// --------------------
extern "C" void
magma_minproduct_ccopyvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr    dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line )
{
    if ( incx == 1 && incy == 1 ) {
        cudaError_t status;
        status = cudaMemcpyAsync(
            dy_dst,
            dx_src,
            n*sizeof(magma_minproductFloatComplex), cudaMemcpyDeviceToDevice, queue );
        check_xerror( status, func, file, line );
    }
    else {
        magma_minproduct_ccopymatrix_async_internal(
            1, n, dx_src, incx, dy_dst, incy, queue, func, file, line );
    }
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C" void
magma_minproduct_csetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex const* hA_src, magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr    dB_dst, magma_minproduct_int_t ldb,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, sizeof(magma_minproductFloatComplex),
        hA_src, lda,
        dB_dst, ldb );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_minproduct_csetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex const* hA_src, magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr    dB_dst, magma_minproduct_int_t ldb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, sizeof(magma_minproductFloatComplex),
        hA_src, lda,
        dB_dst, ldb, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_minproduct_cgetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA_src, magma_minproduct_int_t lda,
    magma_minproductFloatComplex*          hB_dst, magma_minproduct_int_t ldb,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, sizeof(magma_minproductFloatComplex),
        dA_src, lda,
        hB_dst, ldb );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_minproduct_cgetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA_src, magma_minproduct_int_t lda,
    magma_minproductFloatComplex*          hB_dst, magma_minproduct_int_t ldb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, sizeof(magma_minproductFloatComplex),
        dA_src, lda,
        hB_dst, ldb, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_minproduct_ccopymatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA_src, magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr       dB_dst, magma_minproduct_int_t ldb,
    const char* func, const char* file, int line )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, ldb*sizeof(magma_minproductFloatComplex),
        dA_src, lda*sizeof(magma_minproductFloatComplex),
        m*sizeof(magma_minproductFloatComplex), n, cudaMemcpyDeviceToDevice );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_minproduct_ccopymatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA_src, magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr       dB_dst, magma_minproduct_int_t ldb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, ldb*sizeof(magma_minproductFloatComplex),
        dA_src, lda*sizeof(magma_minproductFloatComplex),
        m*sizeof(magma_minproductFloatComplex), n, cudaMemcpyDeviceToDevice, queue );
    check_xerror( status, func, file, line );
}

#endif // HAVE_CUBLAS
