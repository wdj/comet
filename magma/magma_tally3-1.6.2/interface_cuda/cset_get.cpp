/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mark Gates
       @generated from zset_get.cpp normal z -> c, Fri Jan 30 19:00:20 2015
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma_tally3.h"
#include "error.h"

#ifdef HAVE_CUBLAS

// ========================================
// copying vectors
extern "C" void
magma_tally3_csetvector_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex const* hx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr    dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, sizeof(magma_tally3FloatComplex),
        hx_src, incx,
        dy_dst, incy );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally3_csetvector_async_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex const* hx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr    dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, sizeof(magma_tally3FloatComplex),
        hx_src, incx,
        dy_dst, incy, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally3_cgetvector_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex*       hy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, sizeof(magma_tally3FloatComplex),
        dx_src, incx,
        hy_dst, incy );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally3_cgetvector_async_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex*       hy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, sizeof(magma_tally3FloatComplex),
        dx_src, incx,
        hy_dst, incy, queue );
    check_xerror( status, func, file, line );
}

// --------------------
// TODO compare performance with cublasCcopy BLAS function.
// But this implementation can handle any element size, not just [sdcz] precisions.
extern "C" void
magma_tally3_ccopyvector_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr    dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{
    if ( incx == 1 && incy == 1 ) {
        cudaError_t status;
        status = cudaMemcpy(
            dy_dst,
            dx_src,
            n*sizeof(magma_tally3FloatComplex), cudaMemcpyDeviceToDevice );
        check_xerror( status, func, file, line );
    }
    else {
        magma_tally3_ccopymatrix_internal(
            1, n, dx_src, incx, dy_dst, incy, func, file, line );
    }
}

// --------------------
extern "C" void
magma_tally3_ccopyvector_async_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr    dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line )
{
    if ( incx == 1 && incy == 1 ) {
        cudaError_t status;
        status = cudaMemcpyAsync(
            dy_dst,
            dx_src,
            n*sizeof(magma_tally3FloatComplex), cudaMemcpyDeviceToDevice, queue );
        check_xerror( status, func, file, line );
    }
    else {
        magma_tally3_ccopymatrix_async_internal(
            1, n, dx_src, incx, dy_dst, incy, queue, func, file, line );
    }
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C" void
magma_tally3_csetmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex const* hA_src, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr    dB_dst, magma_tally3_int_t ldb,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, sizeof(magma_tally3FloatComplex),
        hA_src, lda,
        dB_dst, ldb );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally3_csetmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex const* hA_src, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr    dB_dst, magma_tally3_int_t ldb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, sizeof(magma_tally3FloatComplex),
        hA_src, lda,
        dB_dst, ldb, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally3_cgetmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA_src, magma_tally3_int_t lda,
    magma_tally3FloatComplex*          hB_dst, magma_tally3_int_t ldb,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, sizeof(magma_tally3FloatComplex),
        dA_src, lda,
        hB_dst, ldb );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally3_cgetmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA_src, magma_tally3_int_t lda,
    magma_tally3FloatComplex*          hB_dst, magma_tally3_int_t ldb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, sizeof(magma_tally3FloatComplex),
        dA_src, lda,
        hB_dst, ldb, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally3_ccopymatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA_src, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr       dB_dst, magma_tally3_int_t ldb,
    const char* func, const char* file, int line )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, ldb*sizeof(magma_tally3FloatComplex),
        dA_src, lda*sizeof(magma_tally3FloatComplex),
        m*sizeof(magma_tally3FloatComplex), n, cudaMemcpyDeviceToDevice );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally3_ccopymatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA_src, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr       dB_dst, magma_tally3_int_t ldb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, ldb*sizeof(magma_tally3FloatComplex),
        dA_src, lda*sizeof(magma_tally3FloatComplex),
        m*sizeof(magma_tally3FloatComplex), n, cudaMemcpyDeviceToDevice, queue );
    check_xerror( status, func, file, line );
}

#endif // HAVE_CUBLAS
