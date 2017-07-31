/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mark Gates
       @precisions normal z -> s d c
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma_tally2.h"
#include "error.h"

#ifdef HAVE_CUBLAS

// generic, type-independent routines to copy data.
// type-safe versions which avoid the user needing sizeof(...) are in [sdcz]set_get.cpp

// ========================================
// copying vectors
extern "C" void
magma_tally2_setvector_internal(
    magma_tally2_int_t n, magma_tally2_int_t elemSize,
    void const* hx_src, magma_tally2_int_t incx,
    magma_tally2_ptr   dy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, elemSize,
        hx_src, incx,
        dy_dst, incy );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_setvector_async_internal(
    magma_tally2_int_t n, magma_tally2_int_t elemSize,
    void const* hx_src, magma_tally2_int_t incx,
    magma_tally2_ptr   dy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, elemSize,
        hx_src, incx,
        dy_dst, incy, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_getvector_internal(
    magma_tally2_int_t n, magma_tally2_int_t elemSize,
    magma_tally2_const_ptr dx_src, magma_tally2_int_t incx,
    void*           hy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, elemSize,
        dx_src, incx,
        hy_dst, incy );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_getvector_async_internal(
    magma_tally2_int_t n, magma_tally2_int_t elemSize,
    magma_tally2_const_ptr dx_src, magma_tally2_int_t incx,
    void*           hy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, elemSize,
        dx_src, incx,
        hy_dst, incy, queue );
    check_xerror( status, func, file, line );
}

// --------------------
// TODO compare performance with cublasZcopy BLAS function.
// But this implementation can handle any element size, not just [sdcz] precisions.
extern "C" void
magma_tally2_copyvector_internal(
    magma_tally2_int_t n, magma_tally2_int_t elemSize,
    magma_tally2_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2_ptr       dy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line )
{
    if ( incx == 1 && incy == 1 ) {
        cudaError_t status;
        status = cudaMemcpy(
            dy_dst,
            dx_src,
            n*elemSize, cudaMemcpyDeviceToDevice );
        check_xerror( status, func, file, line );
    }
    else {
        magma_tally2_copymatrix_internal(
            1, n, elemSize, dx_src, incx, dy_dst, incy, func, file, line );
    }
}

// --------------------
extern "C" void
magma_tally2_copyvector_async_internal(
    magma_tally2_int_t n, magma_tally2_int_t elemSize,
    magma_tally2_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2_ptr       dy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    if ( incx == 1 && incy == 1 ) {
        cudaError_t status;
        status = cudaMemcpyAsync(
            dy_dst,
            dx_src,
            n*elemSize, cudaMemcpyDeviceToDevice, queue );
        check_xerror( status, func, file, line );
    }
    else {
        magma_tally2_copymatrix_async_internal(
            1, n, elemSize, dx_src, incx, dy_dst, incy, queue, func, file, line );
    }
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C" void
magma_tally2_setmatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t elemSize,
    void const* hA_src, magma_tally2_int_t ldha,
    magma_tally2_ptr   dB_dst, magma_tally2_int_t lddb,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, elemSize,
        hA_src, ldha,
        dB_dst, lddb );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_setmatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t elemSize,
    void const* hA_src, magma_tally2_int_t ldha,
    magma_tally2_ptr   dB_dst, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, elemSize,
        hA_src, ldha,
        dB_dst, lddb, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_getmatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t elemSize,
    magma_tally2_const_ptr dA_src, magma_tally2_int_t ldda,
    void*           hB_dst, magma_tally2_int_t ldhb,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, elemSize,
        dA_src, ldda,
        hB_dst, ldhb );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_getmatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t elemSize,
    magma_tally2_const_ptr dA_src, magma_tally2_int_t ldda,
    void*           hB_dst, magma_tally2_int_t ldhb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, elemSize,
        dA_src, ldda,
        hB_dst, ldhb, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_copymatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t elemSize,
    magma_tally2_const_ptr dA_src, magma_tally2_int_t ldda,
    magma_tally2_ptr       dB_dst, magma_tally2_int_t lddb,
    const char* func, const char* file, int line )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, lddb*elemSize,
        dA_src, ldda*elemSize,
        m*elemSize, n, cudaMemcpyDeviceToDevice );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_copymatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t elemSize,
    magma_tally2_const_ptr dA_src, magma_tally2_int_t ldda,
    magma_tally2_ptr       dB_dst, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, lddb*elemSize,
        dA_src, ldda*elemSize,
        m*elemSize, n, cudaMemcpyDeviceToDevice, queue );
    check_xerror( status, func, file, line );
}

#endif // HAVE_CUBLAS
