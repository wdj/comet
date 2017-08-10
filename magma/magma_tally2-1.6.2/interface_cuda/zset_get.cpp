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

// ========================================
// copying vectors
extern "C" void
magma_tally2_zsetvector_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex const* hx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr    dy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetVector(
        n, sizeof(magma_tally2DoubleComplex),
        hx_src, incx,
        dy_dst, incy );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_zsetvector_async_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex const* hx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr    dy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetVectorAsync(
        n, sizeof(magma_tally2DoubleComplex),
        hx_src, incx,
        dy_dst, incy, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_zgetvector_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex*       hy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetVector(
        n, sizeof(magma_tally2DoubleComplex),
        dx_src, incx,
        hy_dst, incy );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_zgetvector_async_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex*       hy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetVectorAsync(
        n, sizeof(magma_tally2DoubleComplex),
        dx_src, incx,
        hy_dst, incy, queue );
    check_xerror( status, func, file, line );
}

// --------------------
// TODO compare performance with cublasZcopy BLAS function.
// But this implementation can handle any element size, not just [sdcz] precisions.
extern "C" void
magma_tally2_zcopyvector_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr    dy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line )
{
    if ( incx == 1 && incy == 1 ) {
        cudaError_t status;
        status = cudaMemcpy(
            dy_dst,
            dx_src,
            n*sizeof(magma_tally2DoubleComplex), cudaMemcpyDeviceToDevice );
        check_xerror( status, func, file, line );
    }
    else {
        magma_tally2_zcopymatrix_internal(
            1, n, dx_src, incx, dy_dst, incy, func, file, line );
    }
}

// --------------------
extern "C" void
magma_tally2_zcopyvector_async_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr    dy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    if ( incx == 1 && incy == 1 ) {
        cudaError_t status;
        status = cudaMemcpyAsync(
            dy_dst,
            dx_src,
            n*sizeof(magma_tally2DoubleComplex), cudaMemcpyDeviceToDevice, queue );
        check_xerror( status, func, file, line );
    }
    else {
        magma_tally2_zcopymatrix_async_internal(
            1, n, dx_src, incx, dy_dst, incy, queue, func, file, line );
    }
}


// ========================================
// copying sub-matrices (contiguous columns)
extern "C" void
magma_tally2_zsetmatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex const* hA_src, magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr    dB_dst, magma_tally2_int_t ldb,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetMatrix(
        m, n, sizeof(magma_tally2DoubleComplex),
        hA_src, lda,
        dB_dst, ldb );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_zsetmatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex const* hA_src, magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr    dB_dst, magma_tally2_int_t ldb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasSetMatrixAsync(
        m, n, sizeof(magma_tally2DoubleComplex),
        hA_src, lda,
        dB_dst, ldb, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_zgetmatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA_src, magma_tally2_int_t lda,
    magma_tally2DoubleComplex*          hB_dst, magma_tally2_int_t ldb,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetMatrix(
        m, n, sizeof(magma_tally2DoubleComplex),
        dA_src, lda,
        hB_dst, ldb );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_zgetmatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA_src, magma_tally2_int_t lda,
    magma_tally2DoubleComplex*          hB_dst, magma_tally2_int_t ldb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cublasStatus_t status;
    status = cublasGetMatrixAsync(
        m, n, sizeof(magma_tally2DoubleComplex),
        dA_src, lda,
        hB_dst, ldb, queue );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_zcopymatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA_src, magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr       dB_dst, magma_tally2_int_t ldb,
    const char* func, const char* file, int line )
{
    cudaError_t status;
    status = cudaMemcpy2D(
        dB_dst, ldb*sizeof(magma_tally2DoubleComplex),
        dA_src, lda*sizeof(magma_tally2DoubleComplex),
        m*sizeof(magma_tally2DoubleComplex), n, cudaMemcpyDeviceToDevice );
    check_xerror( status, func, file, line );
}

// --------------------
extern "C" void
magma_tally2_zcopymatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA_src, magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr       dB_dst, magma_tally2_int_t ldb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    cudaError_t status;
    status = cudaMemcpy2DAsync(
        dB_dst, ldb*sizeof(magma_tally2DoubleComplex),
        dA_src, lda*sizeof(magma_tally2DoubleComplex),
        m*sizeof(magma_tally2DoubleComplex), n, cudaMemcpyDeviceToDevice, queue );
    check_xerror( status, func, file, line );
}

#endif // HAVE_CUBLAS