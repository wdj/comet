/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mark Gates
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma_tally2.h"
#include "error.h"

#ifdef HAVE_CUBLAS

// ========================================
// memory allocation
// Allocate size bytes on GPU, returning pointer in ptrPtr.
extern "C" magma_tally2_int_t
magma_tally2_malloc( magma_tally2_ptr* ptrPtr, size_t size )
{
    // CUDA can't allocate 0 bytes, so allocate some minimal size
    if ( size == 0 )
        size = sizeof(magma_tally2DoubleComplex);
    if ( cudaSuccess != cudaMalloc( ptrPtr, size )) {
        return MAGMA_tally2_ERR_DEVICE_ALLOC;
    }
    return MAGMA_tally2_SUCCESS;
}

// --------------------
// Free GPU memory allocated by magma_tally2_malloc.
extern "C" magma_tally2_int_t
magma_tally2_free_internal( magma_tally2_ptr ptr,
    const char* func, const char* file, int line )
{
    cudaError_t err = cudaFree( ptr );
    check_xerror( err, func, file, line );
    if ( err != cudaSuccess ) {
        return MAGMA_tally2_ERR_INVALID_PTR;
    }
    return MAGMA_tally2_SUCCESS;
}

// --------------------
// Allocate size bytes on CPU, returning pointer in ptrPtr.
// The purpose of using this instead of malloc() is to properly align arrays
// for vector (SSE) instructions. The default implementation uses
// posix_memalign (on Linux, MacOS, etc.) or _aligned_malloc (on Windows)
// to align memory to a 32 byte boundary.
// Use magma_tally2_free_cpu() to free this memory.
extern "C" magma_tally2_int_t
magma_tally2_malloc_cpu( void** ptrPtr, size_t size )
{
    // malloc and free sometimes don't work for size=0, so allocate some minimal size
    if ( size == 0 )
        size = sizeof(magma_tally2DoubleComplex);
#if 1
#if defined( _WIN32 ) || defined( _WIN64 )
    *ptrPtr = _aligned_malloc( size, 32 );
    if ( *ptrPtr == NULL ) {
        return MAGMA_tally2_ERR_HOST_ALLOC;
    }
#else
    int err = posix_memalign( ptrPtr, 32, size );
    if ( err != 0 ) {
        *ptrPtr = NULL;
        return MAGMA_tally2_ERR_HOST_ALLOC;
    }
#endif
#else
    *ptrPtr = malloc( size );
    if ( *ptrPtr == NULL ) {
        return MAGMA_tally2_ERR_HOST_ALLOC;
    }
#endif
    return MAGMA_tally2_SUCCESS;
}

// --------------------
// Free CPU pinned memory previously allocated by magma_tally2_malloc_pinned.
// The default implementation uses free(), which works for both malloc and posix_memalign.
// For Windows, _aligned_free() is used.
extern "C" magma_tally2_int_t
magma_tally2_free_cpu( void* ptr )
{
#if defined( _WIN32 ) || defined( _WIN64 )
    _aligned_free( ptr );
#else
    free( ptr );
#endif
    return MAGMA_tally2_SUCCESS;
}

// --------------------
// Allocate size bytes on CPU in pinned memory, returning pointer in ptrPtr.
extern "C" magma_tally2_int_t
magma_tally2_malloc_pinned( void** ptrPtr, size_t size )
{
    // CUDA can't allocate 0 bytes, so allocate some minimal size
    // (for pinned memory, the error is detected in free)
    if ( size == 0 )
        size = sizeof(magma_tally2DoubleComplex);
    if ( cudaSuccess != cudaMallocHost( ptrPtr, size )) {
        return MAGMA_tally2_ERR_HOST_ALLOC;
    }
    return MAGMA_tally2_SUCCESS;
}

// --------------------
// Free CPU pinned memory previously allocated by magma_tally2_malloc_pinned.
extern "C" magma_tally2_int_t
magma_tally2_free_pinned_internal( void* ptr,
    const char* func, const char* file, int line )
{
    cudaError_t err = cudaFreeHost( ptr );
    check_xerror( err, func, file, line );
    if ( cudaSuccess != err ) {
        return MAGMA_tally2_ERR_INVALID_PTR;
    }
    return MAGMA_tally2_SUCCESS;
}

#endif // HAVE_CUBLAS
