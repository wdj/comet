/*
    -- MAGMA_minproduct (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date
*/

#ifndef MAGMA_minproduct_H
#define MAGMA_minproduct_H

/* ------------------------------------------------------------
 * MAGMA_minproduct BLAS Functions
 * --------------------------------------------------------- */
#include "magma_minproductblas.h"
#include "magma_minproduct_batched.h"

/* ------------------------------------------------------------
 * MAGMA_minproduct functions
 * --------------------------------------------------------- */
#include "magma_minproduct_z.h"
#include "magma_minproduct_c.h"
#include "magma_minproduct_d.h"
#include "magma_minproduct_s.h"
#include "magma_minproduct_zc.h"
#include "magma_minproduct_ds.h"
#include "auxiliary.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct Auxiliary functions to get the NB used
*/
magma_minproduct_int_t magma_minproduct_get_smlsize_divideconquer();

// ========================================
// initialization
magma_minproduct_int_t
magma_minproduct_init( void );

magma_minproduct_int_t
magma_minproduct_finalize( void );

void magma_minproduct_version( magma_minproduct_int_t* major, magma_minproduct_int_t* minor, magma_minproduct_int_t* micro );


// ========================================
// memory allocation
magma_minproduct_int_t
magma_minproduct_malloc( magma_minproduct_ptr *ptrPtr, size_t bytes );

magma_minproduct_int_t
magma_minproduct_malloc_cpu( void **ptrPtr, size_t bytes );

magma_minproduct_int_t
magma_minproduct_malloc_pinned( void **ptrPtr, size_t bytes );

magma_minproduct_int_t
magma_minproduct_free_cpu( void *ptr );

#define magma_minproduct_free( ptr ) \
        magma_minproduct_free_internal( ptr, __func__, __FILE__, __LINE__ )

#define magma_minproduct_free_pinned( ptr ) \
        magma_minproduct_free_pinned_internal( ptr, __func__, __FILE__, __LINE__ )

magma_minproduct_int_t
magma_minproduct_free_internal(
    magma_minproduct_ptr ptr,
    const char* func, const char* file, int line );

magma_minproduct_int_t
magma_minproduct_free_pinned_internal(
    void *ptr,
    const char* func, const char* file, int line );


// type-safe convenience functions to avoid using (void**) cast and sizeof(...)
// here n is the number of elements (floats, doubles, etc.) not the number of bytes.
static inline magma_minproduct_int_t magma_minproduct_imalloc( magma_minproductInt_ptr           *ptrPtr, size_t n ) { return magma_minproduct_malloc( (magma_minproduct_ptr*) ptrPtr, n*sizeof(magma_minproduct_int_t)        ); }
static inline magma_minproduct_int_t magma_minproduct_index_malloc( magma_minproductIndex_ptr    *ptrPtr, size_t n ) { return magma_minproduct_malloc( (magma_minproduct_ptr*) ptrPtr, n*sizeof(magma_minproduct_index_t)      ); }
static inline magma_minproduct_int_t magma_minproduct_smalloc( magma_minproductFloat_ptr         *ptrPtr, size_t n ) { return magma_minproduct_malloc( (magma_minproduct_ptr*) ptrPtr, n*sizeof(float)              ); }
static inline magma_minproduct_int_t magma_minproduct_dmalloc( magma_minproductDouble_ptr        *ptrPtr, size_t n ) { return magma_minproduct_malloc( (magma_minproduct_ptr*) ptrPtr, n*sizeof(double)             ); }
static inline magma_minproduct_int_t magma_minproduct_cmalloc( magma_minproductFloatComplex_ptr  *ptrPtr, size_t n ) { return magma_minproduct_malloc( (magma_minproduct_ptr*) ptrPtr, n*sizeof(magma_minproductFloatComplex)  ); }
static inline magma_minproduct_int_t magma_minproduct_zmalloc( magma_minproductDoubleComplex_ptr *ptrPtr, size_t n ) { return magma_minproduct_malloc( (magma_minproduct_ptr*) ptrPtr, n*sizeof(magma_minproductDoubleComplex) ); }

static inline magma_minproduct_int_t magma_minproduct_imalloc_cpu( magma_minproduct_int_t        **ptrPtr, size_t n ) { return magma_minproduct_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_minproduct_int_t)        ); }
static inline magma_minproduct_int_t magma_minproduct_index_malloc_cpu( magma_minproduct_index_t **ptrPtr, size_t n ) { return magma_minproduct_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_minproduct_index_t)      ); }
static inline magma_minproduct_int_t magma_minproduct_smalloc_cpu( float              **ptrPtr, size_t n ) { return magma_minproduct_malloc_cpu( (void**) ptrPtr, n*sizeof(float)              ); }
static inline magma_minproduct_int_t magma_minproduct_dmalloc_cpu( double             **ptrPtr, size_t n ) { return magma_minproduct_malloc_cpu( (void**) ptrPtr, n*sizeof(double)             ); }
static inline magma_minproduct_int_t magma_minproduct_cmalloc_cpu( magma_minproductFloatComplex  **ptrPtr, size_t n ) { return magma_minproduct_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_minproductFloatComplex)  ); }
static inline magma_minproduct_int_t magma_minproduct_zmalloc_cpu( magma_minproductDoubleComplex **ptrPtr, size_t n ) { return magma_minproduct_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_minproductDoubleComplex) ); }

static inline magma_minproduct_int_t magma_minproduct_imalloc_pinned( magma_minproduct_int_t        **ptrPtr, size_t n ) { return magma_minproduct_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_minproduct_int_t)        ); }
static inline magma_minproduct_int_t magma_minproduct_index_malloc_pinned( magma_minproduct_index_t **ptrPtr, size_t n ) { return magma_minproduct_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_minproduct_index_t)      ); }
static inline magma_minproduct_int_t magma_minproduct_smalloc_pinned( float              **ptrPtr, size_t n ) { return magma_minproduct_malloc_pinned( (void**) ptrPtr, n*sizeof(float)              ); }
static inline magma_minproduct_int_t magma_minproduct_dmalloc_pinned( double             **ptrPtr, size_t n ) { return magma_minproduct_malloc_pinned( (void**) ptrPtr, n*sizeof(double)             ); }
static inline magma_minproduct_int_t magma_minproduct_cmalloc_pinned( magma_minproductFloatComplex  **ptrPtr, size_t n ) { return magma_minproduct_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_minproductFloatComplex)  ); }
static inline magma_minproduct_int_t magma_minproduct_zmalloc_pinned( magma_minproductDoubleComplex **ptrPtr, size_t n ) { return magma_minproduct_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_minproductDoubleComplex) ); }


// ========================================
// device support
magma_minproduct_int_t magma_minproduct_getdevice_arch();

void magma_minproduct_getdevices(
    magma_minproduct_device_t* devices,
    magma_minproduct_int_t     size,
    magma_minproduct_int_t*    numPtr );

void magma_minproduct_getdevice( magma_minproduct_device_t* dev );

void magma_minproduct_setdevice( magma_minproduct_device_t dev );

void magma_minproduct_device_sync();


// ========================================
// queue support
#define magma_minproduct_queue_create( /*device,*/ queuePtr ) \
        magma_minproduct_queue_create_internal( queuePtr, __func__, __FILE__, __LINE__ )

#define magma_minproduct_queue_destroy( queue ) \
        magma_minproduct_queue_destroy_internal( queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_queue_sync( queue ) \
        magma_minproduct_queue_sync_internal( queue, __func__, __FILE__, __LINE__ )

void magma_minproduct_queue_create_internal(
    /*magma_minproduct_device_t device,*/ magma_minproduct_queue_t* queuePtr,
    const char* func, const char* file, int line );

void magma_minproduct_queue_destroy_internal(
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void magma_minproduct_queue_sync_internal(
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

/// Currently, magma_minproduct_queue_t == cudaStream_t.
/// Almost certainly this will change in the future,
/// so these get & set the associated stream in a forward-compatible manner.
/// @see magma_minproduct_queue_set_cuda_stream
#define magma_minproduct_queue_get_cuda_stream( queue ) (queue)

/// @see magma_minproduct_queue_get_cuda_stream
#define magma_minproduct_queue_set_cuda_stream( queue, stream ) ((queue) = (stream))


// ========================================
// event support
void magma_minproduct_event_create( magma_minproduct_event_t* eventPtr );

void magma_minproduct_event_destroy( magma_minproduct_event_t event );

void magma_minproduct_event_record( magma_minproduct_event_t event, magma_minproduct_queue_t queue );

void magma_minproduct_event_query( magma_minproduct_event_t event );

// blocks CPU until event occurs
void magma_minproduct_event_sync( magma_minproduct_event_t event );

// blocks queue (but not CPU) until event occurs
void magma_minproduct_queue_wait_event( magma_minproduct_queue_t queue, magma_minproduct_event_t event );


// ========================================
// error handler
void magma_minproduct_xerbla( const char *name, magma_minproduct_int_t info );

const char* magma_minproduct_strerror( magma_minproduct_int_t error );


// ========================================
/// For integers x >= 0, y > 0, returns ceil( x/y ).
/// For x == 0, this is 0.
__host__ __device__
static inline magma_minproduct_int_t magma_minproduct_ceildiv( magma_minproduct_int_t x, magma_minproduct_int_t y )
{
    return (x + y - 1)/y;
}

/// For integers x >= 0, y > 0, returns x rounded up to multiple of y.
/// For x == 0, this is 0.
/// This implementation does not assume y is a power of 2.
__host__ __device__
static inline magma_minproduct_int_t magma_minproduct_roundup( magma_minproduct_int_t x, magma_minproduct_int_t y )
{
    return magma_minproduct_ceildiv( x, y ) * y;
}


// ========================================
// real and complex square root
// sqrt alone cannot be caught by the generation script because of tsqrt
static inline float  magma_minproduct_ssqrt( float  x ) { return sqrtf( x ); }
static inline double magma_minproduct_dsqrt( double x ) { return sqrt( x ); }
magma_minproductFloatComplex    magma_minproduct_csqrt( magma_minproductFloatComplex  x );
magma_minproductDoubleComplex   magma_minproduct_zsqrt( magma_minproductDoubleComplex x );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproduct_H */
