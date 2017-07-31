/*
    -- MAGMA_tally2 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date
*/

#ifndef MAGMA_tally2_H
#define MAGMA_tally2_H

/* ------------------------------------------------------------
 * MAGMA_tally2 BLAS Functions
 * --------------------------------------------------------- */
#include "magma_tally2blas.h"
#include "magma_tally2_batched.h"

/* ------------------------------------------------------------
 * MAGMA_tally2 functions
 * --------------------------------------------------------- */
#include "magma_tally2_z.h"
#include "magma_tally2_c.h"
#include "magma_tally2_d.h"
#include "magma_tally2_s.h"
#include "magma_tally2_zc.h"
#include "magma_tally2_ds.h"
#include "auxiliary.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 Auxiliary functions to get the NB used
*/
magma_tally2_int_t magma_tally2_get_smlsize_divideconquer();

// ========================================
// initialization
magma_tally2_int_t
magma_tally2_init( void );

magma_tally2_int_t
magma_tally2_finalize( void );

void magma_tally2_version( magma_tally2_int_t* major, magma_tally2_int_t* minor, magma_tally2_int_t* micro );


// ========================================
// memory allocation
magma_tally2_int_t
magma_tally2_malloc( magma_tally2_ptr *ptrPtr, size_t bytes );

magma_tally2_int_t
magma_tally2_malloc_cpu( void **ptrPtr, size_t bytes );

magma_tally2_int_t
magma_tally2_malloc_pinned( void **ptrPtr, size_t bytes );

magma_tally2_int_t
magma_tally2_free_cpu( void *ptr );

#define magma_tally2_free( ptr ) \
        magma_tally2_free_internal( ptr, __func__, __FILE__, __LINE__ )

#define magma_tally2_free_pinned( ptr ) \
        magma_tally2_free_pinned_internal( ptr, __func__, __FILE__, __LINE__ )

magma_tally2_int_t
magma_tally2_free_internal(
    magma_tally2_ptr ptr,
    const char* func, const char* file, int line );

magma_tally2_int_t
magma_tally2_free_pinned_internal(
    void *ptr,
    const char* func, const char* file, int line );


// type-safe convenience functions to avoid using (void**) cast and sizeof(...)
// here n is the number of elements (floats, doubles, etc.) not the number of bytes.
static inline magma_tally2_int_t magma_tally2_imalloc( magma_tally2Int_ptr           *ptrPtr, size_t n ) { return magma_tally2_malloc( (magma_tally2_ptr*) ptrPtr, n*sizeof(magma_tally2_int_t)        ); }
static inline magma_tally2_int_t magma_tally2_index_malloc( magma_tally2Index_ptr    *ptrPtr, size_t n ) { return magma_tally2_malloc( (magma_tally2_ptr*) ptrPtr, n*sizeof(magma_tally2_index_t)      ); }
static inline magma_tally2_int_t magma_tally2_smalloc( magma_tally2Float_ptr         *ptrPtr, size_t n ) { return magma_tally2_malloc( (magma_tally2_ptr*) ptrPtr, n*sizeof(float)              ); }
static inline magma_tally2_int_t magma_tally2_dmalloc( magma_tally2Double_ptr        *ptrPtr, size_t n ) { return magma_tally2_malloc( (magma_tally2_ptr*) ptrPtr, n*sizeof(double)             ); }
static inline magma_tally2_int_t magma_tally2_cmalloc( magma_tally2FloatComplex_ptr  *ptrPtr, size_t n ) { return magma_tally2_malloc( (magma_tally2_ptr*) ptrPtr, n*sizeof(magma_tally2FloatComplex)  ); }
static inline magma_tally2_int_t magma_tally2_zmalloc( magma_tally2DoubleComplex_ptr *ptrPtr, size_t n ) { return magma_tally2_malloc( (magma_tally2_ptr*) ptrPtr, n*sizeof(magma_tally2DoubleComplex) ); }

static inline magma_tally2_int_t magma_tally2_imalloc_cpu( magma_tally2_int_t        **ptrPtr, size_t n ) { return magma_tally2_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_tally2_int_t)        ); }
static inline magma_tally2_int_t magma_tally2_index_malloc_cpu( magma_tally2_index_t **ptrPtr, size_t n ) { return magma_tally2_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_tally2_index_t)      ); }
static inline magma_tally2_int_t magma_tally2_smalloc_cpu( float              **ptrPtr, size_t n ) { return magma_tally2_malloc_cpu( (void**) ptrPtr, n*sizeof(float)              ); }
static inline magma_tally2_int_t magma_tally2_dmalloc_cpu( double             **ptrPtr, size_t n ) { return magma_tally2_malloc_cpu( (void**) ptrPtr, n*sizeof(double)             ); }
static inline magma_tally2_int_t magma_tally2_cmalloc_cpu( magma_tally2FloatComplex  **ptrPtr, size_t n ) { return magma_tally2_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_tally2FloatComplex)  ); }
static inline magma_tally2_int_t magma_tally2_zmalloc_cpu( magma_tally2DoubleComplex **ptrPtr, size_t n ) { return magma_tally2_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_tally2DoubleComplex) ); }

static inline magma_tally2_int_t magma_tally2_imalloc_pinned( magma_tally2_int_t        **ptrPtr, size_t n ) { return magma_tally2_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_tally2_int_t)        ); }
static inline magma_tally2_int_t magma_tally2_index_malloc_pinned( magma_tally2_index_t **ptrPtr, size_t n ) { return magma_tally2_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_tally2_index_t)      ); }
static inline magma_tally2_int_t magma_tally2_smalloc_pinned( float              **ptrPtr, size_t n ) { return magma_tally2_malloc_pinned( (void**) ptrPtr, n*sizeof(float)              ); }
static inline magma_tally2_int_t magma_tally2_dmalloc_pinned( double             **ptrPtr, size_t n ) { return magma_tally2_malloc_pinned( (void**) ptrPtr, n*sizeof(double)             ); }
static inline magma_tally2_int_t magma_tally2_cmalloc_pinned( magma_tally2FloatComplex  **ptrPtr, size_t n ) { return magma_tally2_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_tally2FloatComplex)  ); }
static inline magma_tally2_int_t magma_tally2_zmalloc_pinned( magma_tally2DoubleComplex **ptrPtr, size_t n ) { return magma_tally2_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_tally2DoubleComplex) ); }


// ========================================
// device support
magma_tally2_int_t magma_tally2_getdevice_arch();

void magma_tally2_getdevices(
    magma_tally2_device_t* devices,
    magma_tally2_int_t     size,
    magma_tally2_int_t*    numPtr );

void magma_tally2_getdevice( magma_tally2_device_t* dev );

void magma_tally2_setdevice( magma_tally2_device_t dev );

void magma_tally2_device_sync();


// ========================================
// queue support
#define magma_tally2_queue_create( /*device,*/ queuePtr ) \
        magma_tally2_queue_create_internal( queuePtr, __func__, __FILE__, __LINE__ )

#define magma_tally2_queue_destroy( queue ) \
        magma_tally2_queue_destroy_internal( queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_queue_sync( queue ) \
        magma_tally2_queue_sync_internal( queue, __func__, __FILE__, __LINE__ )

void magma_tally2_queue_create_internal(
    /*magma_tally2_device_t device,*/ magma_tally2_queue_t* queuePtr,
    const char* func, const char* file, int line );

void magma_tally2_queue_destroy_internal(
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void magma_tally2_queue_sync_internal(
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

/// Currently, magma_tally2_queue_t == cudaStream_t.
/// Almost certainly this will change in the future,
/// so these get & set the associated stream in a forward-compatible manner.
/// @see magma_tally2_queue_set_cuda_stream
#define magma_tally2_queue_get_cuda_stream( queue ) (queue)

/// @see magma_tally2_queue_get_cuda_stream
#define magma_tally2_queue_set_cuda_stream( queue, stream ) ((queue) = (stream))


// ========================================
// event support
void magma_tally2_event_create( magma_tally2_event_t* eventPtr );

void magma_tally2_event_destroy( magma_tally2_event_t event );

void magma_tally2_event_record( magma_tally2_event_t event, magma_tally2_queue_t queue );

void magma_tally2_event_query( magma_tally2_event_t event );

// blocks CPU until event occurs
void magma_tally2_event_sync( magma_tally2_event_t event );

// blocks queue (but not CPU) until event occurs
void magma_tally2_queue_wait_event( magma_tally2_queue_t queue, magma_tally2_event_t event );


// ========================================
// error handler
void magma_tally2_xerbla( const char *name, magma_tally2_int_t info );

const char* magma_tally2_strerror( magma_tally2_int_t error );


// ========================================
/// For integers x >= 0, y > 0, returns ceil( x/y ).
/// For x == 0, this is 0.
__host__ __device__
static inline magma_tally2_int_t magma_tally2_ceildiv( magma_tally2_int_t x, magma_tally2_int_t y )
{
    return (x + y - 1)/y;
}

/// For integers x >= 0, y > 0, returns x rounded up to multiple of y.
/// For x == 0, this is 0.
/// This implementation does not assume y is a power of 2.
__host__ __device__
static inline magma_tally2_int_t magma_tally2_roundup( magma_tally2_int_t x, magma_tally2_int_t y )
{
    return magma_tally2_ceildiv( x, y ) * y;
}


// ========================================
// real and complex square root
// sqrt alone cannot be caught by the generation script because of tsqrt
static inline float  magma_tally2_ssqrt( float  x ) { return sqrtf( x ); }
static inline double magma_tally2_dsqrt( double x ) { return sqrt( x ); }
magma_tally2FloatComplex    magma_tally2_csqrt( magma_tally2FloatComplex  x );
magma_tally2DoubleComplex   magma_tally2_zsqrt( magma_tally2DoubleComplex x );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally2_H */
