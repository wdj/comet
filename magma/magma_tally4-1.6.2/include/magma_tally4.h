/*
    -- MAGMA_tally4 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date
*/

#ifndef MAGMA_tally4_H
#define MAGMA_tally4_H

/* ------------------------------------------------------------
 * MAGMA_tally4 BLAS Functions
 * --------------------------------------------------------- */
#include "magma_tally4blas.h"
#include "magma_tally4_batched.h"

/* ------------------------------------------------------------
 * MAGMA_tally4 functions
 * --------------------------------------------------------- */
#include "magma_tally4_z.h"
#include "magma_tally4_c.h"
#include "magma_tally4_d.h"
#include "magma_tally4_s.h"
#include "magma_tally4_zc.h"
#include "magma_tally4_ds.h"
#include "auxiliary.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 Auxiliary functions to get the NB used
*/
magma_tally4_int_t magma_tally4_get_smlsize_divideconquer();

// ========================================
// initialization
magma_tally4_int_t
magma_tally4_init( void );

magma_tally4_int_t
magma_tally4_finalize( void );

void magma_tally4_version( magma_tally4_int_t* major, magma_tally4_int_t* minor, magma_tally4_int_t* micro );


// ========================================
// memory allocation
magma_tally4_int_t
magma_tally4_malloc( magma_tally4_ptr *ptrPtr, size_t bytes );

magma_tally4_int_t
magma_tally4_malloc_cpu( void **ptrPtr, size_t bytes );

magma_tally4_int_t
magma_tally4_malloc_pinned( void **ptrPtr, size_t bytes );

magma_tally4_int_t
magma_tally4_free_cpu( void *ptr );

#define magma_tally4_free( ptr ) \
        magma_tally4_free_internal( ptr, __func__, __FILE__, __LINE__ )

#define magma_tally4_free_pinned( ptr ) \
        magma_tally4_free_pinned_internal( ptr, __func__, __FILE__, __LINE__ )

magma_tally4_int_t
magma_tally4_free_internal(
    magma_tally4_ptr ptr,
    const char* func, const char* file, int line );

magma_tally4_int_t
magma_tally4_free_pinned_internal(
    void *ptr,
    const char* func, const char* file, int line );


// type-safe convenience functions to avoid using (void**) cast and sizeof(...)
// here n is the number of elements (floats, doubles, etc.) not the number of bytes.
static inline magma_tally4_int_t magma_tally4_imalloc( magma_tally4Int_ptr           *ptrPtr, size_t n ) { return magma_tally4_malloc( (magma_tally4_ptr*) ptrPtr, n*sizeof(magma_tally4_int_t)        ); }
static inline magma_tally4_int_t magma_tally4_index_malloc( magma_tally4Index_ptr    *ptrPtr, size_t n ) { return magma_tally4_malloc( (magma_tally4_ptr*) ptrPtr, n*sizeof(magma_tally4_index_t)      ); }
static inline magma_tally4_int_t magma_tally4_smalloc( magma_tally4Float_ptr         *ptrPtr, size_t n ) { return magma_tally4_malloc( (magma_tally4_ptr*) ptrPtr, n*sizeof(float)              ); }
static inline magma_tally4_int_t magma_tally4_dmalloc( magma_tally4Double_ptr        *ptrPtr, size_t n ) { return magma_tally4_malloc( (magma_tally4_ptr*) ptrPtr, n*sizeof(double)             ); }
static inline magma_tally4_int_t magma_tally4_cmalloc( magma_tally4FloatComplex_ptr  *ptrPtr, size_t n ) { return magma_tally4_malloc( (magma_tally4_ptr*) ptrPtr, n*sizeof(magma_tally4FloatComplex)  ); }
static inline magma_tally4_int_t magma_tally4_zmalloc( magma_tally4DoubleComplex_ptr *ptrPtr, size_t n ) { return magma_tally4_malloc( (magma_tally4_ptr*) ptrPtr, n*sizeof(magma_tally4DoubleComplex) ); }

static inline magma_tally4_int_t magma_tally4_imalloc_cpu( magma_tally4_int_t        **ptrPtr, size_t n ) { return magma_tally4_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_tally4_int_t)        ); }
static inline magma_tally4_int_t magma_tally4_index_malloc_cpu( magma_tally4_index_t **ptrPtr, size_t n ) { return magma_tally4_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_tally4_index_t)      ); }
static inline magma_tally4_int_t magma_tally4_smalloc_cpu( float              **ptrPtr, size_t n ) { return magma_tally4_malloc_cpu( (void**) ptrPtr, n*sizeof(float)              ); }
static inline magma_tally4_int_t magma_tally4_dmalloc_cpu( double             **ptrPtr, size_t n ) { return magma_tally4_malloc_cpu( (void**) ptrPtr, n*sizeof(double)             ); }
static inline magma_tally4_int_t magma_tally4_cmalloc_cpu( magma_tally4FloatComplex  **ptrPtr, size_t n ) { return magma_tally4_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_tally4FloatComplex)  ); }
static inline magma_tally4_int_t magma_tally4_zmalloc_cpu( magma_tally4DoubleComplex **ptrPtr, size_t n ) { return magma_tally4_malloc_cpu( (void**) ptrPtr, n*sizeof(magma_tally4DoubleComplex) ); }

static inline magma_tally4_int_t magma_tally4_imalloc_pinned( magma_tally4_int_t        **ptrPtr, size_t n ) { return magma_tally4_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_tally4_int_t)        ); }
static inline magma_tally4_int_t magma_tally4_index_malloc_pinned( magma_tally4_index_t **ptrPtr, size_t n ) { return magma_tally4_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_tally4_index_t)      ); }
static inline magma_tally4_int_t magma_tally4_smalloc_pinned( float              **ptrPtr, size_t n ) { return magma_tally4_malloc_pinned( (void**) ptrPtr, n*sizeof(float)              ); }
static inline magma_tally4_int_t magma_tally4_dmalloc_pinned( double             **ptrPtr, size_t n ) { return magma_tally4_malloc_pinned( (void**) ptrPtr, n*sizeof(double)             ); }
static inline magma_tally4_int_t magma_tally4_cmalloc_pinned( magma_tally4FloatComplex  **ptrPtr, size_t n ) { return magma_tally4_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_tally4FloatComplex)  ); }
static inline magma_tally4_int_t magma_tally4_zmalloc_pinned( magma_tally4DoubleComplex **ptrPtr, size_t n ) { return magma_tally4_malloc_pinned( (void**) ptrPtr, n*sizeof(magma_tally4DoubleComplex) ); }


// ========================================
// device support
magma_tally4_int_t magma_tally4_getdevice_arch();

void magma_tally4_getdevices(
    magma_tally4_device_t* devices,
    magma_tally4_int_t     size,
    magma_tally4_int_t*    numPtr );

void magma_tally4_getdevice( magma_tally4_device_t* dev );

void magma_tally4_setdevice( magma_tally4_device_t dev );

void magma_tally4_device_sync();


// ========================================
// queue support
#define magma_tally4_queue_create( /*device,*/ queuePtr ) \
        magma_tally4_queue_create_internal( queuePtr, __func__, __FILE__, __LINE__ )

#define magma_tally4_queue_destroy( queue ) \
        magma_tally4_queue_destroy_internal( queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_queue_sync( queue ) \
        magma_tally4_queue_sync_internal( queue, __func__, __FILE__, __LINE__ )

void magma_tally4_queue_create_internal(
    /*magma_tally4_device_t device,*/ magma_tally4_queue_t* queuePtr,
    const char* func, const char* file, int line );

void magma_tally4_queue_destroy_internal(
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void magma_tally4_queue_sync_internal(
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

/// Currently, magma_tally4_queue_t == cudaStream_t.
/// Almost certainly this will change in the future,
/// so these get & set the associated stream in a forward-compatible manner.
/// @see magma_tally4_queue_set_cuda_stream
#define magma_tally4_queue_get_cuda_stream( queue ) (queue)

/// @see magma_tally4_queue_get_cuda_stream
#define magma_tally4_queue_set_cuda_stream( queue, stream ) ((queue) = (stream))


// ========================================
// event support
void magma_tally4_event_create( magma_tally4_event_t* eventPtr );

void magma_tally4_event_destroy( magma_tally4_event_t event );

void magma_tally4_event_record( magma_tally4_event_t event, magma_tally4_queue_t queue );

void magma_tally4_event_query( magma_tally4_event_t event );

// blocks CPU until event occurs
void magma_tally4_event_sync( magma_tally4_event_t event );

// blocks queue (but not CPU) until event occurs
void magma_tally4_queue_wait_event( magma_tally4_queue_t queue, magma_tally4_event_t event );


// ========================================
// error handler
void magma_tally4_xerbla( const char *name, magma_tally4_int_t info );

const char* magma_tally4_strerror( magma_tally4_int_t error );


// ========================================
/// For integers x >= 0, y > 0, returns ceil( x/y ).
/// For x == 0, this is 0.
__host__ __device__
static inline magma_tally4_int_t magma_tally4_ceildiv( magma_tally4_int_t x, magma_tally4_int_t y )
{
    return (x + y - 1)/y;
}

/// For integers x >= 0, y > 0, returns x rounded up to multiple of y.
/// For x == 0, this is 0.
/// This implementation does not assume y is a power of 2.
__host__ __device__
static inline magma_tally4_int_t magma_tally4_roundup( magma_tally4_int_t x, magma_tally4_int_t y )
{
    return magma_tally4_ceildiv( x, y ) * y;
}


// ========================================
// real and complex square root
// sqrt alone cannot be caught by the generation script because of tsqrt
static inline float  magma_tally4_ssqrt( float  x ) { return sqrtf( x ); }
static inline double magma_tally4_dsqrt( double x ) { return sqrt( x ); }
magma_tally4FloatComplex    magma_tally4_csqrt( magma_tally4FloatComplex  x );
magma_tally4DoubleComplex   magma_tally4_zsqrt( magma_tally4DoubleComplex x );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally4_H */
