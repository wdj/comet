/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_tally3BLAS_H
#define MAGMA_tally3BLAS_H

#include "magma_tally3blas_z.h"
#include "magma_tally3blas_c.h"
#include "magma_tally3blas_d.h"
#include "magma_tally3blas_s.h"
#include "magma_tally3blas_zc.h"
#include "magma_tally3blas_ds.h"

#include "magma_tally3blas_z_q.h"
#include "magma_tally3blas_c_q.h"
#include "magma_tally3blas_d_q.h"
#include "magma_tally3blas_s_q.h"
#include "magma_tally3blas_zc_q.h"
#include "magma_tally3blas_ds_q.h"

#ifdef __cplusplus
extern "C" {
#endif

// ========================================
// Define magma_tally3 streams
extern magma_tally3_queue_t magma_tally3_stream;

cublasStatus_t magma_tally3blasSetKernelStream( magma_tally3_queue_t stream );
cublasStatus_t magma_tally3blasGetKernelStream( magma_tally3_queue_t *stream );


// ========================================
// copying vectors
// set copies host to device
// get copies device to host
// Add the function, file, and line for error-reporting purposes.

#define magma_tally3_setvector(           n, elemSize, hx_src, incx, dy_dst, incy ) \
        magma_tally3_setvector_internal(  n, elemSize, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_getvector(           n, elemSize, dx_src, incx, hy_dst, incy ) \
        magma_tally3_getvector_internal(  n, elemSize, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_copyvector(          n, elemSize, dx_src, incx, dy_dst, incy ) \
        magma_tally3_copyvector_internal( n, elemSize, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_setvector_async(           n, elemSize, hx_src, incx, dy_dst, incy, stream ) \
        magma_tally3_setvector_async_internal(  n, elemSize, hx_src, incx, dy_dst, incy, stream, __func__, __FILE__, __LINE__ )

#define magma_tally3_getvector_async(           n, elemSize, dx_src, incx, hy_dst, incy, stream ) \
        magma_tally3_getvector_async_internal(  n, elemSize, dx_src, incx, hy_dst, incy, stream, __func__, __FILE__, __LINE__ )

#define magma_tally3_copyvector_async(          n, elemSize, dx_src, incx, dy_dst, incy, stream ) \
        magma_tally3_copyvector_async_internal( n, elemSize, dx_src, incx, dy_dst, incy, stream, __func__, __FILE__, __LINE__ )

void magma_tally3_setvector_internal(
    magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *hx_src, magma_tally3_int_t incx,
    void       *dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void magma_tally3_getvector_internal(
    magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *dx_src, magma_tally3_int_t incx,
    void       *hy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void magma_tally3_copyvector_internal(
    magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *dx_src, magma_tally3_int_t incx,
    void       *dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void magma_tally3_setvector_async_internal(
    magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *hx_src, magma_tally3_int_t incx,
    void       *dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line );

void magma_tally3_getvector_async_internal(
    magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *dx_src, magma_tally3_int_t incx,
    void       *hy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line );

void magma_tally3_copyvector_async_internal(
    magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *dx_src, magma_tally3_int_t incx,
    void       *dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns )
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_tally3_setmatrix(           m, n, elemSize, hA_src, lda, dB_dst, lddb ) \
        magma_tally3_setmatrix_internal(  m, n, elemSize, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally3_getmatrix(           m, n, elemSize, dA_src, ldda, hB_dst, ldb ) \
        magma_tally3_getmatrix_internal(  m, n, elemSize, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_tally3_copymatrix(          m, n, elemSize, dA_src, ldda, dB_dst, lddb ) \
        magma_tally3_copymatrix_internal( m, n, elemSize, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally3_setmatrix_async(           m, n, elemSize, hA_src, lda, dB_dst, lddb, stream ) \
        magma_tally3_setmatrix_async_internal(  m, n, elemSize, hA_src, lda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

#define magma_tally3_getmatrix_async(           m, n, elemSize, dA_src, ldda, hB_dst, ldb, stream ) \
        magma_tally3_getmatrix_async_internal(  m, n, elemSize, dA_src, ldda, hB_dst, ldb, stream, __func__, __FILE__, __LINE__ )

#define magma_tally3_copymatrix_async(          m, n, elemSize, dA_src, ldda, dB_dst, lddb, stream ) \
        magma_tally3_copymatrix_async_internal( m, n, elemSize, dA_src, ldda, dB_dst, lddb, stream, __func__, __FILE__, __LINE__ )

void magma_tally3_setmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *hA_src, magma_tally3_int_t lda,
    void       *dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line );

void magma_tally3_getmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *dA_src, magma_tally3_int_t ldda,
    void       *hB_dst, magma_tally3_int_t ldb,
    const char* func, const char* file, int line );

void magma_tally3_copymatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *dA_src, magma_tally3_int_t ldda,
    void       *dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line );

void magma_tally3_setmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *hA_src, magma_tally3_int_t lda,
    void       *dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line );

void magma_tally3_getmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *dA_src, magma_tally3_int_t ldda,
    void       *hB_dst, magma_tally3_int_t ldb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line );

void magma_tally3_copymatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t elemSize,
    const void *dA_src, magma_tally3_int_t ldda,
    void       *dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line );


// ========================================
// copying vectors - version for magma_tally3_int_t
// TODO to make these truly type-safe, would need intermediate inline
//      magma_tally3_i* functions that call the generic magma_tally3_* functions.
//      Could do the same with magma_tally3_[sdcz]* set/get functions.

#define magma_tally3_isetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_tally3_isetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_igetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_tally3_igetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )
                                       
#define magma_tally3_icopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_tally3_icopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_isetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_tally3_isetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )
                                             
#define magma_tally3_igetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_tally3_igetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )
                                             
#define magma_tally3_icopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_tally3_icopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

static inline void magma_tally3_isetvector_internal(
    magma_tally3_int_t n,
    const magma_tally3_int_t *hx_src, magma_tally3_int_t incx,
    magma_tally3_int_t       *dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{ magma_tally3_setvector_internal( n, sizeof(magma_tally3_int_t), hx_src, incx, dy_dst, incy, func, file, line ); }

static inline void magma_tally3_igetvector_internal(
    magma_tally3_int_t n,
    const magma_tally3_int_t *dx_src, magma_tally3_int_t incx,
    magma_tally3_int_t       *hy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{ magma_tally3_getvector_internal( n, sizeof(magma_tally3_int_t), dx_src, incx, hy_dst, incy, func, file, line ); }

static inline void magma_tally3_icopyvector_internal(
    magma_tally3_int_t n,
    const magma_tally3_int_t *dx_src, magma_tally3_int_t incx,
    magma_tally3_int_t       *dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{ magma_tally3_copyvector_internal( n, sizeof(magma_tally3_int_t), dx_src, incx, dy_dst, incy, func, file, line ); }

static inline void magma_tally3_isetvector_async_internal(
    magma_tally3_int_t n,
    const magma_tally3_int_t *hx_src, magma_tally3_int_t incx,
    magma_tally3_int_t       *dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_setvector_async_internal( n, sizeof(magma_tally3_int_t), hx_src, incx, dy_dst, incy, stream, func, file, line ); }

static inline void magma_tally3_igetvector_async_internal(
    magma_tally3_int_t n,
    const magma_tally3_int_t *dx_src, magma_tally3_int_t incx,
    magma_tally3_int_t       *hy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_getvector_async_internal( n, sizeof(magma_tally3_int_t), dx_src, incx, hy_dst, incy, stream, func, file, line ); }

static inline void magma_tally3_icopyvector_async_internal(
    magma_tally3_int_t n,
    const magma_tally3_int_t *dx_src, magma_tally3_int_t incx,
    magma_tally3_int_t       *dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_copyvector_async_internal( n, sizeof(magma_tally3_int_t), dx_src, incx, dy_dst, incy, stream, func, file, line ); }


// ========================================
// copying sub-matrices - version for magma_tally3_int_t

#define magma_tally3_isetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_tally3_isetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )
                                          
#define magma_tally3_igetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_tally3_igetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )
                                          
#define magma_tally3_icopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_tally3_icopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally3_isetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_tally3_isetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )
                                                
#define magma_tally3_igetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_tally3_igetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )
                                                
#define magma_tally3_icopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_tally3_icopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

static inline void magma_tally3_isetmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_int_t *hA_src, magma_tally3_int_t lda,
    magma_tally3_int_t       *dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line )
{ magma_tally3_setmatrix_internal( m, n, sizeof(magma_tally3_int_t), hA_src, lda, dB_dst, lddb, func, file, line ); }

static inline void magma_tally3_igetmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_int_t *dA_src, magma_tally3_int_t ldda,
    magma_tally3_int_t       *hB_dst, magma_tally3_int_t ldb,
    const char* func, const char* file, int line )
{ magma_tally3_getmatrix_internal( m, n, sizeof(magma_tally3_int_t), dA_src, ldda, hB_dst, ldb, func, file, line ); }

static inline void magma_tally3_icopymatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_int_t *dA_src, magma_tally3_int_t ldda,
    magma_tally3_int_t       *dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line )
{ magma_tally3_copymatrix_internal( m, n, sizeof(magma_tally3_int_t), dA_src, ldda, dB_dst, lddb, func, file, line ); }

static inline void magma_tally3_isetmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_int_t *hA_src, magma_tally3_int_t lda,
    magma_tally3_int_t       *dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_setmatrix_async_internal( m, n, sizeof(magma_tally3_int_t), hA_src, lda, dB_dst, lddb, stream, func, file, line ); }

static inline void magma_tally3_igetmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_int_t *dA_src, magma_tally3_int_t ldda,
    magma_tally3_int_t       *hB_dst, magma_tally3_int_t ldb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_getmatrix_async_internal( m, n, sizeof(magma_tally3_int_t), dA_src, ldda, hB_dst, ldb, stream, func, file, line ); }

static inline void magma_tally3_icopymatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_int_t *dA_src, magma_tally3_int_t ldda,
    magma_tally3_int_t       *dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_copymatrix_async_internal( m, n, sizeof(magma_tally3_int_t), dA_src, ldda, dB_dst, lddb, stream, func, file, line ); }


// ========================================
// copying vectors - version for magma_tally3_index_t

#define magma_tally3_index_setvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_tally3_index_setvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )
                                            
#define magma_tally3_index_getvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_tally3_index_getvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )
                                            
#define magma_tally3_index_copyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_tally3_index_copyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )
        
#define magma_tally3_index_setvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_tally3_index_setvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )
                                                  
#define magma_tally3_index_getvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_tally3_index_getvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )
                                                  
#define magma_tally3_index_copyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_tally3_index_copyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

static inline void magma_tally3_index_setvector_internal(
    magma_tally3_int_t n,
    const magma_tally3_index_t *hx_src, magma_tally3_int_t incx,
    magma_tally3_index_t       *dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{ magma_tally3_setvector_internal( n, sizeof(magma_tally3_index_t), hx_src, incx, dy_dst, incy, func, file, line ); }

static inline void magma_tally3_index_getvector_internal(
    magma_tally3_int_t n,
    const magma_tally3_index_t *dx_src, magma_tally3_int_t incx,
    magma_tally3_index_t       *hy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{ magma_tally3_getvector_internal( n, sizeof(magma_tally3_index_t), dx_src, incx, hy_dst, incy, func, file, line ); }

static inline void magma_tally3_index_copyvector_internal(
    magma_tally3_int_t n,
    const magma_tally3_index_t *dx_src, magma_tally3_int_t incx,
    magma_tally3_index_t       *dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line )
{ magma_tally3_copyvector_internal( n, sizeof(magma_tally3_index_t), dx_src, incx, dy_dst, incy, func, file, line ); }

static inline void magma_tally3_index_setvector_async_internal(
    magma_tally3_int_t n,
    const magma_tally3_index_t *hx_src, magma_tally3_int_t incx,
    magma_tally3_index_t       *dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_setvector_async_internal( n, sizeof(magma_tally3_index_t), hx_src, incx, dy_dst, incy, stream, func, file, line ); }

static inline void magma_tally3_index_getvector_async_internal(
    magma_tally3_int_t n,
    const magma_tally3_index_t *dx_src, magma_tally3_int_t incx,
    magma_tally3_index_t       *hy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_getvector_async_internal( n, sizeof(magma_tally3_index_t), dx_src, incx, hy_dst, incy, stream, func, file, line ); }

static inline void magma_tally3_index_copyvector_async_internal(
    magma_tally3_int_t n,
    const magma_tally3_index_t *dx_src, magma_tally3_int_t incx,
    magma_tally3_index_t       *dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_copyvector_async_internal( n, sizeof(magma_tally3_index_t), dx_src, incx, dy_dst, incy, stream, func, file, line ); }


// ========================================
// copying sub-matrices - version for magma_tally3_index_t

#define magma_tally3_index_setmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_tally3_index_setmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )
                                               
#define magma_tally3_index_getmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_tally3_index_getmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )
                                               
#define magma_tally3_index_copymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_tally3_index_copymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally3_index_setmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_tally3_index_setmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )
                                                                         
#define magma_tally3_index_getmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_tally3_index_getmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )
                                                     
#define magma_tally3_index_copymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_tally3_index_copymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )


static inline void magma_tally3_index_setmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_index_t *hA_src, magma_tally3_int_t lda,
    magma_tally3_index_t       *dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line )
{ magma_tally3_setmatrix_internal( m, n, sizeof(magma_tally3_index_t), hA_src, lda, dB_dst, lddb, func, file, line ); }

static inline void magma_tally3_index_getmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_index_t *dA_src, magma_tally3_int_t ldda,
    magma_tally3_index_t       *hB_dst, magma_tally3_int_t ldb,
    const char* func, const char* file, int line )
{ magma_tally3_getmatrix_internal( m, n, sizeof(magma_tally3_index_t), dA_src, ldda, hB_dst, ldb, func, file, line ); }

static inline void magma_tally3_index_copymatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_index_t *dA_src, magma_tally3_int_t ldda,
    magma_tally3_index_t       *dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line )
{ magma_tally3_copymatrix_internal( m, n, sizeof(magma_tally3_index_t), dA_src, ldda, dB_dst, lddb, func, file, line ); }

static inline void magma_tally3_index_setmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_index_t *hA_src, magma_tally3_int_t lda,
    magma_tally3_index_t       *dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_setmatrix_async_internal( m, n, sizeof(magma_tally3_index_t), hA_src, lda, dB_dst, lddb, stream, func, file, line ); }

static inline void magma_tally3_index_getmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_index_t *dA_src, magma_tally3_int_t ldda,
    magma_tally3_index_t       *hB_dst, magma_tally3_int_t ldb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_getmatrix_async_internal( m, n, sizeof(magma_tally3_index_t), dA_src, ldda, hB_dst, ldb, stream, func, file, line ); }

static inline void magma_tally3_index_copymatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3_index_t *dA_src, magma_tally3_int_t ldda,
    magma_tally3_index_t       *dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t stream,
    const char* func, const char* file, int line )
{ magma_tally3_copymatrix_async_internal( m, n, sizeof(magma_tally3_index_t), dA_src, ldda, dB_dst, lddb, stream, func, file, line ); }


#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally3BLAS_H */
