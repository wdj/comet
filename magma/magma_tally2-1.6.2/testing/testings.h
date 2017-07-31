#ifndef TESTINGS_H
#define TESTINGS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef HAVE_CUBLAS
#include <cublas_v2.h>  // before magma_tally2.h
#endif

#include "magma_tally2.h"


/***************************************************************************//**
 *  For portability to Windows
 */
#if defined( _WIN32 ) || defined( _WIN64 )
    // functions where Microsoft fails to provide C99 or POSIX standard
    // (only with Microsoft, not with nvcc on Windows)
    // in both common_magma_tally2.h and testings.h
    #ifndef __NVCC__
    
        #include <float.h>
        #define copysign(x,y) _copysign(x,y)
        #define isnan(x)      _isnan(x)
        #define isinf(x)      ( ! _finite(x) && ! _isnan(x) )
        #define isfinite(x)   _finite(x)
        // note _snprintf has slightly different semantics than snprintf
        #define snprintf      _snprintf
        #define unlink        _unlink
        
    #endif
#endif


#ifdef __cplusplus
extern "C" {
#endif

void flops_init();

/***************************************************************************//**
 *  Global utilities
 *  in both common_magma_tally2.h and testings.h
 **/
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

// for integers a >  0, b > 0, returns ceil( a/b ).
// for integers a == 0, b > 0, returns 1.
#ifndef ceildiv
#define ceildiv(a, b) ((a - 1)/b + 1)
#endif

// for integers a >  0, b > 0, returns a rounded up to multiple of b.
// for integers a == 0, b > 0, returns b.
// old implementation assumes b is power of 2:
// (b <= 0) ? (a) : (((a) + (b)-1) & ~((b)-1))
#ifndef roundup
#define roundup(a, b) (ceildiv((a), (b)) * (b))
#endif

// suppress "warning: unused variable" in a portable fashion
#define MAGMA_tally2_UNUSED(var)  ((void)var)


/***************************************************************************//**
 * Macros to handle error checking.
 */

#define TESTING_INIT()                                                     \
    magma_tally2_init();                                                          \
    flops_init();                                                          \
    magma_tally2_print_environment();

#define TESTING_FINALIZE()                                                 \
    magma_tally2_finalize();


/******************* CPU memory */
#define TESTING_MALLOC_CPU( ptr, type, size )                              \
    if ( MAGMA_tally2_SUCCESS !=                                                  \
            magma_tally2_malloc_cpu( (void**) &ptr, (size)*sizeof(type) )) {      \
        fprintf( stderr, "!!!! magma_tally2_malloc_cpu failed for: %s\n", #ptr ); \
        magma_tally2_finalize();                                                  \
        exit(-1);                                                          \
    }

#define TESTING_FREE_CPU( ptr ) magma_tally2_free_cpu( ptr )


/******************* Pinned CPU memory */
#ifdef HAVE_CUBLAS
    // In CUDA, this allocates pinned memory.
    #define TESTING_MALLOC_PIN( ptr, type, size )                                 \
        if ( MAGMA_tally2_SUCCESS !=                                                     \
                magma_tally2_malloc_pinned( (void**) &ptr, (size)*sizeof(type) )) {      \
            fprintf( stderr, "!!!! magma_tally2_malloc_pinned failed for: %s\n", #ptr ); \
            magma_tally2_finalize();                                                     \
            exit(-1);                                                             \
        }
    
    #define TESTING_FREE_PIN( ptr ) magma_tally2_free_pinned( ptr )
#else
    // For OpenCL, we don't support pinned memory yet.
    #define TESTING_MALLOC_PIN( ptr, type, size )                              \
        if ( MAGMA_tally2_SUCCESS !=                                                  \
                magma_tally2_malloc_cpu( (void**) &ptr, (size)*sizeof(type) )) {      \
            fprintf( stderr, "!!!! magma_tally2_malloc_cpu failed for: %s\n", #ptr ); \
            magma_tally2_finalize();                                                  \
            exit(-1);                                                          \
        }
    
    #define TESTING_FREE_PIN( ptr ) magma_tally2_free_cpu( ptr )
#endif


/******************* GPU memory */
#ifdef HAVE_CUBLAS
    // In CUDA, this has (void**) cast.
    #define TESTING_MALLOC_DEV( ptr, type, size )                              \
        if ( MAGMA_tally2_SUCCESS !=                                                  \
                magma_tally2_malloc( (void**) &ptr, (size)*sizeof(type) )) {          \
            fprintf( stderr, "!!!! magma_tally2_malloc failed for: %s\n", #ptr );     \
            magma_tally2_finalize();                                                  \
            exit(-1);                                                          \
        }
#else
    // For OpenCL, ptr is cl_mem* and there is no cast.
    #define TESTING_MALLOC_DEV( ptr, type, size )                              \
        if ( MAGMA_tally2_SUCCESS !=                                                  \
                magma_tally2_malloc( &ptr, (size)*sizeof(type) )) {                   \
            fprintf( stderr, "!!!! magma_tally2_malloc failed for: %s\n", #ptr );     \
            magma_tally2_finalize();                                                  \
            exit(-1);                                                          \
        }
#endif

#define TESTING_FREE_DEV( ptr ) magma_tally2_free( ptr )


/***************************************************************************//**
 * Functions and data structures used for testing.
 */
void magma_tally2_zmake_hermitian( magma_tally2_int_t N, magma_tally2DoubleComplex* A, magma_tally2_int_t lda );
void magma_tally2_cmake_hermitian( magma_tally2_int_t N, magma_tally2FloatComplex*  A, magma_tally2_int_t lda );
void magma_tally2_dmake_symmetric( magma_tally2_int_t N, double*             A, magma_tally2_int_t lda );
void magma_tally2_smake_symmetric( magma_tally2_int_t N, float*              A, magma_tally2_int_t lda );

void magma_tally2_zmake_hpd( magma_tally2_int_t N, magma_tally2DoubleComplex* A, magma_tally2_int_t lda );
void magma_tally2_cmake_hpd( magma_tally2_int_t N, magma_tally2FloatComplex*  A, magma_tally2_int_t lda );
void magma_tally2_dmake_hpd( magma_tally2_int_t N, double*             A, magma_tally2_int_t lda );
void magma_tally2_smake_hpd( magma_tally2_int_t N, float*              A, magma_tally2_int_t lda );

void magma_tally2_assert( bool condition, const char* msg, ... );

void magma_tally2_assert_warn( bool condition, const char* msg, ... );

#define MAX_NTEST 1050

typedef struct magma_tally2_opts
{
    // matrix size
    magma_tally2_int_t ntest;
    magma_tally2_int_t msize[ MAX_NTEST ];
    magma_tally2_int_t nsize[ MAX_NTEST ];
    magma_tally2_int_t ksize[ MAX_NTEST ];
    magma_tally2_int_t mmax;
    magma_tally2_int_t nmax;
    magma_tally2_int_t kmax;
    magma_tally2_int_t batchcount;
    
    // scalars
    magma_tally2_int_t device;
    magma_tally2_int_t roundup;
    magma_tally2_int_t nb;
    magma_tally2_int_t nrhs;
    magma_tally2_int_t nstream;
    magma_tally2_int_t ngpu;
    magma_tally2_int_t niter;
    magma_tally2_int_t nthread;
    magma_tally2_int_t offset;
    magma_tally2_int_t itype;     // hegvd: problem type
    magma_tally2_int_t svd_work;  // gesvd
    magma_tally2_int_t version;   // hemm_mgpu, hetrd
    double      fraction;  // hegvdx
    double      tolerance;
    magma_tally2_int_t panel_nthread; //in magma_tally2_amc: first dimension for a 2D big panel
    double fraction_dcpu; //in magma_tally2_amc: fraction of the work for the cpu 
    // boolean arguments
    int check;
    int lapack;
    int warmup;
    int all;
    int verbose;
    
    // lapack flags
    magma_tally2_uplo_t    uplo;
    magma_tally2_trans_t   transA;
    magma_tally2_trans_t   transB;
    magma_tally2_side_t    side;
    magma_tally2_diag_t    diag;
    magma_tally2_vec_t     jobu;    // gesvd:  no left  singular vectors
    magma_tally2_vec_t     jobvt;   // gesvd:  no right singular vectors
    magma_tally2_vec_t     jobz;    // heev:   no eigen vectors
    magma_tally2_vec_t     jobvr;   // geev:   no right eigen vectors
    magma_tally2_vec_t     jobvl;   // geev:   no left  eigen vectors
    
    // queue for default device
    magma_tally2_queue_t   queue;
    magma_tally2_queue_t   queues2[3];  // 2 queues + 1 extra NULL entry to catch errors
    
    #ifdef HAVE_CUBLAS
    // handle for directly calling cublas
    cublasHandle_t  handle;
    #endif
    
    // misc
    int flock_op;   // shared or exclusive lock
    int flock_fd;   // lock file
} magma_tally2_opts;

void parse_opts( int argc, char** argv, magma_tally2_opts *opts );

extern const char* g_platform_str;

#ifdef __cplusplus
}
#endif

#endif /* TESTINGS_H */
