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

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(MAGMA_tally2_WITH_MKL)
#include <mkl_service.h>
#endif

#if defined(MAGMA_tally2_WITH_ACML)
#include <acml.h>
#endif

// defining MAGMA_tally2_LAPACK_H is a hack to NOT include magma_tally2_lapack.h
// via common_magma_tally2.h here, since it conflicts with acml.h and we don't
// need lapack here, but we want acml.h for the acmlversion() function.
#define MAGMA_tally2_LAPACK_H
#include "common_magma_tally2.h"
#include "error.h"

#ifdef HAVE_CUBLAS

// --------------------
// subset of the CUDA device properties, set by magma_tally2_init()
struct magma_tally2_device
{
    size_t memory;
    magma_tally2_int_t cuda_arch;
};

int g_magma_tally2_devices_cnt = 0;
struct magma_tally2_device* g_magma_tally2_devices = NULL;


// ========================================
// initialization
// --------------------
// Caches information about available CUDA devices.
// When renumbering devices after calling magma_tally2_init,
// call magma_tally2_finalize, then cudaSetValidDevices, then magma_tally2_init again.
// Ideally magma_tally2_init is paired with magma_tally2_finalize, but this implementation
// ensures there isn't a memory leak if magma_tally2_init is called multiple times
// without calling magma_tally2_finalize.
extern "C"
magma_tally2_int_t magma_tally2_init()
{
    if ( g_magma_tally2_devices == NULL ) {
        cudaError_t err;
        err = cudaGetDeviceCount( &g_magma_tally2_devices_cnt );
        check_error( err );
        g_magma_tally2_devices = (struct magma_tally2_device*) malloc( g_magma_tally2_devices_cnt * sizeof(struct magma_tally2_device) );
        for( int i = 0; i < g_magma_tally2_devices_cnt; ++i ) {
            cudaDeviceProp prop;
            err = cudaGetDeviceProperties( &prop, i );
            check_error( err );
            g_magma_tally2_devices[i].memory = prop.totalGlobalMem;
            g_magma_tally2_devices[i].cuda_arch  = prop.major*100 + prop.minor*10;
        }
    }
    return MAGMA_tally2_SUCCESS;
}

// --------------------
// Frees information about CUDA devices.
extern "C"
magma_tally2_int_t magma_tally2_finalize()
{
    free( g_magma_tally2_devices );
    g_magma_tally2_devices = NULL;
    return MAGMA_tally2_SUCCESS;
}

// --------------------
// Print the available GPU devices. Used in testing.
extern "C"
void magma_tally2_print_environment()
{
    magma_tally2_int_t major, minor, micro;
    magma_tally2_version( &major, &minor, &micro );
    printf( "MAGMA_tally2 %d.%d.%d %s compiled for CUDA capability >= %.1f\n",
            (int) major, (int) minor, (int) micro, MAGMA_tally2_VERSION_STAGE, MIN_CUDA_ARCH/100. );
    
    int cuda_runtime, cuda_driver;
    cudaError_t err;
    err = cudaDriverGetVersion( &cuda_driver );
    check_error( err );
    err = cudaRuntimeGetVersion( &cuda_runtime );
    check_error( err );
    printf( "CUDA runtime %d, driver %d. ", cuda_runtime, cuda_driver );
    
#if defined(_OPENMP)
    int omp_threads = 0;
    #pragma omp parallel
    {
        omp_threads = omp_get_num_threads();
    }
    printf( "OpenMP threads %d. ", omp_threads );
#else
    printf( "MAGMA_tally2 not compiled with OpenMP. " );
#endif
    
#if defined(MAGMA_tally2_WITH_MKL)
    MKLVersion mkl_version;
    mkl_get_version( &mkl_version );
    printf( "MKL %d.%d.%d, MKL threads %d. ",
            mkl_version.MajorVersion,
            mkl_version.MinorVersion,
            mkl_version.UpdateVersion,
            mkl_get_max_threads() );
#endif
    
#if defined(MAGMA_tally2_WITH_ACML)
    // ACML 4 doesn't have acml_build parameter
    int acml_major, acml_minor, acml_patch, acml_build;
    acmlversion( &acml_major, &acml_minor, &acml_patch, &acml_build );
    printf( "ACML %d.%d.%d.%d ", acml_major, acml_minor, acml_patch, acml_build );
#endif
    
    printf( "\n" );
    
    int ndevices = 0;
    err = cudaGetDeviceCount( &ndevices );
printf( "ndevices %d\n", ndevices );
    check_error( err );
    for( int dev = 0; dev < ndevices; dev++ ) {
        cudaDeviceProp prop;
        err = cudaGetDeviceProperties( &prop, dev );
        check_error( err );
        printf( "device %d: %s, %.1f MHz clock, %.1f MB memory, capability %d.%d\n",
                dev,
                prop.name,
                prop.clockRate / 1000.,
                prop.totalGlobalMem / (1024.*1024.),
                prop.major,
                prop.minor );
        
        int arch = prop.major*100 + prop.minor*10;
        if ( arch < MIN_CUDA_ARCH ) {
            printf("\n"
                   "==============================================================================\n"
                   "WARNING: MAGMA_tally2 was compiled only for CUDA capability %.1f and higher;\n"
                   "device %d has only capability %.1f; some routines will not run correctly!\n"
                   "==============================================================================\n\n",
                   MIN_CUDA_ARCH/100., dev, arch/100. );
        }
    }
}


// ========================================
// device support
// ---------------------------------------------
// Returns CUDA architecture capability for the current device.
// This requires magma_tally2_init to be called first (to cache the information).
// Version is integer xyz, where x is major, y is minor, and z is micro,
// the same as __CUDA_ARCH__. Thus for architecture 1.3 it returns 130.
extern "C"
magma_tally2_int_t magma_tally2_getdevice_arch()
{
    int dev;
    cudaError_t err;
    err = cudaGetDevice( &dev );
    check_error( err );
    if ( g_magma_tally2_devices == NULL || dev < 0 || dev >= g_magma_tally2_devices_cnt ) {
        fprintf( stderr, "Error in %s: MAGMA_tally2 not initialized (call magma_tally2_init() first) or bad device\n", __func__ );
        return 0;
    }
    return g_magma_tally2_devices[dev].cuda_arch;
}


// --------------------
extern "C"
void magma_tally2_getdevices(
    magma_tally2_device_t* devices,
    magma_tally2_int_t     size,
    magma_tally2_int_t*    numPtr )
{
    cudaError_t err;
    int cnt;
    err = cudaGetDeviceCount( &cnt );
    check_error( err );
    
    cnt = min( cnt, size );
    for( int i = 0; i < cnt; ++i ) {
        devices[i] = i;
    }
    *numPtr = cnt;
}

// --------------------
extern "C"
void magma_tally2_getdevice( magma_tally2_device_t* device )
{
    cudaError_t err;
    err = cudaGetDevice( device );
    check_error( err );
}

// --------------------
extern "C"
void magma_tally2_setdevice( magma_tally2_device_t device )
{
    cudaError_t err;
    err = cudaSetDevice( device );
    check_error( err );
}

// --------------------
extern "C"
void magma_tally2_device_sync()
{
    cudaError_t err;
    err = cudaDeviceSynchronize();
    check_error( err );
}


// ========================================
// queue support
// At the moment, MAGMA_tally2 queue == CUDA stream.
// In the future, MAGMA_tally2 queue may be CUBLAS handle.
// --------------------
extern "C"
void magma_tally2_queue_create_internal(
    /*magma_tally2_device_t device,*/ magma_tally2_queue_t* queuePtr,
    const char* func, const char* file, int line )    
{
    //cudaStream_t   stream;
    //cublasStatus_t stat;
    cudaError_t    err;
    //err  = cudaSetDevice( device );
    //stat = cublasCreate( queuePtr );
    err  = cudaStreamCreate( queuePtr );  //&stream );
    //stat = cublasSetStream( *queuePtr, stream );
    check_xerror( err, func, file, line );
}

// --------------------
extern "C"
void magma_tally2_queue_destroy_internal(
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    if ( queue != NULL ) {
        cudaError_t err = cudaStreamDestroy( queue );
        check_xerror( err, func, file, line );
    }
}

// --------------------
extern "C"
void magma_tally2_queue_sync_internal(
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line )
{
    //cudaStream_t   stream;
    //cublasStatus_t stat;
    cudaError_t    err;
    //stat = cublasGetStream( queue, &stream );
    err  = cudaStreamSynchronize( queue );  //stream );
    check_xerror( err, func, file, line );
}


// ========================================
// event support
// --------------------
extern "C"
void magma_tally2_event_create( magma_tally2_event_t* event )
{
    cudaError_t err;
    err = cudaEventCreate( event );
    check_error( err );
}

// --------------------
extern "C"
void magma_tally2_event_destroy( magma_tally2_event_t event )
{
    if ( event != NULL ) {
        cudaError_t err;
        err = cudaEventDestroy( event );
        check_error( err );
    }
}

// --------------------
extern "C"
void magma_tally2_event_record( magma_tally2_event_t event, magma_tally2_queue_t queue )
{
    cudaError_t err;
    err = cudaEventRecord( event, queue );
    check_error( err );
}

// --------------------
// blocks CPU until event occurs
extern "C"
void magma_tally2_event_sync( magma_tally2_event_t event )
{
    cudaError_t err;
    err = cudaEventSynchronize( event );
    check_error( err );
}

// --------------------
// blocks queue (but not CPU) until event occurs
extern "C"
void magma_tally2_queue_wait_event( magma_tally2_queue_t queue, magma_tally2_event_t event )
{
    cudaError_t err;
    err = cudaStreamWaitEvent( queue, event, 0 );
    check_error( err );
}

#endif // HAVE_CUBLAS
