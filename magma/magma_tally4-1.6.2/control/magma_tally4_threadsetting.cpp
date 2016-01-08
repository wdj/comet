/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Mark Gates
*/
#include "common_magma_tally4.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(MAGMA_tally4_WITH_MKL)
#include <mkl_service.h>
#endif

#if defined(HAVE_HWLOC)
#include <hwloc.h>
#endif


/***************************************************************************//**
    Purpose
    -------
    Returns the number of threads to use for parallel sections of MAGMA_tally4.
    Typically, it is initially set by the environment variables
    OMP_NUM_THREADS or MAGMA_tally4_NUM_THREADS.

    If MAGMA_tally4_NUM_THREADS is set, this returns
        min( num_cores, MAGMA_tally4_NUM_THREADS );
    else if MAGMA_tally4 is compiled with OpenMP, this queries OpenMP and returns
        min( num_cores, OMP_NUM_THREADS );
    else this returns num_cores.

    For the number of cores, if MAGMA_tally4 is compiled with hwloc, this queries hwloc;
    else it queries sysconf (on Unix) or GetSystemInfo (on Windows).

    @sa magma_tally4_get_lapack_numthreads
    @sa magma_tally4_set_lapack_numthreads
    @ingroup magma_tally4_util
    ********************************************************************/
extern "C"
magma_tally4_int_t magma_tally4_get_parallel_numthreads()
{
    // query number of cores
    magma_tally4_int_t ncores = 0;

#ifdef HAVE_HWLOC
    // hwloc gives physical cores, not hyperthreads
    // from http://stackoverflow.com/questions/12483399/getting-number-of-cores-not-ht-threads
    hwloc_topology_t topology;
    hwloc_topology_init( &topology );
    hwloc_topology_load( topology );
    magma_tally4_int_t depth = hwloc_get_type_depth( topology, HWLOC_OBJ_CORE );
    if (depth != HWLOC_TYPE_DEPTH_UNKNOWN) {
        ncores = hwloc_get_nbobjs_by_depth( topology, depth );
    }
    hwloc_topology_destroy( topology );
#endif

    if ( ncores == 0 ) {
        #ifdef _MSC_VER  // Windows
        SYSTEM_INFO sysinfo;
        GetSystemInfo( &sysinfo );
        ncores = sysinfo.dwNumberOfProcessors;
        #else
        ncores = sysconf( _SC_NPROCESSORS_ONLN );
        #endif
    }

    // query MAGMA_tally4_NUM_THREADS or OpenMP
    const char *threads_str = getenv("MAGMA_tally4_NUM_THREADS");
    magma_tally4_int_t threads = 0;
    if ( threads_str != NULL ) {
        char* endptr;
        threads = strtol( threads_str, &endptr, 10 );
        if ( threads < 1 || *endptr != '\0' ) {
            threads = 1;
            fprintf( stderr, "$MAGMA_tally4_NUM_THREADS='%s' is an invalid number; using %d thread.\n",
                     threads_str, (int) threads );
        }
    }
    else {
        #if defined(_OPENMP)
        #pragma omp parallel
        {
            threads = omp_get_num_threads();
        }
        #else
            threads = ncores;
        #endif
    }

    // limit to range [1, number of cores]
    threads = max( 1, min( ncores, threads ));
    return threads;
}


/***************************************************************************//**
    Purpose
    -------
    Returns the number of threads currently used for LAPACK and BLAS.
    Typically, the number of threads is initially set by the environment variables
    OMP_NUM_THREADS or MKL_NUM_THREADS.

    If MAGMA_tally4 is compiled with MAGMA_tally4_WITH_MKL, this queries MKL;
    else if MAGMA_tally4 is compiled with OpenMP, this queries OpenMP;
    else this returns 1.
    
    @sa magma_tally4_get_parallel_numthreads
    @sa magma_tally4_set_lapack_numthreads
    @ingroup magma_tally4_util
    ********************************************************************/
extern "C"
magma_tally4_int_t magma_tally4_get_lapack_numthreads()
{
    magma_tally4_int_t threads = 1;

#if defined(MAGMA_tally4_WITH_MKL)
    threads = mkl_get_max_threads();
#elif defined(_OPENMP)
    #pragma omp parallel
    {
        threads = omp_get_num_threads();
    }
#endif

    return threads;
}


/***************************************************************************//**
    Purpose
    -------
    Sets the number of threads to use for LAPACK and BLAS.
    This is often used to set BLAS to be single-threaded during sections
    where MAGMA_tally4 uses explicit pthread parallelism. Example use:

        nthread_save = magma_tally4_get_lapack_numthreads();
        magma_tally4_set_lapack_numthreads( 1 );

        ... launch pthreads, do work, terminate pthreads ...

        magma_tally4_set_lapack_numthreads( nthread_save );

    If MAGMA_tally4 is compiled with MAGMA_tally4_WITH_MKL, this sets MKL threads;
    else if MAGMA_tally4 is compiled with OpenMP, this sets OpenMP threads;
    else this does nothing.

    Arguments
    ---------
    @param[in]
    threads INTEGER
            Number of threads to use. threads >= 1.
            If threads < 1, this silently does nothing.

    @sa magma_tally4_get_parallel_numthreads
    @sa magma_tally4_get_lapack_numthreads
    @ingroup magma_tally4_util
    ********************************************************************/
extern "C"
void magma_tally4_set_lapack_numthreads(magma_tally4_int_t threads)
{
    if ( threads < 1 ) {
        return;
    }

#if defined(MAGMA_tally4_WITH_MKL)
    mkl_set_num_threads( threads );
#elif defined(_OPENMP)
    omp_set_num_threads( threads );
#endif
}





/*
 #else
    #ifndef NO_WARNING
    printf("==============================================================================\n"
           "WARNING: a parallel section of MAGMA_tally4 could not be run in parallel because\n"
           "OpenMP was not enabled; add -fopenmp (for gcc) or -openmp (for icc) to\n"
           "both compilation and linkage flags (CFLAGS and LDFLAGS in make.inc).\n"
           "Also, if using multi-threaded MKL, please compile MAGMA_tally4 with -DMAGMA_tally4_WITH_MKL.\n"
           "To disable this warning, compile MAGMA_tally4 with -DNO_WARNING.\n"
           "==============================================================================\n");
    #endif
   #ifndef NO_WARNING
    printf("==============================================================================\n"
           "WARNING: a parallel section of MAGMA_tally4 could not be run in parallel because\n"
           "OpenMP was not enabled; add -fopenmp (for gcc) or -openmp (for icc) to\n"
           "both compilation and linkage flags (CFLAGS and LDFLAGS in make.inc).\n"
           "To disable this warning, compile MAGMA_tally4 with -DNO_WARNING.\n"
           "==============================================================================\n");
    #endif

*/

