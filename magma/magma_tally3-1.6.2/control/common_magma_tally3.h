/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mathieu Faverge
 
       Based on PLASMA common.h
*/

/***************************************************************************//**
 *  MAGMA_tally3 facilities of interest to both src and magma_tally3blas directories
 **/
#ifndef MAGMA_tally3_COMMON_H
#define MAGMA_tally3_COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#if defined( _WIN32 ) || defined( _WIN64 )

    #include "magma_tally3winthread.h"
    #include <windows.h>
    #include <limits.h>
    #include <io.h>

    // functions where Microsoft fails to provide C99 standard
    // (only with Microsoft, not with nvcc on Windows)
    // in both common_magma_tally3.h and testings.h
    #ifndef __NVCC__
    
        #include <float.h>
        #define copysign(x,y) _copysign(x,y)
        #define isnan(x)      _isnan(x)
        #define isinf(x)      ( ! _finite(x) && ! _isnan(x) )
        #define isfinite(x)   _finite(x)
        // note _snprintf has slightly different semantics than snprintf
        #define snprintf _snprintf
    
    #endif

#else

    #include <pthread.h>
    #include <unistd.h>
    #include <inttypes.h>

#endif

// provide our own support for pthread_barrier on MacOS and Windows
#include "pthread_barrier.h"

#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "magma_tally3_operators.h"
#include "transpose.h"
#include "magma_tally3_threadsetting.h"

/** ****************************************************************************
 *  Determine if weak symbols are allowed
 */
#if defined(linux) || defined(__linux) || defined(__linux__)
#if defined(__GNUC_EXCL__) || defined(__GNUC__)
#define MAGMA_tally3_HAVE_WEAK    1
#endif
#endif

/***************************************************************************//**
 *  Global utilities
 *  in both common_magma_tally3.h and testings.h
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
#define MAGMA_tally3_UNUSED(var)  ((void)var)


/** ****************************************************************************
 *  Define magma_tally3_[sd]sqrt functions
 *    - sqrt alone cannot be caught by the generation script because of tsqrt
 */

#define magma_tally3_dsqrt sqrt
#define magma_tally3_ssqrt sqrtf

#endif /* MAGMA_tally3_COMMON_H */
