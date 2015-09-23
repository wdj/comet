/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_minproduct_MANGLING_H
#define MAGMA_minproduct_MANGLING_H

#include "mangling.h"

/* Define how to name mangle Fortran names.
 * If using CMake, it defines MAGMA_minproduct_GLOBAL in mangling.h
 * Otherwise, the make.inc file should have one of -DADD_, -DNOCHANGE, or -DUPCASE.
 * If using outside of MAGMA_minproduct, put one of those in your compiler flags (e.g., CFLAGS).
 * These macros are used in:
 *   include/magma_minproduct_*lapack.h
 *   control/magma_minproduct_*f77.cpp
 */
#ifndef MAGMA_minproduct_FORTRAN_NAME
    #if defined(MAGMA_minproduct_GLOBAL)
        #define FORTRAN_NAME(lcname, UCNAME)  MAGMA_minproduct_GLOBAL( lcname, UCNAME )
    #elif defined(ADD_)
        #define FORTRAN_NAME(lcname, UCNAME)  lcname##_
    #elif defined(NOCHANGE)
        #define FORTRAN_NAME(lcname, UCNAME)  lcname
    #elif defined(UPCASE)
        #define FORTRAN_NAME(lcname, UCNAME)  UCNAME
    #else
        #error "One of ADD_, NOCHANGE, or UPCASE must be defined to set how Fortran functions are name mangled. For example, in MAGMA_minproduct, add -DADD_ to CFLAGS, FFLAGS, etc. in make.inc. If using CMake, it defines MAGMA_minproduct_GLOBAL instead."
    #endif
#endif

#endif        //  #ifndef MAGMA_minproduct_MANGLING_H
