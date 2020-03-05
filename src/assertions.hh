//-----------------------------------------------------------------------------
/*!
 * \file   assertions.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Macros and code for assertions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_assertions_hh_
#define _comet_assertions_hh_

#include <type_traits>

#if defined COMET_USE_CUDA
# ifdef __CUDA_ARCH__
#   define COMET_DEVICE_COMPILE
# endif
#elif defined COMET_USE_HIP
# ifdef __HIP_DEVICE_COMPILE__
#   define COMET_DEVICE_COMPILE
# endif
#endif

//=============================================================================
// Macros.

//-----------------------------------------------------------------------------
/// \brief Insist macro - assert a condition even for release builds.
///
///        This should be used only for non-performance-sensitive
///        code locations -- e.g., not in a deep loop nest.

#define COMET_INSIST(condition) \
    (void)((condition) || \
      (comet::insist(#condition, __FILE__, __LINE__, true), 0))

//-----------------------------------------------------------------------------
/// \brief Insist macro specifically for a user-caused error condition.

# define COMET_INSIST_INTERFACE(env, condition) \
    (void)((condition) || \
      (comet::insist_interface(env, #condition, __FILE__, __LINE__), 0))

//-----------------------------------------------------------------------------
/// \brief Assertion macro (for debug builds only).

#ifndef NDEBUG
# define COMET_ASSERTIONS_ON
# define COMET_ASSERT(condition) COMET_INSIST(condition)
#else
# define COMET_ASSERT(condition)
#endif

//-----------------------------------------------------------------------------
/// \brief Static assert - use macro to ensure removal for release build.

#ifndef NDEBUG
# define COMET_STATIC_ASSERT(condition) static_assert(condition, "")
#else
# define COMET_STATIC_ASSERT(condition)
#endif

//=============================================================================
// Declarations.

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Function to support the COMET_ASSERT, COMET_INSIST macros.

//        Use trailing underscore to avoid collision with C assert macro.

__host__ __device__
static void insist(const char* condition_string, const char* file, int line,
                   bool do_print);

//-----------------------------------------------------------------------------
/// \brief Function to support the COMET_INSIST_INTERFACE macro.

void insist_interface(void const * const env,
                      const char* condition_string,
                      const char* file,
                      int line);

//-----------------------------------------------------------------------------
/// \brief For a test build, cause test harness to throw an error.

void trigger_test_harness_failure();

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#include "assertions.i.hh"

//-----------------------------------------------------------------------------

#endif // _comet_assertions_hh_

//-----------------------------------------------------------------------------
