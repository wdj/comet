//-----------------------------------------------------------------------------
/*!
 * \file   assertions.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Macros and code for assertions.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#ifndef _COMET_ASSERTIONS_HH_
#define _COMET_ASSERTIONS_HH_

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

#if ! defined(NDEBUG)
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
// Implementation include files.

#include "assertions.i.hh"

#endif // _COMET_ASSERTIONS_HH_

//-----------------------------------------------------------------------------
