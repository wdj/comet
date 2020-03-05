//-----------------------------------------------------------------------------
/*!
 * \file   assertions.i.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Macros and code for assertions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_assertions_i_hh_
#define _comet_assertions_i_hh_

#include "stdio.h"
#include "stdlib.h"
#include "errno.h"

#include "assertions.hh"

//# if defined(TESTING) && ! defined(COMET_DEVICE_COMPILE)
//FIX#include "gtest/gtest.h"
//#endif

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief For a test build, cause test harness to throw an error.

//__host__ __device__
//static void trigger_test_harness_failure() {
//# if defined(TESTING) && ! defined(COMET_DEVICE_COMPILE)
////FIX    ASSERT_TRUE(0);
//# endif
//}

//-----------------------------------------------------------------------------
/// \brief Function to support the COMET_ASSERT, COMET_INSIST macros.

__host__ __device__
static void insist(const char* condition_string, const char* file, int line,
                   bool do_print) {

  if (do_print) {
#   ifdef COMET_DEVICE_COMPILE
      printf("%s: \"%s\". At file %s, line %i.\n", "Assertion error",
             condition_string, file, line);
#   else
      fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Assertion error",
              condition_string, file, line);
#   endif
  }

# ifndef COMET_DEVICE_COMPILE
    trigger_test_harness_failure();
# endif

# ifdef __CUDA_ARCH__
#   ifdef NDEBUG
      // See https://stackoverflow.com/questions/50755717/triggering-a-runtime-error-within-a-cuda-kernel
      asm("trap;");
#   else
      // CUDA kernel assertion.
      assert(false);
#   endif
# elif __HIP_DEVICE_COMPILE__
    // Crude method until kernel assertions implemented later.
    abort();
# else // ! defined(COMET_DEVICE_COMPILE)
    //raise(SIGUSR1);
    exit(EXIT_FAILURE);
#endif
}

//=============================================================================

#endif // _comet_assertions_i_hh_

} // namespace comet

//-----------------------------------------------------------------------------
