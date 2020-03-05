//-----------------------------------------------------------------------------
/*!
 * \file   assertions.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Macros and code for assertions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdio.h"
#include "stdlib.h"
#include "errno.h"

#ifdef TESTING
#include "gtest/gtest.h"
#endif

#include "env.hh"
#include "assertions.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief For a test build, cause test harness to throw an error.

void trigger_test_harness_failure() {
# ifdef TESTING
    ASSERT_TRUE(0);
# endif
}

#if 0
//-----------------------------------------------------------------------------
/// \brief For a test build, cause test harness to throw an error.

__host__ __device__
static void trigger_test_harness_failure() {
# if defined(TESTING) && ! defined(COMET_DEVICE_COMPILE)
    ASSERT_TRUE(0);
# endif
}

//-----------------------------------------------------------------------------
/// \brief Function to support the COMET_ASSERT, COMET_INSIST macros.

__host__ __device__
void insist(const char* condition_string, const char* file, int line,
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

  trigger_test_harness_failure();

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
# else
    //raise(SIGUSR1);
    exit(EXIT_FAILURE);
#endif
}
#endif

//-----------------------------------------------------------------------------
/// \brief Function to support the COMET_INSIST_INTERFACE macro.

void insist_interface(const void* const env,
                         const char* condition_string,
                         const char* file,
                         int line) {

  const bool do_print = ((CEnv*)env)->proc_num() == 0;

  insist(condition_string, file, line, do_print);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
