//-----------------------------------------------------------------------------
/*!
 * \file   assertions.i.hh
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

#ifndef _COMET_ASSERTIONS_I_HH_
#define _COMET_ASSERTIONS_I_HH_

#include "stdio.h"
#include "stdlib.h"
#include "errno.h"
#include "unistd.h"

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
      const size_t hnlen = 256;
      char hn[hnlen];
      gethostname(hn, hnlen);
      fprintf(stderr, "%s: \"%s\". At file %s, line %i, host \"%s\".\n", "Assertion error",
              condition_string, file, line, hn);
#   endif
  }
  //abort();

# ifndef COMET_DEVICE_COMPILE
    trigger_test_harness_failure();
# endif

# ifdef __CUDA_ARCH__
#   ifdef NDEBUG
      // See https://stackoverflow.com/questions/50755717/triggering-a-runtime-error-within-a-cuda-kernel
      asm("trap;");
#   else
      // CUDA kernel assertion.
      // NOTE: standard assert in cuda is no-op if NDEBUG defined
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

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_ASSERTIONS_I_HH_

//-----------------------------------------------------------------------------
