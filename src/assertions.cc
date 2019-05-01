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

namespace CoMet {

//-----------------------------------------------------------------------------
/// \brief For a test build, cause test harness to throw an error.

static void make_test_harness_failure() {
#ifdef TESTING
  ASSERT_TRUE(0);
#endif
}

//-----------------------------------------------------------------------------
/// \brief Function to support the GMAssert, GMInsist macros.

void assert_(const char* condition_string, const char* file, int line) {
  fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Assertion error",
          condition_string, file, line);
  make_test_harness_failure();
#ifdef GM_ASSERTIONS_ON
  //raise(SIGUSR1);
  exit(EXIT_FAILURE);
#else
  exit(EXIT_FAILURE);
#endif
}

//-----------------------------------------------------------------------------
/// \brief Function to support the GMInsistInterface macro.

void insist_interface(const void* const env,
                         const char* condition_string,
                         const char* file,
                         int line) {
  if (GMEnv_proc_num((const GMEnv* const)env) == 0) {
    fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Interface error",
            condition_string, file, line);
  }
  make_test_harness_failure();
  exit(EXIT_FAILURE);
}

//=============================================================================

} // namespace CoMet

//-----------------------------------------------------------------------------
