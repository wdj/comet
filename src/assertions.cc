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
