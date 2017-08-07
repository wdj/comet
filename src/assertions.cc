/*---------------------------------------------------------------------------*/
/*!
 * \file   assertions.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Assertions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <errno.h>

#ifdef TESTING
#include "gtest/gtest.h"
#endif

#include "assertions.hh"
#include "env.hh"

/*---------------------------------------------------------------------------*/

static void gm_test_wrapper() {
#ifdef TESTING
  ASSERT_TRUE(0);
#endif
}

/*===========================================================================*/
/*---Assertions---*/

void gm_assert(const char* condition_string, const char* file, int line) {
  fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Assertion error",
          condition_string, file, line);
  gm_test_wrapper();
#ifdef GM_ASSERTIONS_ON
  raise(SIGUSR1);
#else
  exit(EXIT_FAILURE);
#endif
}

/*---------------------------------------------------------------------------*/

void gm_insist(const void* const env,
               const char* condition_string,
               const char* file,
               int line) {
  if (GMEnv_proc_num((const GMEnv* const)env) == 0) {
    fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Interface error",
            condition_string, file, line);
  }
  gm_test_wrapper();
  exit(EXIT_FAILURE);
}

/*---------------------------------------------------------------------------*/
