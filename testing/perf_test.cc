//-----------------------------------------------------------------------------
/*!
 * \file   driver_test.cc
 * \author Wayne Joubert
 * \date   Fri Nov  6 18:18:21 EST 2015
 * \brief  Tester for driver.
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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "string"

#ifdef COMET_USE_GTEST
# include "gtest/gtest.h"
#else
# define GTEST_API_
# define EXPECT_EQ(a, b) COMET_INSIST((a) == (b));
#endif

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "checksum.hh"
#include "compute_metrics.hh"

#include "driver.hh"
#include "vectors_io.hh"
#include "metrics_io.hh"
#include "test_problems.hh"

//=============================================================================

void DriverTest_perf_() {
#ifdef COMET_PLATFORM_IBM_AC922

  using namespace comet;

  const bool is_proc_num_0 = System::is_proc_num_0();

  const char* options[2] = {
    "--num_field_local 393216 --num_vector_local 9984 --metric_type ccc "
    "--sparse yes --all2all yes --compute_method GPU --problem_type random "
    "--checksum yes --num_proc_vector 12 "
    "--verbosity 1 --tc 1 --num_tc_steps 4",
    //
    "--num_way 3 --num_field_local 400000 --num_vector_local 5952 "
    "--metric_type ccc --sparse yes --all2all yes --compute_method GPU "
    "--problem_type random --checksum yes --num_proc_vector 2 "
    "--num_proc_repl 2 --num_stage 248 --stage_min 247 --verbosity 1 "
    "--tc 1 --num_tc_steps 4 "};

  // NOTE: running debug version of code, slightly slower than release version.
  const double perf_target[2] = {70e12, 70e+12};

  bool is_passed = true;

  for (int i=0; i<2; ++i) {

    if (is_proc_num_0)
      printf("%s\n", options[i]);

    CEnv env(MPI_COMM_WORLD, options[i]);

    perform_run(options[i], MPI_COMM_WORLD, &env);

    const double ops = env.ops();
    const double ops_rate_proc = ops / (env.ctime() * env.num_proc());

     is_passed = is_passed && is_proc_num_0 ?
       ops_rate_proc >= perf_target[i] : true;;

  }

  EXPECT_EQ(is_passed, true);

#else
  EXPECT_EQ(true, true);
#endif
}

//=============================================================================

#ifdef COMET_USE_GTEST

TEST(DriverTest, perf) {
  DriverTest_perf_();
}

#else

int RUN_ALL_TESTS() {
  DriverTest_perf_();
  return 0;
}

#endif

//=============================================================================

GTEST_API_ int main(int argc, char** argv) {

# ifdef COMET_USE_GTEST
    ::testing::InitGoogleTest(&argc, argv);
# endif

  COMET_MPI_SAFE_CALL(MPI_Init(&argc, &argv));

  int comm_rank = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank));

  if (comm_rank != 0) {
#   ifdef COMET_USE_GTEST
      ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
      delete listeners.Release(listeners.default_result_printer());
#   endif
  }

  int result = RUN_ALL_TESTS();

  int result_g = 11;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&result, &result_g, 1, MPI_INT, MPI_MAX,
    MPI_COMM_WORLD));

  COMET_MPI_SAFE_CALL(MPI_Finalize());
  return result_g;
}

//=============================================================================

//-----------------------------------------------------------------------------
