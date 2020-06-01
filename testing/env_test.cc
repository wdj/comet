//-----------------------------------------------------------------------------
/*!
 * \file   env_test.cc
 * \author Wayne Joubert
 * \date   Thu Oct 24 11:49:15 EDT 2019
 * \brief  Tester for env class.
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

#include "cstdlib"
#include "cstdio"
#include "cstdint"
#include "inttypes.h"
#include "string.h"
#include "math.h"

#include "gtest/gtest.h"

#include "env.hh"

enum {PROCS_MAX = TEST_PROCS_MAX};

//=============================================================================

void CEnvTest_general_() {

  int comm_rank = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank));

  const int num_proc = 20;
  COMET_INSIST(num_proc <= PROCS_MAX);

  const char options_template[] = "--num_proc_vector %i --verbosity 0 "
                                  "--num_vector_local 3 --num_field 1 "
                                  "--all2all yes --compute_method GPU "
                                  "--num_way 3 ";




  char options[1024];
  sprintf(options, options_template, num_proc);

  comet::CEnv env(options, num_proc, comm_rank);




  EXPECT_EQ(true, true);

}

//=============================================================================

TEST(CEnvTest, general) {
  CEnvTest_general_();
}

//=============================================================================

GTEST_API_ int main(int argc, char** argv) {

  ::testing::InitGoogleTest(&argc, argv);

  COMET_MPI_SAFE_CALL(MPI_Init(&argc, &argv));

  int comm_rank = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank));

  if (comm_rank != 0) {
    ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
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
