//-----------------------------------------------------------------------------
/*!
 * \file   metrics_test.cc
 * \author Wayne Joubert
 * \date   Fri Feb 10 17:35:20 EST 2017
 * \brief  Tester for metrics class.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gtest/gtest.h"

#include "env.hh"
#include "metrics.hh"

enum {PROCS_MAX = TEST_PROCS_MAX};

//=============================================================================

void MetricsTest_3way_num_elts_local_() {

  int comm_rank = 0;
  int mpi_code = 0;
  mpi_code *= 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  const int num_proc_mock = 20;
  const int nvl = 12;
  const int num_stage = 2;

  const GMUInt64 num_vector = nvl * (GMUInt64)num_proc_mock;

  const GMUInt64 num_elts_local_expected = num_vector * (num_vector-1) *
                                           (num_vector-2) / 6;

  size_t num_elts_local = 0;

  const char options_template[] = "--num_proc_vector %i --verbosity 0 "
                                  "--num_vector_local %i --num_field 1 "
                                  "--all2all yes --compute_method GPU "
                                  "--num_way 3 ";

  if (comm_rank == 0) {

    for (int proc_num_mock=0; proc_num_mock<num_proc_mock; ++proc_num_mock) {
      for (int stage_num=0; stage_num<num_stage; ++stage_num) {

        char options[1024];
        sprintf(options, options_template, num_proc_mock, nvl);

        GMEnv env = GMEnv_null();
        GMEnv_create_no_comms(&env, options, NULL, num_proc_mock,
                              proc_num_mock);
        env.num_stage = num_stage;
        env.stage_num = stage_num;

        GMMetrics metrics = GMMetrics_null();
        GMMetrics_3way_num_elts_local(&metrics, nvl, &env);
        num_elts_local += metrics.num_elts_local;

        GMEnv_destroy(&env);
      } //---for
    } //---for

    printf("--------------------------------------- "
           "Metrics elements found %lu, expected %llu."
           " ---------------------------------------\n",
           num_elts_local, num_elts_local_expected);

  EXPECT_EQ(true, num_elts_local == num_elts_local_expected);
  } else {

    num_elts_local = num_elts_local_expected;

  } //---if

}

//=============================================================================

TEST(MetricsTest, 3way_num_elts_local) {
  MetricsTest_3way_num_elts_local_();
}

//=============================================================================

GTEST_API_ int main(int argc, char** argv) {

  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);

  int comm_rank = 0;
  int mpi_code = 0;
  mpi_code *= 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  if (comm_rank != 0) {
    ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
  }

  int result = RUN_ALL_TESTS();
  int result_g = 11;

  mpi_code = MPI_Allreduce(&result, &result_g,
                           1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  MPI_Finalize();
  return result_g;
}

//=============================================================================

//-----------------------------------------------------------------------------
