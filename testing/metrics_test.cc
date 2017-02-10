/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics_test.cc
 * \author Wayne Joubert
 * \date   Fri Feb 10 17:35:20 EST 2017
 * \brief  Tester for metrics class.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gtest/gtest.h"

#include "env.hh"
#include "metrics.hh"

enum {PROCS_MAX = TEST_PROCS_MAX};

/*===========================================================================*/

void MetricsTest_3way_num_elts_local_() {


}

/*===========================================================================*/

TEST(MetricsTest, 3way_num_elts_local) {
  MetricsTest_3way_num_elts_local_();
}

/*===========================================================================*/

GTEST_API_ int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  MPI_Finalize();
  return result;
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
