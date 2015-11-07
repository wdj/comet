/*---------------------------------------------------------------------------*/
/*!
 * \file   basic_test.c
 * \author Wayne Joubert
 * \date   Fri Nov  6 18:18:21 EST 2015
 * \brief  Perform basic tests.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>

#include "gtest/gtest.h"

#include "env.h"

/*===========================================================================*/

TEST(BasicTest,One) {
  EXPECT_EQ(1, 1);
}

/*===========================================================================*/

GTEST_API_ int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  testing::InitGoogleTest(&argc, argv);

  GMEnv env = GMEnv_null();
  GMEnv_create(&env);

  if (env.proc_num == 0) {
    printf("///////// Running main()  /////////\n");
  }

  MPI_Finalize();

  return RUN_ALL_TESTS();
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
