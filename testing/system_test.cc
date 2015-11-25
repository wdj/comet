/*---------------------------------------------------------------------------*/
/*!
 * \file   system_test.cc
 * \author Wayne Joubert
 * \date   Fri Nov  6 18:18:21 EST 2015
 * \brief  Perform high-level system tests.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

#include "gtest/gtest.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics.h"

#include "driver_utils.h"

/*===========================================================================*/

_Bool compare_runs(const char* options1, const char* options2) {

  /*---Convert options strings to args---*/

  size_t len1 = strlen(options1);
  char* argstring1 = (char*)malloc((len1+1)*sizeof(char));
  const char* argv1[len1+1];
  int argc1 = 0;
  strcpy(argstring1, options1);
  create_args(argstring1, &argc1, argv1);

  size_t len2 = strlen(options2);
  char* argstring2 = (char*)malloc((len2+1)*sizeof(char));
  const char* argv2[len2+1];
  int argc2 = 0;
  strcpy(argstring2, options2);
  create_args(argstring2, &argc2, argv2);

  /*---Do runs---*/

  GMChecksum checksum1 = perform_run(argc1, argv1);
  GMChecksum checksum2 = perform_run(argc2, argv2);

  /*---Need test result only on proc 0---*/

  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);
  if (proc_num != 0) {
    int i = 0;
    for (i = 0; i < GM_CHECKSUM_SIZE; ++i ) {
      checksum1.data[i] = 0;
      checksum2.data[i] = 0;
    }
  }

  free(argstring1);
  free(argstring2);

  return gm_are_checksums_equal(checksum1, checksum2);
}

/*===========================================================================*/

TEST(SystemTest,One) {

  //----------
  //---2-way, all2all no
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                      "--compute_method CPU",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                      "--compute_method GPU"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method GPU"));

  //----------
  //---2-way, all2all yes, small
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                       "--compute_method CPU ",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                       "--compute_method CPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                       "--compute_method CPU ",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                       "--compute_method GPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                                       "--compute_method CPU ",
    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                                       "--compute_method CPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                                       "--compute_method CPU ",
    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                                       "--compute_method GPU --all2all yes"));

  //----------
  //---2-way, all2all yes, large
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method GPU"
      " --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                                        "--compute_method CPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                                       "--compute_method GPU --all2all yes"));

  //----------
  //---3-way, all2all no
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                                      "--compute_method CPU "
      "--num_way 3",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                                      "--compute_method GPU --num_way 3"));
  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                        "--compute_method GPU --num_way 3"));

  //----------
  //---3-way, all2all yes, small
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                        "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                        "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                        "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                        "--compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---3-way, all2all yes, large
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                        "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                        "--compute_method GPU --num_way 3",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                        "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                        "--compute_method CPU --num_way 3 --all2all yes",
    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                        "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                        "--compute_method GPU --num_way 3 --all2all yes",
    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                        "--compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---num_proc_field
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_proc_field 1 --num_field 2 --num_vector_local 2 "
                        "--compute_method CPU",
    "--num_proc_vector 1 --num_proc_field 2 --num_field 2 --num_vector_local 2 "
                        "--compute_method GPU"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_proc_field 1 --num_field 2 --num_vector_local 4 "
                        "--compute_method CPU",
    "--num_proc_vector 2 --num_proc_field 2 --num_field 2 --num_vector_local 2 "
                        "--compute_method GPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_proc_field 1 --num_field 2 --num_vector_local 3 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_proc_field 2 --num_field 2 --num_vector_local 3 "
                        "--compute_method GPU --num_way 3"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc_vector 1 --num_proc_field 1 --num_field 2 --num_vector_local 18"
                       " --compute_method CPU --num_way 3",
    "--num_proc_vector 3 --num_proc_field 2 --num_field 2 --num_vector_local 6 "
                       " --compute_method GPU --num_way 3 --all2all yes"));
}

/*===========================================================================*/

GTEST_API_ int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  MPI_Finalize();
  return result;
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
