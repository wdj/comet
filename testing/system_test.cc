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
#include <math.h>

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

TEST(SystemTest,General) {

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

TEST(SystemTest,CCCSimple) {

  const int num_field = 5;
  const int num_vector_local = 2;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 2;
  env->all2all_ = GM_BOOL_FALSE;
  Env_set_compute_method(env, GM_COMPUTE_METHOD_REF);
  Env_set_num_proc(env, 1, 1);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, Env_data_type_vectors(env),
                   num_field, num_vector_local, env);

  const int G0 = 0;
  const int T0 = 1;
  const int G1 = 0;
  const int A1 = 1;

  if (Env_is_proc_active(env)) {
    int i = 0;
    GMVectors_bits2_set(vectors, i++, 0, 2*G0 + 1*T0 , env);
    GMVectors_bits2_set(vectors, i++, 0, 2*T0 + 1*T0 , env);
    GMVectors_bits2_set(vectors, i++, 0, 2*T0 + 1*T0 , env);
    GMVectors_bits2_set(vectors, i++, 0, 2*T0 + 1*T0 , env);
    GMVectors_bits2_set(vectors, i++, 0, 2*T0 + 1*T0 , env);
    i = 0;
    GMVectors_bits2_set(vectors, i++, 1, 2*G1 + 1*G1 , env);
    GMVectors_bits2_set(vectors, i++, 1, 2*A1 + 1*G1 , env);
    GMVectors_bits2_set(vectors, i++, 1, 2*G1 + 1*G1 , env);
    GMVectors_bits2_set(vectors, i++, 1, 2*G1 + 1*G1 , env);
    GMVectors_bits2_set(vectors, i++, 1, 2*G1 + 1*A1 , env);
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetrics_create(metrics, Env_data_type_metrics(env),
                   num_field, num_vector_local, env);

  gm_compute_metrics(metrics, vectors, env);

  if (Env_is_proc_active(env)) {
    const double result00 = GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 0,
                                                             env);
    const double result01 = GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 1,
                                                             env);
    const double result10 = GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 0,
                                                             env);
    const double result11 = GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 1,
                                                             env);

    printf("G G  %.5f\n", result00);
    printf("G A  %.5f\n", result01);
    printf("T G  %.5f\n", result10);
    printf("T A  %.5f\n", result11);

    const double ref00 = .196;
    const double ref01 = .000;
    const double ref10 = .588;
    const double ref11 = .312;

    const double eps = 1.e-5;

    EXPECT_EQ(GM_BOOL_TRUE, fabs(result00-ref00)<eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result01-ref01)<eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result10-ref10)<eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result11-ref11)<eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);

  GMEnv_destroy(env);
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
