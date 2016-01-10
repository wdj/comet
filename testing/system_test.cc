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

_Bool compare_2runs(const char* options1, const char* options2) {

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

_Bool compare_3runs(const char* options1, const char* options2,
                    const char* options3) {

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

  size_t len3 = strlen(options3);
  char* argstring3 = (char*)malloc((len3+1)*sizeof(char));
  const char* argv3[len3+1];
  int argc3 = 0;
  strcpy(argstring3, options3);
  create_args(argstring3, &argc3, argv3);

  /*---Do runs---*/

  GMChecksum checksum1 = perform_run(argc1, argv1);
  GMChecksum checksum2 = perform_run(argc2, argv2);
  GMChecksum checksum3 = perform_run(argc3, argv3);

  /*---Need test result only on proc 0---*/

  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);
  if (proc_num != 0) {
    int i = 0;
    for (i = 0; i < GM_CHECKSUM_SIZE; ++i ) {
      checksum1.data[i] = 0;
      checksum2.data[i] = 0;
      checksum3.data[i] = 0;
    }
  }

  free(argstring1);
  free(argstring2);
  free(argstring3);

  return gm_are_checksums_equal(checksum1, checksum2) &&
         gm_are_checksums_equal(checksum1, checksum3);
}

/*===========================================================================*/

void SystemTest_czekanowski_() {

  //----------
  //---2-way, all2all no
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                      "--compute_method CPU",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                      "--compute_method GPU"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method GPU"));

  //----------
  //---2-way, all2all yes, small
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                       "--compute_method CPU ",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                       "--compute_method CPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                       "--compute_method CPU ",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                                       "--compute_method GPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                                       "--compute_method CPU ",
    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                                       "--compute_method CPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                                       "--compute_method CPU ",
    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                                       "--compute_method GPU --all2all yes"));

  //----------
  //---2-way, all2all yes, large
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method GPU"
      " --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                                        "--compute_method CPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                       "--compute_method CPU",
    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                                       "--compute_method GPU --all2all yes"));

  //----------
  //---3-way, all2all no
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                                      "--compute_method CPU "
      "--num_way 3",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                                      "--compute_method GPU --num_way 3"));
  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                        "--compute_method GPU --num_way 3"));

  //----------
  //---3-way, all2all yes, small
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                        "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                        "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                        "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                        "--compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---3-way, all2all yes, large
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                        "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                                        "--compute_method GPU --num_way 3",
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                        "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                        "--compute_method CPU --num_way 3 --all2all yes",
    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                        "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                        "--compute_method GPU --num_way 3 --all2all yes",
    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                        "--compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---num_proc_field
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_proc_field 1 --num_field 2 --num_vector_local 2 "
                        "--compute_method CPU",
    "--num_proc_vector 1 --num_proc_field 2 --num_field 2 --num_vector_local 2 "
                        "--compute_method GPU"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_proc_field 1 --num_field 2 --num_vector_local 4 "
                        "--compute_method CPU",
    "--num_proc_vector 2 --num_proc_field 2 --num_field 2 --num_vector_local 2 "
                        "--compute_method GPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_proc_field 1 --num_field 2 --num_vector_local 3 "
                        "--compute_method CPU --num_way 3",
    "--num_proc_vector 1 --num_proc_field 2 --num_field 2 --num_vector_local 3 "
                        "--compute_method GPU --num_way 3"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(
    "--num_proc_vector 1 --num_proc_field 1 --num_field 2 --num_vector_local 18"
                       " --compute_method CPU --num_way 3",
    "--num_proc_vector 3 --num_proc_field 2 --num_field 2 --num_vector_local 6 "
                       " --compute_method GPU --num_way 3 --all2all yes"));
}

/*===========================================================================*/

void SystemTest_ccc_simple_() {

  const int num_field = 5;
  const int num_vector_local = 2;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 2;
  env->all2all_ = GM_BOOL_FALSE;
  Env_set_compute_method(env, GM_COMPUTE_METHOD_GPU);
  Env_set_num_proc(env, 1, 1);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, Env_data_type_vectors(env),
                   num_field, num_vector_local, env);

  if (Env_is_proc_active(env)) {
    const int G0 = 0;
    const int T0 = 1;
    int i = 0;
    GMVectors_bits2_set(vectors, i++, 0, 2*G0 + 1*T0 , env);
    GMVectors_bits2_set(vectors, i++, 0, 2*T0 + 1*T0 , env);
    GMVectors_bits2_set(vectors, i++, 0, 2*T0 + 1*T0 , env);
    GMVectors_bits2_set(vectors, i++, 0, 2*T0 + 1*T0 , env);
    GMVectors_bits2_set(vectors, i++, 0, 2*T0 + 1*T0 , env);
    const int G1 = 0;
    const int A1 = 1;
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

void SystemTest_ccc_() {

  char options1[1024];
  char options2[1024];
  char options3[1024];
  int i;

  //----------
  //---2-way, all2all no
  //----------

  char options_template_1[] = "--metric_type ccc --verbosity %i "
    "--num_proc_vector 1 --num_field %i --num_vector_local %i "
                                       "--compute_method %s";

  sprintf(options1, options_template_1,  2,  1, 2, "REF");
  sprintf(options2, options_template_1,  2,  1, 2, "CPU");
  sprintf(options3, options_template_1,  2,  1, 2, "GPU");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_1,  2,  100, 2, "REF");
  sprintf(options2, options_template_1,  2,  100, 2, "CPU");
  sprintf(options3, options_template_1,  2,  100, 2, "GPU");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_1,  1,  100, 48, "REF");
  sprintf(options2, options_template_1,  1,  100, 48, "CPU");
  sprintf(options3, options_template_1,  1,  100, 48, "GPU");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  for (i=1; i<=200; ++i) {
    sprintf(options1, options_template_1,  0,  i, 48, "REF");
    sprintf(options2, options_template_1,  0,  i, 48, "CPU");
    sprintf(options3, options_template_1,  0,  i, 48, "GPU");
    EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));
  }

  //----------
  //---2-way, all2all yes, small
  //----------

  char options_template_2[] = "--metric_type ccc --verbosity %i "
    "--num_proc_vector %i --num_field %i --num_vector_local %i "
                                       "--compute_method %s --all2all %s";

  sprintf(options1, options_template_2,  1,  1, 1, 2, "REF", "no");
  sprintf(options2, options_template_2,  1,  1, 1, 2, "CPU", "yes");
  sprintf(options3, options_template_2,  1,  1, 1, 2, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_2,  1,  1, 1, 4, "REF", "no");
  sprintf(options2, options_template_2,  1,  2, 1, 2, "CPU", "yes");
  sprintf(options3, options_template_2,  1,  2, 1, 2, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  //----------
  //---2-way, all2all yes, large
  //----------

  sprintf(options1, options_template_2,  1,  1, 100, 48, "REF", "no");
  sprintf(options2, options_template_2,  1,  1, 100, 48, "CPU", "yes");
  sprintf(options3, options_template_2,  1,  1, 100, 48, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_2,  1,  1, 100, 48, "REF", "no");
  sprintf(options2, options_template_2,  1,  2, 100, 24, "CPU", "yes");
  sprintf(options3, options_template_2,  1,  2, 100, 24, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  //----------
  //---num_proc_field
  //----------

  char options_template_3[] = "--metric_type ccc --verbosity %i "
    "--num_proc_vector %i --num_proc_field %i "
    " --num_field %i --num_vector_local %i "
                                       "--compute_method %s --all2all %s";

  for (i=1; i<=3; ++i) {
    sprintf(options1, options_template_3,  1,  1, 1,     60, 48, "REF", "no");
    sprintf(options2, options_template_3,  1,  1, 2*i-1, 60, 48, "GPU", "yes");
    sprintf(options3, options_template_3,  1,  1, 2*i,   60, 48, "GPU", "yes");
    EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));
    sprintf(options1, options_template_3,  1,  1, 1, 60, 48, "REF", "no");
    sprintf(options2, options_template_3,  1,  1, i, 60, 48, "GPU", "yes");
    sprintf(options3, options_template_3,  1,  2, i, 60, 24, "GPU", "yes");
    EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));
  }
}

/*===========================================================================*/

TEST(SystemTest,czekanowski) {
  SystemTest_czekanowski_();
}

TEST(SystemTest,ccc_simple) {
  SystemTest_ccc_simple_();
}

TEST(SystemTest,ccc) {
  SystemTest_ccc_();
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