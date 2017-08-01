/*---------------------------------------------------------------------------*/
/*!
 * \file   driver_test.cc
 * \author Wayne Joubert
 * \date   Fri Nov  6 18:18:21 EST 2015
 * \brief  Tester for driver.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gtest/gtest.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics.hh"

#include "driver.hh"

enum {PROCS_MAX = TEST_PROCS_MAX};

/*===========================================================================*/

_Bool compare_2runs(const char* options1, const char* options2) {
  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  /*---Do runs---*/

  if (proc_num == 0) {
    printf("%s\n", options1);
  }
  GMChecksum checksum1 = perform_run(options1);

  if (proc_num == 0) {
    printf("%s\n", options2);
  }
  GMChecksum checksum2 = perform_run(options2);

  /*---Need test result only on proc 0---*/

  const _Bool is_passed = proc_num != 0 ? true :
                          gm_are_checksums_equal(checksum1, checksum2);

  return is_passed;
}

/*===========================================================================*/

_Bool compare_3runs(const char* options1,
                    const char* options2,
                    const char* options3) {
  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  /*---Do runs---*/

  if (proc_num == 0) {
    printf("%s\n", options1);
  }
  GMChecksum checksum1 = perform_run(options1);

  if (proc_num == 0) {
    printf("%s\n", options2);
  }
  GMChecksum checksum2 = perform_run(options2);

  if (proc_num == 0) {
    printf("%s\n", options3);
  }
  GMChecksum checksum3 = perform_run(options3);

  /*---Need test result only on proc 0---*/

  const _Bool is_passed = proc_num != 0 ? true :
                          gm_are_checksums_equal(checksum1, checksum2) &&
                          gm_are_checksums_equal(checksum1, checksum3);
  return is_passed;
}

/*===========================================================================*/

void test_2runs(const char* options1,
                const char* options2) {
  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(options1, options2));
}

/*===========================================================================*/

void DriverTest_czekanowski_() {

  //----------
  //---2-way, all2all no
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU"));

  //----------
  //---2-way, all2all yes, small
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU ",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU ",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                    "--compute_method CPU ",
                    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 4 "
                    "--compute_method CPU ",
                    "--num_proc_vector 2 --num_field 1 --num_vector_local 2 "
                    "--compute_method GPU --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_3runs("--num_proc_vector 1 --num_field 2 --num_vector 5 "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 2 --num_field 2 --num_vector 5 "
                    "--compute_method CPU --all2all yes",
                    "--num_proc_vector 2 --num_field 2 --num_vector 5 "
                    "--compute_method GPU --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_proc_field 1 "
                    "--num_field 7 --num_vector 5 "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_field 3 "
                    "--num_field 7 --num_vector 5 "
                    "--compute_method GPU --all2all yes"));

  //----------
  //---2-way, all2all yes, large
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU"
                    " --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                    "--compute_method CPU --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU",
                    "--num_proc_vector 2 --num_field 100 --num_vector_local 24 "
                    "--compute_method GPU --all2all yes"));

  //----------
  //---3-way, all2all no
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU "
                    "--num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method GPU --num_way 3"));
  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3"));

  //----------
  //---3-way, all2all yes, small
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 1 --num_vector_local 3 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 1 --num_vector_local 18 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 3 --num_field 1 --num_vector_local 6 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_3runs("--num_proc_vector 1 --num_field 2 --num_vector 16 "
                    "--compute_method REF --all2all yes --num_way 3",
                    "--num_proc_vector 3 --num_field 2 --num_vector 16 "
                    "--compute_method CPU --all2all yes --num_way 3",
                    "--num_proc_vector 3 --num_field 2 --num_vector 16 "
                    "--compute_method GPU --all2all yes --num_way 3"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_proc_field 1 "
                    "--num_field 13 --num_vector 6 "
                    "--compute_method REF --all2all yes --num_way 3",
                    "--num_proc_vector 1 --num_proc_field 5 "
                    "--num_field 13 --num_vector 6 "
                    "--compute_method GPU --all2all yes --num_way 3"));

  //----------
  //---3-way, all2all yes, large
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3",
                    "--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method CPU --num_way 3 --all2all yes",
                    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                    "--compute_method CPU --num_way 3 --all2all yes"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 100 --num_vector_local 48 "
                    "--compute_method GPU --num_way 3 --all2all yes",
                    "--num_proc_vector 4 --num_field 100 --num_vector_local 12 "
                    "--compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---num_proc_field
  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                        "1 --num_field 2 --num_vector_local 2 "
                                        "--compute_method CPU",
                                        "--num_proc_vector 1 --num_proc_field "
                                        "2 --num_field 2 --num_vector_local 2 "
                                        "--compute_method GPU"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                        "1 --num_field 2 --num_vector_local 4 "
                                        "--compute_method CPU",
                                        "--num_proc_vector 2 --num_proc_field "
                                        "2 --num_field 2 --num_vector_local 2 "
                                        "--compute_method GPU --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs("--num_proc_vector 1 --num_proc_field "
                                        "1 --num_field 2 --num_vector_local 3 "
                                        "--compute_method CPU --num_way 3",
                                        "--num_proc_vector 1 --num_proc_field "
                                        "2 --num_field 2 --num_vector_local 3 "
                                        "--compute_method GPU --num_way 3"));

  EXPECT_EQ(GM_BOOL_TRUE,
            compare_2runs("--num_proc_vector 1 --num_proc_field 1 --num_field "
                          "2 --num_vector_local 18"
                          " --compute_method CPU --num_way 3",
                          "--num_proc_vector 3 --num_proc_field 2 --num_field "
                          "2 --num_vector_local 6 "
                          " --compute_method GPU --num_way 3 --all2all yes"));

  //----------
  //---num_repl, 2-way
  //----------

  char options1[1024];
  char options2[1024];

  char options_template_1[] =
      "--metric_type czekanowski "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_stage %i";

  int num_vector_local = 0;
  int num_proc_vector = 0;
  int num_proc_repl = 0;
  int gpu = 0;
  int num_stage = 1;

  for (gpu=0; gpu<=1; ++gpu) {
    for (num_vector_local=4; num_vector_local<=5; ++num_vector_local) {
      for (num_proc_vector=1; num_proc_vector<=6; ++num_proc_vector) {
        for (num_proc_repl=2; num_proc_repl<=6; ++num_proc_repl) {
          const int num_proc_field = gpu ? 2 : 1;
          if (num_proc_vector * num_proc_field * num_proc_repl > PROCS_MAX) {
            continue;
          }
          const int num_way = 2;
          sprintf(options1, options_template_1,
                  num_vector_local*num_proc_vector,
                  "GPU", 1, 1,
                  1, num_way, 1);
          sprintf(options2, options_template_1, num_vector_local,
                  gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                  num_proc_field, num_way, num_stage);
          test_2runs(options1, options2);
        }
      }
    }
  }

  //----------
  //---num_repl, num_stage, 3-way
  //----------

  for (gpu=0; gpu<=1; ++gpu) {
    for (num_vector_local=6; num_vector_local<=18; num_vector_local+=12) {
      for (num_proc_vector=1; num_proc_vector<=6; ++num_proc_vector) {
        for (num_proc_repl=2; num_proc_repl<=6; ++num_proc_repl) {
          for (num_stage=1; num_stage<=6; num_stage+=4) {
          const int num_proc_field = gpu ? 2 : 1;
            if (num_proc_vector * num_proc_field * num_proc_repl > PROCS_MAX) {
              continue;
            }
            const int num_way = 3;
            sprintf(options1, options_template_1,
                    num_vector_local*num_proc_vector,
                    "GPU", 1, 1,
                   1, num_way, 1);
            sprintf(options2, options_template_1, num_vector_local,
                    gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                    num_proc_field, num_way, num_stage);
            test_2runs(options1, options2);
          }
        }
      }
    }
  }

  //----------
  //---num_phase, 2-way
  //----------

  int num_phase = 0;

  for (num_proc_vector=1; num_proc_vector<=8; ++num_proc_vector) {
    for (num_proc_repl=1; num_proc_repl<=8; ++num_proc_repl) {
      for (num_phase=2; num_phase<=8; ++num_phase) {
        if (!(num_phase <= 1 + num_proc_vector/2)) {
          continue;
        }
        char options_template[] =
          "--metric_type czekanowski "
          "--num_field 7 --num_vector 12 --compute_method GPU --all2all yes "
          "--num_proc_vector %i --num_proc_repl %i --num_phase %i "
          "--num_way 2";
        sprintf(options1, options_template, num_proc_vector, num_proc_repl,
                num_phase);
        sprintf(options2, options_template, 1, 1, 1);
        test_2runs(options1, options2);
      }
    }
  }

  //----------
  //---file output, 2-way
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_proc_vector 1 "
                    "--num_field 7 --num_vector 10 "
                    "--num_way 2 --metric_type czekanowski "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_vector 3 "
                    "--num_field 7 --num_vector 10 "
                    "--num_way 2 --metric_type czekanowski "
                    "--compute_method GPU --all2all yes "
                    "--verbosity 2 "
                    "--output_file_stub test_czek_2way"));

  //----------
  //---file output, 3-way
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_proc_vector 1 "
                    "--num_field 7 --num_vector 18 "
                    "--num_way 3 --metric_type czekanowski "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_vector 3 "
                    "--num_field 7 --num_vector 18 "
                    "--num_way 3 --metric_type czekanowski "
                    "--compute_method GPU --all2all yes "
                    "--verbosity 2 "
                    "--output_file_stub test_czek_3way"));

  //----------
  //---Misc options
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 30 --num_vector 3 "
                    "--verbosity 2 --all2all yes "
                    "--compute_method GPU",
                    "--num_proc_vector 1 --num_field 30 --num_vector 3 "
                    "--verbosity 2 --all2all yes "
                    "--compute_method GPU --threshold .65"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_field 3 --num_vector 3 "
                    "--compute_method GPU",
                    "--num_proc_vector 1 --num_field 3 --num_vector 3 "
                    "--compute_method GPU --checksum no"));

}

/*===========================================================================*/

void DriverTest_ccc2_simple_compute_method(int compute_method) {
  const int num_field = 5;
  const int num_vector_local = 2;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 2;
  env->all2all_ = GM_BOOL_FALSE;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), num_field,
                   num_field, num_vector_local, env);

  if (GMEnv_is_proc_active(env)) {
    {
      const int G = 0;
      const int T = 1;
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
    }
    {
      const int G = 0;
      const int A = 1;
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * G, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * G, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * G, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * G, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * G + 1 * A, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), num_field, num_field,
                   num_vector_local, num_vector_local, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result00 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 0, env);
    const double result01 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 1, env);
    const double result10 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 0, env);
    const double result11 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 1, env);

    printf("G G  %.5f\n", result00);
    printf("G A  %.5f\n", result01);
    printf("T G  %.5f\n", result10);
    printf("T A  %.5f\n", result11);
    printf("\n");

    const double ref00 = .196;
    const double ref01 = .000;
    const double ref10 = .588;
    const double ref11 = .312;

    const double eps = 1.e-5;

    EXPECT_EQ(GM_BOOL_TRUE, fabs(result00 - ref00) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result01 - ref01) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result10 - ref10) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result11 - ref11) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);

  GMEnv_destroy(env);
}

/*===========================================================================*/

void DriverTest_ccc2_simple_() {
  DriverTest_ccc2_simple_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_ccc2_simple_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_ccc2_simple_compute_method(GM_COMPUTE_METHOD_GPU);
}

/*===========================================================================*/

void DriverTest_ccc2_simple_sparse_compute_method(int compute_method) {
  const int num_field = 5;
  const int num_vector_local = 2;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 2;
  env->all2all_ = GM_BOOL_FALSE;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);
  env->sparse = true;

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), num_field,
                   num_field, num_vector_local, env);

  if (GMEnv_is_proc_active(env)) {
    const int UN = 2 * 1 + 1 * 0;
    {
      const int G = 0;
      const int T = 1;
    //const int GG =  2 * G + 1 * G;
      const int GT =  2 * G + 1 * T;
    //const int TG =  2 * G + 1 * T;
      const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, GT, env);
      GMVectors_bits2_set(vectors, f++, i, TT, env);
      GMVectors_bits2_set(vectors, f++, i, TT, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
    }
    {
      const int G = 0;
      const int A = 1;
      const int GG =  2 * G + 1 * G;
      const int GA =  2 * G + 1 * A;
      const int AG =  2 * G + 1 * A;
    //const int AA =  2 * Q + 1 * A;
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, GG, env);
      GMVectors_bits2_set(vectors, f++, i, AG, env);
      GMVectors_bits2_set(vectors, f++, i, GG, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, GA, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), num_field, num_field,
                   num_vector_local, num_vector_local, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result00 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 0, env);
    const double result01 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 0, 1, env);
    const double result10 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 0, env);
    const double result11 =
        GMMetrics_ccc_get_from_index_2(metrics, 0, 1, 1, env);

    printf("G G  %.5f\n", result00);
    printf("G A  %.5f\n", result01);
    printf("T G  %.5f\n", result10);
    printf("T A  %.5f\n", result11);
    printf("\n");

    const double s0_0 = 1;
    const double s0_1 = 5;

    const double s1_0 = 6;
    const double s1_1 = 2;

    const double c0 = s0_0 + s0_1;
    const double c1 = s1_0 + s1_1;

    const double f0_0 = s0_0 / c0;
    const double f0_1 = s0_1 / c0;

    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;

    const double r0_00 = 2;
    const double r0_01 = 0;
    const double r0_10 = 2;
    const double r0_11 = 0;

    const double r1_00 = 0;
    const double r1_01 = 0;
    const double r1_10 = 2;
    const double r1_11 = 2;

    const double r2_00 = 0;
    const double r2_01 = 0;
    const double r2_10 = 4;
    const double r2_11 = 0;

    const double r_00 = r0_00 + r1_00 + r2_00;
    const double r_01 = r0_01 + r1_01 + r2_01;
    const double r_10 = r0_10 + r1_10 + r2_10;
    const double r_11 = r0_11 + r1_11 + r2_11;

    const double c = r_00 + r_01 + r_10 + r_11;

    const double f_00 = r_00 / c;
    const double f_01 = r_01 / c;
    const double f_10 = r_10 / c;
    const double f_11 = r_11 / c;

    const double fm = 9 / (double) 2;
    const double cp = 2 / (double) 3;

    const double ref00 = fm * f_00 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_0 );
    const double ref01 = fm * f_01 * ( 1 - cp * f0_0 ) * ( 1 - cp * f1_1 );
    const double ref10 = fm * f_10 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_0 );
    const double ref11 = fm * f_11 * ( 1 - cp * f0_1 ) * ( 1 - cp * f1_1 );

    const double eps = 1.e-5;

    EXPECT_EQ(GM_BOOL_TRUE, fabs(result00 - ref00) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result01 - ref01) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result10 - ref10) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result11 - ref11) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);

  GMEnv_destroy(env);
}

/*===========================================================================*/

void DriverTest_ccc2_simple_sparse_() {
  DriverTest_ccc2_simple_sparse_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_ccc2_simple_sparse_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_ccc2_simple_sparse_compute_method(GM_COMPUTE_METHOD_GPU);
}

/*===========================================================================*/

void DriverTest_ccc3_simple_compute_method(int compute_method) {
  const int num_field = 10;
  const int num_vector_local = 3;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 3;
  env->all2all_ = GM_BOOL_TRUE;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), num_field,
                   num_field, num_vector_local, env);

  if (GMEnv_is_proc_active(env)) {
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
    }
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
    }
    {
      const int A = 0;
      const int T = 1;
      int f = 0;
      const int i = 2;
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * T + 1 * T, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
      GMVectors_bits2_set(vectors, f++, i, 2 * A + 1 * A, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), num_field, num_field,
                   num_vector_local, num_vector_local, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result000 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 0, 0, env);
    const double result001 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 0, 1, env);
    const double result010 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 1, 0, env);
    const double result011 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 1, 1, env);
    const double result100 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 0, 0, env);
    const double result101 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 0, 1, env);
    const double result110 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 1, 0, env);
    const double result111 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 1, 1, env);

    printf("A A A  %.8f\n", result000);
    printf("A A T  %.5f\n", result001);
    printf("A T A  %.8f\n", result010);
    printf("A T T  %.8f\n", result011);
    printf("T A A  %.8f\n", result100);
    printf("T A T  %.8f\n", result101);
    printf("T T A  %.8f\n", result110);
    printf("T T T  %.8f\n", result111);
    printf("\n");

    const double ref000 = .055;
    const double ref001 = .016;
    const double ref010 = .016;
    const double ref011 = .030;
    const double ref100 = .039;
    const double ref101 = .008;
    const double ref110 = .008;
    const double ref111 = .015;

    const double eps = 1.e-3;

    EXPECT_EQ(GM_BOOL_TRUE, fabs(result000 - ref000) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result001 - ref001) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result010 - ref010) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result011 - ref011) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result100 - ref100) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result101 - ref101) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result110 - ref110) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result111 - ref111) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);

  GMEnv_destroy(env);
}

/*===========================================================================*/

void DriverTest_ccc3_simple_() {
  DriverTest_ccc3_simple_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_ccc3_simple_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_ccc3_simple_compute_method(GM_COMPUTE_METHOD_GPU);
}

/*===========================================================================*/

void DriverTest_ccc3_simple_sparse_compute_method(int compute_method) {
  const int num_field = 10;
  const int num_vector_local = 3;

  GMEnv env_value = GMEnv_null();
  GMEnv* env = &env_value;
  GMEnv_create(env, MPI_COMM_WORLD, 0, NULL, NULL);
  env->metric_type_ = GM_METRIC_TYPE_CCC;
  env->num_way_ = 3;
  env->all2all_ = GM_BOOL_TRUE;
  GMEnv_set_compute_method(env, compute_method);
  GMEnv_set_num_proc(env, 1, 1, 1);
  env->sparse = true;

  GMVectors vectors_value = GMVectors_null();
  GMVectors* vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), num_field,
                   num_field, num_vector_local, env);

  if (GMEnv_is_proc_active(env)) {
    const int UN = 2 * 1 + 1 * 0;
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T;
      const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 0;
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, TT, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
    }
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T;
    //const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 1;
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
    }
    {
      const int A = 0;
      const int T = 1;
      const int AA =  2 * A + 1 * A;
      const int AT =  2 * A + 1 * T;
    //const int TA =  2 * A + 1 * T;
    //const int TT =  2 * T + 1 * T;
      int f = 0;
      const int i = 2;
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AT, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, AA, env);
      GMVectors_bits2_set(vectors, f++, i, UN, env);
    }
  }

  GMMetrics metrics_value = GMMetrics_null();
  GMMetrics* metrics = &metrics_value;
  GMMetrics_create(metrics, GMEnv_data_type_metrics(env), num_field, num_field,
                   num_vector_local, num_vector_local, env);

  gm_compute_metrics(metrics, vectors, env);

  if (GMEnv_is_proc_active(env)) {
    const double result000 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 0, 0, env);
    const double result001 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 0, 1, env);
    const double result010 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 1, 0, env);
    const double result011 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 0, 1, 1, env);
    const double result100 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 0, 0, env);
    const double result101 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 0, 1, env);
    const double result110 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 1, 0, env);
    const double result111 =
        GMMetrics_ccc_get_from_index_3(metrics, 0, 1, 1, 1, env);

    printf("A A A  %.8f\n", result000);
    printf("A A T  %.5f\n", result001);
    printf("A T A  %.8f\n", result010);
    printf("A T T  %.8f\n", result011);
    printf("T A A  %.8f\n", result100);
    printf("T A T  %.8f\n", result101);
    printf("T T A  %.8f\n", result110);
    printf("T T T  %.8f\n", result111);
    printf("\n");

    const double s0_0 = 9;
    const double s0_1 = 5;
    const double c0 = s0_0 + s0_1;

    const double s1_0 = 12;
    const double s1_1 = 4;
    const double c1 = s1_0 + s1_1;

    const double s2_0 = 14;
    const double s2_1 = 2;
    const double c2 = s2_0 + s2_1;

    const double f0_0 = s0_0 / c0;
    const double f0_1 = s0_1 / c0;
    const double f1_0 = s1_0 / c1;
    const double f1_1 = s1_1 / c1;
    const double f2_0 = s2_0 / c2;
    const double f2_1 = s2_1 / c2;

    const double r0_000 = 4;
    const double r0_001 = 0;
    const double r0_010 = 4;
    const double r0_011 = 0;
    const double r0_100 = 0;
    const double r0_101 = 0;
    const double r0_110 = 0;
    const double r0_111 = 0;

    const double r1_000 = 4;
    const double r1_001 = 0;
    const double r1_010 = 0;
    const double r1_011 = 0;
    const double r1_100 = 4;
    const double r1_101 = 0;
    const double r1_110 = 0;
    const double r1_111 = 0;

    const double r2_000 = 0;
    const double r2_001 = 0;
    const double r2_010 = 0;
    const double r2_011 = 0;
    const double r2_100 = 8;
    const double r2_101 = 0;
    const double r2_110 = 0;
    const double r2_111 = 0;

    const double r4_000 = 1;
    const double r4_001 = 1;
    const double r4_010 = 1;
    const double r4_011 = 1;
    const double r4_100 = 1;
    const double r4_101 = 1;
    const double r4_110 = 1;
    const double r4_111 = 1;

    const double r5_000 = 4;
    const double r5_001 = 0;
    const double r5_010 = 0;
    const double r5_011 = 0;
    const double r5_100 = 4;
    const double r5_101 = 0;
    const double r5_110 = 0;
    const double r5_111 = 0;

    const double r_000 = r0_000 + r1_000 + r2_000 + r4_000 + r5_000;
    const double r_001 = r0_001 + r1_001 + r2_001 + r4_001 + r5_001;
    const double r_010 = r0_010 + r1_010 + r2_010 + r4_010 + r5_010;
    const double r_011 = r0_011 + r1_011 + r2_011 + r4_011 + r5_011;
    const double r_100 = r0_100 + r1_100 + r2_100 + r4_100 + r5_100;
    const double r_101 = r0_101 + r1_101 + r2_101 + r4_101 + r5_101;
    const double r_110 = r0_110 + r1_110 + r2_110 + r4_110 + r5_110;
    const double r_111 = r0_111 + r1_111 + r2_111 + r4_111 + r5_111;

    const double c = r_000 + r_001 + r_010 + r_011 +
                     r_100 + r_101 + r_110 + r_111;

    const double f_000 = r_000 / c;
    const double f_001 = r_001 / c;
    const double f_010 = r_010 / c;
    const double f_011 = r_011 / c;
    const double f_100 = r_100 / c;
    const double f_101 = r_101 / c;
    const double f_110 = r_110 / c;
    const double f_111 = r_111 / c;

    const double fm = 1;
    const double cp = 2 / (double) 3;

    const double ref000 = fm * f_000 * (1-cp*f0_0) * (1-cp*f1_0) * (1-cp*f2_0);
    const double ref001 = fm * f_001 * (1-cp*f0_0) * (1-cp*f1_0) * (1-cp*f2_1);
    const double ref010 = fm * f_010 * (1-cp*f0_0) * (1-cp*f1_1) * (1-cp*f2_0);
    const double ref011 = fm * f_011 * (1-cp*f0_0) * (1-cp*f1_1) * (1-cp*f2_1);
    const double ref100 = fm * f_100 * (1-cp*f0_1) * (1-cp*f1_0) * (1-cp*f2_0);
    const double ref101 = fm * f_101 * (1-cp*f0_1) * (1-cp*f1_0) * (1-cp*f2_1);
    const double ref110 = fm * f_110 * (1-cp*f0_1) * (1-cp*f1_1) * (1-cp*f2_0);
    const double ref111 = fm * f_111 * (1-cp*f0_1) * (1-cp*f1_1) * (1-cp*f2_1);

    //const double ref000 = .055;
    //const double ref001 = .016;
    //const double ref010 = .016;
    //const double ref011 = .030;
    //const double ref100 = .039;
    //const double ref101 = .008;
    //const double ref110 = .008;
    //const double ref111 = .015;

    const double eps = 1.e-3;

    EXPECT_EQ(GM_BOOL_TRUE, fabs(result000 - ref000) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result001 - ref001) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result010 - ref010) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result011 - ref011) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result100 - ref100) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result101 - ref101) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result110 - ref110) < eps);
    EXPECT_EQ(GM_BOOL_TRUE, fabs(result111 - ref111) < eps);
  }

  GMMetrics_destroy(metrics, env);
  GMVectors_destroy(vectors, env);

  GMEnv_destroy(env);
}

/*===========================================================================*/

void DriverTest_ccc3_simple_sparse_() {
  DriverTest_ccc3_simple_sparse_compute_method(GM_COMPUTE_METHOD_REF);
  DriverTest_ccc3_simple_sparse_compute_method(GM_COMPUTE_METHOD_CPU);
  DriverTest_ccc3_simple_sparse_compute_method(GM_COMPUTE_METHOD_GPU);
}

/*===========================================================================*/

void DriverTest_ccc_() {
  char options1[1024];
  char options2[1024];
  char options3[1024];
  int i = 0;

  //----------
  //---2-way, all2all no
  //----------

  char options_template_1[] =
      "--metric_type ccc --verbosity %i "
      "--num_proc_vector 1 --num_field %i --num_vector_local %i "
      "--compute_method %s";

  sprintf(options1, options_template_1, 2, 1, 2, "REF");
  sprintf(options2, options_template_1, 2, 1, 2, "CPU");
  sprintf(options3, options_template_1, 2, 1, 2, "GPU");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_1, 2, 100, 2, "REF");
  sprintf(options2, options_template_1, 2, 100, 2, "CPU");
  sprintf(options3, options_template_1, 2, 100, 2, "GPU");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_1, 1, 100, 48, "REF");
  sprintf(options2, options_template_1, 1, 100, 48, "CPU");
  sprintf(options3, options_template_1, 1, 100, 48, "GPU");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  for (i = 1; i <= 100; ++i) {
    if ((i+3) % 32 > 6) {
      continue;
    }
    sprintf(options1, options_template_1, 0, i, 48, "REF");
    sprintf(options2, options_template_1, 0, i, 48, "CPU");
    sprintf(options3, options_template_1, 0, i, 48, "GPU");
    EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));
  }

  //----------
  //---2-way, all2all yes, small
  //----------

  char options_template_2[] =
      "--metric_type ccc --verbosity %i "
      "--num_proc_vector %i --num_field %i --num_vector_local %i "
      "--compute_method %s --all2all %s";

  sprintf(options1, options_template_2, 2, 1, 1, 2, "REF", "no");
  sprintf(options2, options_template_2, 2, 1, 1, 2, "CPU", "yes");
  sprintf(options3, options_template_2, 2, 1, 1, 2, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_2, 1, 1, 1, 4, "REF", "no");
  sprintf(options2, options_template_2, 1, 2, 1, 2, "CPU", "yes");
  sprintf(options3, options_template_2, 1, 2, 1, 2, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_3runs("--num_proc_vector 1 --num_field 2 --num_vector 5 "
                    "--compute_method REF --all2all yes --metric_type ccc",
                    "--num_proc_vector 2 --num_field 2 --num_vector 5 "
                    "--compute_method CPU --all2all yes --metric_type ccc",
                    "--num_proc_vector 2 --num_field 2 --num_vector 5 "
                    "--compute_method GPU --all2all yes --metric_type ccc"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_proc_field 1 "
                    "--num_field 7 --num_vector 2 "
                    "--compute_method REF --all2all yes --metric_type ccc",
                    "--num_proc_vector 1 --num_proc_field 3 "
                    "--num_field 7 --num_vector 2 "
                    "--compute_method GPU --all2all yes --metric_type ccc"));

  //----------
  //---2-way, all2all yes, large
  //----------

  sprintf(options1, options_template_2, 1, 1, 100, 48, "REF", "no");
  sprintf(options2, options_template_2, 1, 1, 100, 48, "CPU", "yes");
  sprintf(options3, options_template_2, 1, 1, 100, 48, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_2, 1, 1, 100, 48, "REF", "no");
  sprintf(options2, options_template_2, 1, 2, 100, 24, "CPU", "yes");
  sprintf(options3, options_template_2, 1, 2, 100, 24, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  //----------
  //---3-way, all2all no
  //----------

  char options_template_3[] =
                    "--metric_type ccc --verbosity %i "
                    "--num_proc_vector 1 --num_field %i --num_vector_local %i "
                    "--compute_method %s --num_way 3";

  sprintf(options1, options_template_3, 2, 1, 3, "REF");
  sprintf(options2, options_template_3, 2, 1, 3, "CPU");
  sprintf(options3, options_template_3, 2, 1, 3, "GPU");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_3, 1, 100, 48, "REF");
  sprintf(options2, options_template_3, 1, 100, 48, "CPU");
  sprintf(options3, options_template_3, 1, 100, 48, "GPU");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  for (i = 1; i <= 100; ++i) {
    if ((i+3) % 32 > 6) {
      continue;
    }
    sprintf(options1, options_template_3, 0, i, 24, "REF");
    sprintf(options2, options_template_3, 0, i, 24, "CPU");
    sprintf(options3, options_template_3, 0, i, 24, "GPU");
    EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));
  }

  //----------
  //---2-way, all2all yes, small
  //----------

  char options_template_4[] =
      "--metric_type ccc --verbosity %i "
      "--num_proc_vector %i --num_field %i --num_vector_local %i "
      "--compute_method %s --all2all %s --num_way 3";

  sprintf(options1, options_template_4, 2, 1, 1, 3, "REF", "no");
  sprintf(options2, options_template_4, 2, 1, 1, 3, "CPU", "yes");
  sprintf(options3, options_template_4, 2, 1, 1, 3, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_4, 1, 1, 1, 6, "REF", "no");
  sprintf(options2, options_template_4, 1, 2, 1, 3, "CPU", "yes");
  sprintf(options3, options_template_4, 1, 2, 1, 3, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  //----------
  //---2-way, all2all yes, large
  //----------

  sprintf(options1, options_template_4, 1, 1, 100, 48, "REF", "no");
  sprintf(options2, options_template_4, 1, 1, 100, 48, "CPU", "yes");
  sprintf(options3, options_template_4, 1, 1, 100, 48, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

  sprintf(options1, options_template_4, 1, 1, 100, 48, "REF", "no");
  sprintf(options2, options_template_4, 1, 4, 100, 12, "CPU", "yes");
  sprintf(options3, options_template_4, 1, 4, 100, 12, "GPU", "yes");
  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));

//  sprintf(options1, options_template_4, 1, 3, 1, 6, "REF", "yes");
//  sprintf(options2, options_template_4, 1, 1, 1, 18, "GPU", "yes");
//  sprintf(options3, options_template_4, 1, 3, 1, 6, "GPU", "yes");
//  EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));
  //EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(options1, options2));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_3runs("--num_proc_vector 1 --num_field 2 --num_vector 16 "
                    "--compute_method REF --all2all yes --num_way 3 "
                    "--metric_type ccc",
                    "--num_proc_vector 3 --num_field 2 --num_vector 16 "
                    "--compute_method CPU --all2all yes --num_way 3 "
                    "--metric_type ccc",
                    "--num_proc_vector 3 --num_field 2 --num_vector 16 "
                    "--compute_method GPU --all2all yes --num_way 3 "
                    "--metric_type ccc"));

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_proc_field 1 "
                    "--num_field 13 --num_vector 6 "
                    "--compute_method REF --all2all yes --num_way 3 "
                    "--metric_type ccc",
                    "--num_proc_vector 1 --num_proc_field 5 "
                    "--num_field 13 --num_vector 6 "
                    "--compute_method GPU --all2all yes --num_way 3 "
                    "--metric_type ccc"));

  //----------
  //---num_proc_field
  //----------

  char options_template_5[] =
      "--metric_type ccc --verbosity %i "
      "--num_proc_vector %i --num_proc_field %i "
      " --num_field %i --num_vector_local %i "
      "--compute_method %s --all2all %s";

  for (i = 1; i <= 3; ++i) {
    sprintf(options1, options_template_5, 1, 1, 1, 60, 48, "REF", "no");
    sprintf(options2, options_template_5, 1, 1, 2 * i - 1, 60, 48, "GPU",
            "yes");
    sprintf(options3, options_template_5, 1, 1, 2 * i, 60, 48, "GPU", "yes");
    EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));
    sprintf(options1, options_template_5, 1, 1, 1, 60, 48, "REF", "no");
    sprintf(options2, options_template_5, 1, 1, i, 60, 48, "GPU", "yes");
    sprintf(options3, options_template_5, 1, 2, i, 60, 24, "GPU", "yes");
    EXPECT_EQ(GM_BOOL_TRUE, compare_3runs(options1, options2, options3));
  }

  //----------
  //---num_repl, 2-way
  //----------

  char options_template_10[] =
      "--metric_type ccc "
      "--num_field 4 --num_vector_local %i --compute_method %s --all2all yes "
      "--num_proc_vector %i --num_proc_repl %i "
      "--num_proc_field %i --num_way %i --num_stage %i";

  int num_vector_local = 0;
  int num_proc_vector = 0;
  int num_proc_repl = 0;
  int gpu = 0;
  int num_stage = 1;

  for (gpu=0; gpu<=1; ++gpu) {
    for (num_vector_local=4; num_vector_local<=5; ++num_vector_local) {
      for (num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
        for (num_proc_repl=2; num_proc_repl<=5; ++num_proc_repl) {
          const int num_proc_field = gpu ? 2 : 1;
          const int num_way = 2;
          sprintf(options1, options_template_10,
                  num_vector_local*num_proc_vector,
                  "GPU", 1, 1,
                  1, num_way, 1);
          sprintf(options2, options_template_10, num_vector_local,
                  gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                  num_proc_field, num_way, num_stage);
          test_2runs(options1, options2);
        }
      }
    }
  }

  //----------
  //---num_repl, num_stage, 3-way
  //----------

  for (gpu=0; gpu<=1; ++gpu) {
    for (num_vector_local=6; num_vector_local<=18; num_vector_local+=12) {
      for (num_proc_vector=1; num_proc_vector<=4; ++num_proc_vector) {
        for (num_proc_repl=2; num_proc_repl<=5; ++num_proc_repl) {
          for (num_stage=1; num_stage<=6; num_stage+=4) {
            const int num_proc_field = gpu ? 2 : 1;
            const int num_way = 3;
            sprintf(options1, options_template_10,
                    num_vector_local*num_proc_vector,
                    "GPU", 1, 1,
                    1, num_way, 1);
            sprintf(options2, options_template_10, num_vector_local,
                    gpu ? "GPU" : "CPU", num_proc_vector, num_proc_repl,
                    num_proc_field, num_way, num_stage);
            test_2runs(options1, options2);
          }
        }
      }
    }
  }

  //----------
  //---ccc_param
  //----------

  char options_template_11a[] =
      "--metric_type ccc --verbosity %i "
      "--num_proc_vector 1 --num_field 30 --num_vector_local 40 "
      "--compute_method GPU";
  char options_template_11b[] =
      "--metric_type ccc --verbosity %i "
      "--num_proc_vector 1 --num_field 30 --num_vector_local 40 "
      "--compute_method GPU --ccc_param %.20e";

  sprintf(options1, options_template_11a, 1);
  sprintf(options2, options_template_11b, 1, ((double)2) / ((double)3));
  EXPECT_EQ(GM_BOOL_TRUE, compare_2runs(options1, options2));

  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  sprintf(options1, options_template_11a, 1);
  sprintf(options2, options_template_11b, 1, ((double)1) / ((double)2));
  const int result11 = compare_2runs(options1, options2);
  EXPECT_EQ(GM_BOOL_TRUE, proc_num==0 ? !result11 : GM_BOOL_TRUE);

  //----------
  //---num_phase, 2-way
  //----------

  int num_phase = 0;

  for (num_proc_vector=1; num_proc_vector<=8; ++num_proc_vector) {
    for (num_proc_repl=1; num_proc_repl<=8; ++num_proc_repl) {
      for (num_phase=2; num_phase<=8; ++num_phase) {
        if (!(num_phase <= 1 + num_proc_vector/2)) {
          continue;
        }
        char options_template[] =
          "--metric_type ccc "
          "--num_field 7 --num_vector 12 --compute_method GPU --all2all yes "
          "--num_proc_vector %i --num_proc_repl %i --num_phase %i "
          "--num_way 2";
        sprintf(options1, options_template, num_proc_vector, num_proc_repl,
                num_phase);
        sprintf(options2, options_template, 1, 1, 1);
        test_2runs(options1, options2);
      }
    }
  }

  //----------
  //---file output, 2-way
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_proc_vector 1 "
                    "--num_field 7 --num_vector 10 "
                    "--num_way 2 --metric_type ccc "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_vector 3 "
                    "--num_field 7 --num_vector 10 "
                    "--num_way 2 --metric_type ccc "
                    "--compute_method GPU --all2all yes "
                    "--verbosity 2 "
                    "--output_file_stub test_ccc_2way"));

  //----------
  //---file output, 3-way
  //----------

  EXPECT_EQ(
      GM_BOOL_TRUE,
      compare_2runs("--num_proc_vector 1 --num_proc_vector 1 "
                    "--num_field 7 --num_vector 18 "
                    "--num_way 3 --metric_type ccc "
                    "--compute_method REF --all2all yes",
                    "--num_proc_vector 1 --num_proc_vector 3 "
                    "--num_field 7 --num_vector 18 "
                    "--num_way 3 --metric_type ccc "
                    "--compute_method GPU --all2all yes "
                    "--verbosity 2 "
                    "--output_file_stub test_ccc_3way"));
}

/*===========================================================================*/

#if 1
TEST(DriverTest, czekanowski) {
  DriverTest_czekanowski_();
}

TEST(DriverTest, ccc2_simple) {
  DriverTest_ccc2_simple_();
}

TEST(DriverTest, ccc2_simple_sparse) {
  DriverTest_ccc2_simple_sparse_();
}

TEST(DriverTest, ccc3_simple) {
  DriverTest_ccc3_simple_();
}
#endif

#if 1
TEST(DriverTest, ccc3_simple_sparse) {
  DriverTest_ccc3_simple_sparse_();
}
#endif

#if 1
TEST(DriverTest, ccc) {
  DriverTest_ccc_();
}
#endif

/*===========================================================================*/

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

  mpi_code *= 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(&result, &result_g,
                           1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  MPI_Finalize();
  return result_g;
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
