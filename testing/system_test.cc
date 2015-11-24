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

/*===========================================================================*/
/*---Set the entries of the vectors---*/

void input_vectors(GMVectors* vectors, GMEnv* env) {
  switch (Env_data_type(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        int field;
        for (field = 0; field < vectors->num_field; ++field) {
          /*---compute element unique id---*/
          const size_t uid = field +
                         vectors->num_field * (vector_local +
                                               vectors->num_vector_local *
                                                  ((size_t)Env_proc_num(env)));
          /*---Generate large random number---*/
          size_t rand1 = uid;
          rand1 = gm_randomize(rand1);
          rand1 = gm_randomize(rand1);
          size_t rand2 = uid;
          rand2 = gm_randomize(rand2);
          rand2 = gm_randomize(rand2);
          rand2 = gm_randomize(rand2);
          size_t rand_value = rand1 + gm_randomize_max() * rand2;
          /*---Reduce so that after summing num_field times the integer
               still fully fits in double precision fraction part---*/
          rand_value >>= (64-52) + gm_log2(vectors->num_field);
          /*---Store as floating point value---*/
          GMFloat value = rand_value;
          GMVectors_float_set(vectors, field, vector_local, value, env);
        } /*---field---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BIT: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        int field;
        for (field = 0; field < vectors->num_field; ++field) {
          /*---compute element unique id---*/
          size_t index = field +
                         vectors->num_field * (vector_local +
                                               vectors->num_vector_local *
                                                  ((size_t)Env_proc_num(env)));
          /*---randomize---*/
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 1---*/
          GMFloat rand_value = index / (GMFloat)gm_randomize_max();
          /*---Create single bit value---*/
          _Bool value = rand_value < .5 ? GM_BOOL_FALSE : GM_BOOL_TRUE;
          GMVectors_bit_set(vectors, field, vector_local, value, env);
        } /*---field---*/
      }   /*---vector_local---*/

    } break;
    /*--------------------*/
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  }
}

/*===========================================================================*/
/*---Parse remaining unprocessed arguments---*/

void finish_parsing(int argc,
                    const char** argv,
                    GMEnv* env,
                    int* num_field,
                    int* num_vector_local,
                    int* verbosity) {
  const int uninitialized = -1;
  *num_field = uninitialized;
  *num_vector_local = uninitialized;
  *verbosity = 1;

  int i = 0;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--num_field") == 0) {
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_field." : 0);
      *num_field = atoi(argv[i]);
      GMInsist(env, *num_field >= 0 ? "Invalid setting for num_field." : 0);

    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_vector_local." : 0);
      *num_vector_local = atoi(argv[i]);
      GMInsist(env, *num_vector_local >= 0
                        ? "Invalid setting for num_vector_local."
                        : 0);
    } else if (strcmp(argv[i], "--verbosity") == 0) {
      ++i;
      GMInsist(env, i < argc ? "Missing value for verbosity." : 0);
      *verbosity = atoi(argv[i]);
      GMInsist(env, *verbosity >= 0 ? "Invalid setting for verbosity." : 0);
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else {
      if (Env_proc_num(env) == 0) {
        fprintf(stderr, "Invalid argument \"%s\".", argv[i]);
      }
      GMInsist(env, GM_BOOL_FALSE ? "Error: argument not recognized." : 0);
    } /*---if/else---*/

  } /*---for i---*/

  GMInsist(env, *num_field != uninitialized ? "Error: num_field not set." : 0);
  GMInsist(env, *num_vector_local != uninitialized
                    ? "Error: num_vector_local not set."
                    : 0);
}

/*===========================================================================*/
/*---Parse string to construct arguments---*/

void create_args(char* argstring, int* argc, const char **argv) {

  size_t len = strlen(argstring);

  argv[0] = &argstring[0];
  *argc = 1;
  _Bool is_delim_prev = GM_BOOL_TRUE;
  int i = 0;
  for (i=0; i<(int)len; ++i) {
    const _Bool is_delim = argstring[i] == ' ' || argstring[i] == '\t';
    if (is_delim) {
      argstring[i] = 0;
    }
    if (is_delim_prev && ! is_delim) {
      argv[*argc] = &(argstring[i]);
      (*argc)++;
    }
    is_delim_prev = is_delim;
  }
}

/*===========================================================================*/

GMChecksum perform_run(const char* options) {

  GMChecksum checksum;

  /*---Convert options string to args---*/

  size_t len = strlen(options);
  char* argstring = (char*)malloc((len+1)*sizeof(char));
  strcpy(argstring, options);
  const char* argv[len+1];
  int argc = 0;

  create_args(argstring, &argc, argv);

  /*---Initialize environment---*/

  GMEnv env = GMEnv_null();
  GMEnv_create_from_args(&env, argc, argv);

  /*---Parse remaining unprocessed arguments---*/

  int num_field = 0;
  int num_vector_local = 0;
  int verbosity = 0;
  finish_parsing(argc, argv, &env, &num_field, &num_vector_local, &verbosity);

  if (Env_is_proc_active(&env)) {

    /*---Initialize vectors---*/

    GMVectors vectors = GMVectors_null();
    GMVectors_create(&vectors, Env_data_type(&env),
                     num_field, num_vector_local, &env);

    input_vectors(&vectors, &env);

    /*---Set up metrics container for results---*/

    GMMetrics metrics = GMMetrics_null();
    GMMetrics_create(&metrics, Env_data_type(&env),
                     num_vector_local, &env);

    /*---Calculate metrics---*/

    /*---Run once first, discard timing: first gpu run is sometimes slow---*/
    if (Env_compute_method(&env) == GM_COMPUTE_METHOD_GPU) {
      gm_compute_metrics(&metrics, &vectors, &env);
    }

    double time_begin = GMEnv_get_synced_time(&env);

    gm_compute_metrics(&metrics, &vectors, &env);

    double time_end = GMEnv_get_synced_time(&env);

    /*---Output run information---*/

    double time_compute_metrics = time_end - time_begin;

    checksum = GMMetrics_checksum(&metrics, &env);

    if (Env_proc_num(&env) == 0) {
      printf("metrics checksum ");
      int i = 0;
      for (i = 0; i < GM_CHECKSUM_SIZE; ++i ) {
        printf("%s%li", i==0 ? "" : "-", checksum.data[GM_CHECKSUM_SIZE-1-i]);
      }
      printf(" compute time %.6f\n", time_compute_metrics);
    }

    /*---Finalize---*/

    GMMetrics_destroy(&metrics, &env);
    GMVectors_destroy(&vectors, &env);

  } /*---if (Env_is_proc_active(&env))---*/

  GMEnv_destroy(&env);

  free(argstring);

  return checksum;
}

/*===========================================================================*/

_Bool compare_runs(const char* options1, const char* options2) {

  GMChecksum checksum1 = perform_run(options1);
  GMChecksum checksum2 = perform_run(options2);

  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);
  if (proc_num != 0) {
    int i = 0;
    for (i = 0; i < GM_CHECKSUM_SIZE; ++i ) {
      checksum1.data[i] = 0;
      checksum2.data[i] = 0;
    }
  }

  return gm_are_checksums_equal(checksum1, checksum2);
}

/*===========================================================================*/

TEST(SystemTest,One) {

  //int proc_num = 0;
  //MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 1 --num_vector_local 2 --compute_method CPU",
    "--num_proc 1 --num_field 1 --num_vector_local 2 --compute_method GPU"));

  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU",
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method GPU"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU",
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU"
      " --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU",
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method GPU"
      " --all2all yes"));

  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU",
    "--num_proc 2 --num_field 100 --num_vector_local 24 --compute_method CPU"
      " --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU",
    "--num_proc 2 --num_field 100 --num_vector_local 24 --compute_method GPU"
      " --all2all yes"));

  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 1 --num_vector_local 3 --compute_method CPU "
      "--num_way 3",
    "--num_proc 1 --num_field 1 --num_vector_local 3 --compute_method GPU "
      "--num_way 3"));

  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU "
      "--num_way 3",
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method GPU "
      "--num_way 3"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU "
      "--num_way 3",
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU "
      "--num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method GPU "
      "--num_way 3",
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method GPU "
      "--num_way 3 --all2all yes"));

  //----------

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method CPU "
      "--num_way 3 --all2all yes",
    "--num_proc 4 --num_field 100 --num_vector_local 12 --compute_method CPU "
      "--num_way 3 --all2all yes"));

  EXPECT_EQ(GM_BOOL_TRUE, compare_runs(
    "--num_proc 1 --num_field 100 --num_vector_local 48 --compute_method GPU "
      "--num_way 3 --all2all yes",
    "--num_proc 4 --num_field 100 --num_vector_local 12 --compute_method GPU "
      "--num_way 3 --all2all yes"));

  //----------
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
