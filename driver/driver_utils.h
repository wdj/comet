/*---------------------------------------------------------------------------*/
/*!
 * \file   driver_utils.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _driver_utils_h_
#define _driver_utils_h_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "mpi.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Utility to parse a string to construct arguments---*/

static void create_args(char* argstring, int* argc, const char **argv) {

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
/*---Set the entries of the vectors---*/

static void input_vectors(GMVectors* vectors, int verbosity, GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (!Env_is_proc_active(env)) {
    return;
  }

  switch (Env_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_BITS1: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      //---(design is not complete)

      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        size_t vector = vector_local + vectors->num_vector_local *
                                              (size_t)Env_proc_num_vector(env);
        int field_local;
        for (field_local = 0; field_local < vectors->num_field_local;
             ++field_local) {
          size_t field = field_local + vectors->num_field_local *
                                              (size_t)Env_proc_num_field(env);
          /*---compute element unique id---*/
          const size_t uid = field + vectors->num_field * vector;
          size_t index = uid;
          /*---randomize---*/
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 1---*/
          GMFloat rand_value = index / (GMFloat)gm_randomize_max();
          /*---Create single bit value---*/
          _Bool value = rand_value < .5 ? GM_BOOL_FALSE : GM_BOOL_TRUE;
          /*---Store---*/
          GMVectors_bit_set(vectors, field_local, vector_local, value, env);
          /*---Print---*/
          if (verbosity > 2) {
          }
        } /*---field---*/
      }   /*---vector_local---*/

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        size_t vector = vector_local + vectors->num_vector_local *
                                              (size_t)Env_proc_num_vector(env);
        int field_local = 0;
        for (field_local = 0; field_local < vectors->num_field_local;
             ++field_local) {
          size_t field = field_local + vectors->num_field_local *
                                              (size_t)Env_proc_num_field(env);
          /*---Compute element unique id---*/
          const size_t uid = field + vectors->num_field * vector;
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
          /*---Store---*/
          GMFloat value = rand_value;
          GMVectors_float_set(vectors, field_local, vector_local, value, env);
          /*---Print---*/
          if (verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                   Env_proc_num_vector(env), vector_local,
                   Env_proc_num_field(env), field_local, 
                   value);
          }
        } /*---field---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        size_t vector = vector_local + vectors->num_vector_local *
                                              (size_t)Env_proc_num_vector(env);
        int field_local;
        for (field_local = 0; field_local < vectors->num_field_local;
             ++field_local) {
          size_t field = field_local + vectors->num_field_local *
                                              (size_t)Env_proc_num_field(env);
          /*---Compute element unique id---*/
          const size_t uid = field + vectors->num_field * vector;
          size_t index = uid;
          /*---Randomize---*/
          index = gm_randomize(index);
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 3---*/
          const float float_rand_value = index / (float)gm_randomize_max();
          /*---Create 2-bit value - make extra sure less than 4---*/
          GMBits2 value = (int)( (4.-1e-5) * float_rand_value);
          /*---Store---*/
          GMVectors_bits2_set(vectors, field_local, vector_local, value, env);
          /*---Print---*/
          if (verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %.1i%.1i\n",
                   Env_proc_num_vector(env), vector_local,
                   Env_proc_num_field(env), field_local, 
                   value/2, value%2);
          }
        } /*---field---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Output the result metrics values---*/

static void output_metrics(GMMetrics* metrics, GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);

  if (!Env_is_proc_active(env)) {
    return;
  }

  /*---Due to redundancy, only results from some processors are needed---*/
  if (Env_proc_num_field(env) != 0) {
    return;
  }

  switch (Env_data_type_metrics(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_BITS1: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
      size_t index;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        /*---Coordinate---*/
        printf("element (");
        int coord_num = 0;
        for (coord_num = 0; coord_num < Env_num_way(env); ++coord_num) {
          if (coord_num > 0) {
            printf(",");
          }
          /*---Present to the user as 1-based---*/
          printf("%i", 1 + GMMetrics_coord_global_from_index(metrics, index,
                                                             coord_num, env));
        }
        /*---Value---*/
        printf("): value:");
        printf(" %.17e",
               GMMetrics_czekanowski_get_from_index(metrics, index, env));
        printf("    [from proc %i]\n", Env_proc_num(env));
      } /*---for index---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
      size_t index;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        /*---Coordinate---*/
        printf("element (");
        int coord_num = 0;
        for (coord_num = 0; coord_num < Env_num_way(env); ++coord_num) {
          if (coord_num > 0) {
            printf(",");
          }
          /*---Present to the user as 1-based---*/
          printf("%i", 1 + GMMetrics_coord_global_from_index(metrics, index,
                                                             coord_num, env));
        }
        /*---Value---*/
        printf("): values:");
        int i;
        for (i=0; i<2; ++i) {
          int j;
          for (j=0; j<2; ++j) {
            printf(" %.17e", GMMetrics_ccc_get_from_index_2(metrics, index,
                                                            i, j, env));
          }
        }
        printf("    [from proc %i]\n", Env_proc_num(env));
      } /*---for index---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
      size_t index;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        /*---Coordinate---*/
        printf("element (");
        int coord_num = 0;
        for (coord_num = 0; coord_num < Env_num_way(env); ++coord_num) {
          if (coord_num > 0) {
            printf(",");
          }
          /*---Present to the user as 1-based---*/
          printf("%i", 1 + GMMetrics_coord_global_from_index(metrics, index,
                                                             coord_num, env));
        }
        /*---Value---*/
        printf("): values:");
        int i;
        for (i=0; i<2; ++i) {
          int j;
          for (j=0; j<2; ++j) {
            int k;
            for (k=0; k<2; ++k) {
              printf(" %.17e", GMMetrics_ccc_get_from_index_3(metrics, index,
                                                              i, j, k, env));
            }
          }
        }
        printf("    [from proc %i]\n", Env_proc_num(env));
      } /*---for index---*/
    } break;
    /*--------------------*/
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Parse remaining unprocessed arguments---*/

static void finish_parsing(int argc,
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
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
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
/*---Perform a single metrics computation run---*/

static GMChecksum perform_run(int argc, const char** argv) {

  GMChecksum checksum;

  /*---Initialize environment---*/

  GMEnv env = GMEnv_null();
  GMEnv_create_from_args(&env, argc, argv);

  /*---Parse remaining unprocessed arguments---*/

  int num_field = 0;
  int num_vector_local = 0;
  int verbosity = 0;
  finish_parsing(argc, argv, &env, &num_field, &num_vector_local, &verbosity);

  /*---Initialize vectors---*/

  GMVectors vectors = GMVectors_null();
  GMVectors_create(&vectors, Env_data_type_vectors(&env),
                   num_field, num_vector_local, &env);

  input_vectors(&vectors, verbosity, &env);

  /*---Set up metrics container for results---*/

  GMMetrics metrics = GMMetrics_null();
  GMMetrics_create(&metrics, Env_data_type_metrics(&env),
                   num_field, num_vector_local, &env);

  /*---Calculate metrics---*/

  gm_compute_metrics(&metrics, &vectors, &env);

  /*---Output results---*/

  if (verbosity > 1) {
    output_metrics(&metrics, &env);
  }

  /*---Output run information---*/

  checksum = GMMetrics_checksum(&metrics, &env);

  if (Env_is_proc_active(&env) && Env_proc_num(&env) == 0 && verbosity > 0) {
    printf("metrics checksum ");
    int i = 0;
    for (i = 0; i < GM_CHECKSUM_SIZE; ++i ) {
      printf("%s%li", i==0 ? "" : "-", checksum.data[GM_CHECKSUM_SIZE-1-i]);
    }
    printf(" time %.6f", env.time);
    printf(" ops %e", env.ops);
    if (env.time > 0) {
      printf(" rate %e", env.ops / env.time);
    }
    printf("\n");
  }

  /*---Finalize---*/

  GMMetrics_destroy(&metrics, &env);
  GMVectors_destroy(&vectors, &env);

  GMEnv_destroy(&env);

  return checksum;
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_driver_utils_h_---*/

/*---------------------------------------------------------------------------*/
