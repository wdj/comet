/*---------------------------------------------------------------------------*/
/*!
 * \file   driver_utils.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_driver_utils_hh_
#define _gm_driver_utils_hh_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>

#include "mpi.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Utility to parse a string to construct arguments---*/

static void create_args(char* argstring, int* argc, char** argv) {
  size_t len = strlen(argstring);

  argv[0] = &argstring[0];
  *argc = 1;
  _Bool is_delim_prev = GM_BOOL_TRUE;
  int i = 0;
  for (i = 0; i < (int)len; ++i) {
    const _Bool is_delim = argstring[i] == ' ' || argstring[i] == '\t';
    if (is_delim) {
      argstring[i] = 0;
    }
    if (is_delim_prev && !is_delim) {
      argv[*argc] = &(argstring[i]);
      (*argc)++;
    }
    is_delim_prev = is_delim;
  }
}

/*===========================================================================*/
/*---Set the entries of the vectors---*/

static void input_vectors(GMVectors* vectors, int verbosity, GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_BITS1: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      //---(design is not complete)

      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        size_t vector =
            vector_local +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        int field_local;
        for (field_local = 0; field_local < vectors->num_field_local;
             ++field_local) {
          size_t field =
              field_local +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
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
        size_t vector =
            vector_local +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        int fl = 0;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
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
          const size_t randomize_max = gm_randomize_max();
          size_t rand_value = rand1 + randomize_max * rand2;
          /*---Reduce so that after summing num_field times the integer
               still exactly representable by floating point type---*/
          const size_t rand_max = randomize_max * randomize_max;
          GMAssertAlways(FLT_RADIX == 2);
          const int mant_dig = sizeof(GMFloat) == 8 ? DBL_MANT_DIG :
                                                      FLT_MANT_DIG;
          const int shift_amount = gm_log2(rand_max * vectors->num_field) -
                                   mant_dig;
          rand_value >>= shift_amount > 0 ? shift_amount : 0;
#if 0
          if (sizeof(GMFloat) == 8) {
            const int shift_amount = (64 - 52) + gm_log2(vectors->num_field);
            rand_value >>= shift_amount;
          } else {
            const int shift_amount = gm_log2(rand_value*rand_value) - 22
              + gm_log2(vectors->num_field);


            const int shift_amount = gm_log2(rand_max * vectors->num_field) -
                                     mant_dig;
            rand_value >>= shift_amount > 0 ? shift_amount : 0;
          }
#endif
          /*---Store---*/
          GMFloat float_value = (GMFloat)rand_value;
          GMAssertAlways((size_t)float_value == rand_value);
          GMAssertAlways(float_value * vectors->num_field <=
                         ((size_t)1)<<mant_dig);
          GMVectors_float_set(vectors, fl, vector_local, float_value, env);
          /*---Print---*/
          if (verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                   GMEnv_proc_num_vector_i(env), vector_local,
                   GMEnv_proc_num_field(env), fl, float_value);
          }
        } /*---field---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        size_t vector =
            vector_local +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        int fl;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          /*---Compute element unique id---*/
          const size_t uid = field + vectors->num_field * vector;
          size_t index = uid;
          /*---Randomize---*/
          index = gm_randomize(index);
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 3---*/
          const float float_rand_value = index / (float)gm_randomize_max();
          /*---Create 2-bit value - make extra sure less than 4---*/
          GMBits2 value = (int)((4. - 1e-5) * float_rand_value);
          /*---Store---*/
          GMVectors_bits2_set(vectors, fl, vector_local, value, env);
          /*---Print---*/
          if (verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %.1i%.1i\n",
                   GMEnv_proc_num_vector_i(env), vector_local,
                   GMEnv_proc_num_field(env), fl, value / 2, value % 2);
          }
        } /*---field---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    default:
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Output the result metrics values---*/

static void output_metrics(GMMetrics* metrics, GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(env != NULL);

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  /*---Due to redundancy, only results from some processors are needed---*/
  if (GMEnv_proc_num_field(env) != 0) {
    return;
  }

  switch (GMEnv_data_type_metrics(env)) {
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
        for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
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
        printf("    [from proc %i]\n", GMEnv_proc_num(env));
      } /*---for index---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
      size_t index;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        /*---Coordinate---*/
        printf("element (");
        int coord_num = 0;
        for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
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
        for (i = 0; i < 2; ++i) {
          int j;
          for (j = 0; j < 2; ++j) {
            printf(" %.17e",
                   GMMetrics_ccc_get_from_index_2(metrics, index, i, j, env));
          }
        }
        printf("    [from proc %i]\n", GMEnv_proc_num(env));
      } /*---for index---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
      size_t index;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        /*---Coordinate---*/
        printf("element (");
        int coord_num = 0;
        for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
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
        for (i = 0; i < 2; ++i) {
          int j;
          for (j = 0; j < 2; ++j) {
            int k;
            for (k = 0; k < 2; ++k) {
              printf(" %.17e", GMMetrics_ccc_get_from_index_3(metrics, index, i,
                                                              j, k, env));
            }
          }
        }
        printf("    [from proc %i]\n", GMEnv_proc_num(env));
      } /*---for index---*/
    } break;
    /*--------------------*/
    default:
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Parse remaining unprocessed arguments---*/

static void finish_parsing(int argc,
                           char** argv,
                           GMEnv* env,
                           int* num_field,
                           int* num_vector_local,
                           int* verbosity,
                           int* num_stage) {
  const int uninitialized = -1;
  *num_field = uninitialized;
  *num_vector_local = uninitialized;
  *verbosity = 1;

  int i = 0;
  for (i = 1; i < argc; ++i) {
    /*----------*/
    if (strcmp(argv[i], "--num_field") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_field." : 0);
      //FIX: safe atoi
      *num_field = atoi(argv[i]);
      GMInsist(env, *num_field >= 0 ? "Invalid setting for num_field." : 0);

    /*----------*/
    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_vector_local." : 0);
      //FIX: safe atoi
      *num_vector_local = atoi(argv[i]);
      GMInsist(env, *num_vector_local >= 0
                        ? "Invalid setting for num_vector_local."
                        : 0);
    /*----------*/
    } else if (strcmp(argv[i], "--verbosity") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for verbosity." : 0);
      //FIX: safe atoi
      *verbosity = atoi(argv[i]);
      GMInsist(env, *verbosity >= 0 ? "Invalid setting for verbosity." : 0);
    /*----------*/
    } else if (strcmp(argv[i], "--num_stage") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_stage." : 0);
      //FIX: safe atoi
      *num_stage = atoi(argv[i]);
      GMInsist(env, *num_stage >= 1 ? "Invalid setting for num_stage." : 0);
    /*----------*/
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
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else {
    /*----------*/
      if (GMEnv_proc_num(env) == 0) {
        fprintf(stderr, "Invalid argument \"%s\".", argv[i]);
      }
      GMInsist(env, GM_BOOL_FALSE ? "Error: argument not recognized." : 0);
    /*----------*/
    } /*---if/else---*/

  } /*---for i---*/

  GMInsist(env, *num_field != uninitialized ? "Error: num_field not set." : 0);
  GMInsist(env, *num_vector_local != uninitialized
                    ? "Error: num_vector_local not set."
                    : 0);
}

/*===========================================================================*/
/*---Perform a single metrics computation run---*/

static GMChecksum perform_run(int argc, char** argv,
                              char const * const description) {

  /*---Initialize environment---*/

  GMEnv env = GMEnv_null();
  GMEnv_create_from_args(&env, argc, argv, description);

  /*---Parse remaining unprocessed arguments---*/

  int num_field = 0;
  int num_vector_local = 0;
  int verbosity = 0;
  int num_stage = 1;
  finish_parsing(argc, argv, &env, &num_field, &num_vector_local, &verbosity,
                 &num_stage);
  env.num_stage = num_stage;

  /*---Initialize vectors---*/

  GMVectors vectors = GMVectors_null();
  GMVectors_create(&vectors, GMEnv_data_type_vectors(&env), num_field,
                   num_vector_local, &env);

  input_vectors(&vectors, verbosity, &env);

  GMChecksum checksum = GMChecksum_null();

  /*---Loop over stages---*/

  int stage_num = 0;
  for (stage_num=0; stage_num<env.num_stage; ++stage_num) {

    env.stage_num = stage_num;

    /*---Set up metrics container for results---*/

    GMMetrics metrics = GMMetrics_null();
    GMMetrics_create(&metrics, GMEnv_data_type_metrics(&env), num_field,
                     num_vector_local, &env);

    /*---Calculate metrics---*/

    gm_compute_metrics(&metrics, &vectors, &env);

    /*---Output results---*/

    if (verbosity > 1) {
      output_metrics(&metrics, &env);
    }

    /*---Compute checksum---*/

    GMMetrics_checksum(&metrics, &checksum, &env);

    GMMetrics_destroy(&metrics, &env);

  }

  /*---Output run information---*/

  if (GMEnv_is_proc_active(&env) && GMEnv_proc_num(&env) == 0 && verbosity > 0) {
    printf("metrics checksum ");
    int i = 0;
    for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
      printf("%s%li", i == 0 ? "" : "-",
             checksum.data[GM_CHECKSUM_SIZE - 1 - i]);
    }
    if (checksum.is_overflowed) {
      printf("-OVFL");
      printf("-%e", checksum.value_max);
    }
    printf(" time %.6f", env.time);
    printf(" ops %e", env.ops);
    if (env.time > 0) {
      printf(" rate %e", env.ops / env.time);
    }
    printf("\n");
  }

  /*---Finalize---*/

  GMVectors_destroy(&vectors, &env);

  GMEnv_destroy(&env);

  return checksum;
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_gm_driver_utils_hh_---*/

/*---------------------------------------------------------------------------*/
