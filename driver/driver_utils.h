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
/*---Output the result metrics values---*/

void output_metrics(GMMetrics* metrics, GMEnv* env) {
  switch (Env_data_type(env)) {
    case GM_DATA_TYPE_FLOAT: {
      size_t index;
      for (index = 0; index < metrics->num_elts_local; ++index) {
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
        printf("): value: %.17e    [from proc %i]\n",
               GMMetrics_float_get_from_index(metrics, index, env),
               Env_proc_num(env));
      } /*---for index---*/
    } break;
    case GM_DATA_TYPE_BIT: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
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

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_driver_utils_h_---*/

/*---------------------------------------------------------------------------*/
