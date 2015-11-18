/*---------------------------------------------------------------------------*/
/*!
 * \file   genomics_metric.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code for genomics metric calculation.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

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
/*---Inform user of usage of command---*/

void usage() {
  /* clang-format off */
  printf(
  "genomics_metric: calculation of comparison metrics from genomics data\n"
  "\n"
  "Usage:\n"
  "\n"
  "    genomics_metric <option> ...\n"
  "\n"
  "Options:\n"
  "\n"
  "    --num_field <value>\n"
  "        (Required) the number of elements in each vector\n"
  "\n"
  "    --num_vector_local <value>\n"
  "        (Required) the number of vectors to be processed on each processor\n"
  "\n"
  "    --metric_type <value>\n"
  "        metric type to compute (sorenson=Sorenson,\n"
  "        czekanowski=Czekanowski (default), ccc=CCC)\n"
  "\n"
  "    --num_way <value>\n"
  "        dimension of metric to compute (2=2-way, 3=3-way\n"
  "\n"
  "    --all2all\n"
  "        whether to perform global all-to-all as rather than computing\n"
  "        on each processor separately (yes=yes, no=no (default))\n"
  "\n"
  "    --compute_method\n"
  "        manner of computing the result (CPU=cpu, GPU=gpu,\n"
  "        REF=reference method)\n"
  "\n"
  "    --num_proc\n"
  "        restrict run to this number of the available MPI processes\n"
  "        (used for testing purposes)\n"
  "\n"
  "    ---verbosity <value>\n"
  "        verbosity level of output (0=none, 1=some (default) 2=more)\n"
  "\n"
  );
  /* clang-format on */
}

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
          size_t index = field +
                         vectors->num_field * (vector_local +
                                               vectors->num_vector_local *
                                                  ((size_t)Env_proc_num(env)));
          /*---randomize---*/
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 1---*/
          GMFloat rand_value = index / (GMFloat)gm_randomize_max();
          /*---Create large integer in a specified range, store as float---*/
          GMFloat value = (int)((1 << 27) * rand_value);
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
                    char** argv,
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
/*---Main---*/

int main(int argc, char** argv) {
  /*---Initialize---*/

  MPI_Init(&argc, &argv);

  if (argc == 1) {
    usage();
    MPI_Finalize();
    return 0;
  }

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

    /*---Run once first, discard timing---*/
    if (Env_compute_method(&env) == GM_COMPUTE_METHOD_GPU) {
      gm_compute_metrics(&metrics, &vectors, &env);
    }

    double time_begin = GMEnv_get_synced_time(&env);

    gm_compute_metrics(&metrics, &vectors, &env);

    double time_end = GMEnv_get_synced_time(&env);

    /*---Output run information---*/

    double time_compute_metrics = time_end - time_begin;

    double checksum = GMMetrics_checksum(&metrics, &env);

    if (Env_proc_num(&env) == 0 && verbosity > 0) {
      printf("metrics checksum %.17e compute time %.6f\n", checksum,
             time_compute_metrics);
    }

    /*---Output results---*/

    if (verbosity > 1) {
      output_metrics(&metrics, &env);
    }

    /*---Finalize---*/

    GMMetrics_destroy(&metrics, &env);
    GMVectors_destroy(&vectors, &env);

  } /*---if (Env_is_proc_active(&env))---*/

  GMEnv_destroy(&env);
  MPI_Finalize();
  return 0;
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
