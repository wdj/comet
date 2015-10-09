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
  "        manner of computing the result (0=CPU, 1=GPU)\n"
  "\n"
  "    ---verbosity <value>\n"
  "        verbosity level of output (0=none, 1=some (default) 2=more)\n"
  "\n"
  );
  /* clang-format on */
}

/*===========================================================================*/
/*---Set the entries of the vectors---*/

void input_vectors(Vectors* vectors, Env* env) {
  switch (data_type_id_from_metric_type(env->metric_type, env)) {
    case DATA_TYPE_ID_FLOAT: {
      int vector_local;
      for (vector_local = 0; vector_local < vectors->num_vector_local;
           ++vector_local) {
        int field;
        for (field = 0; field < vectors->num_field; ++field) {
          /*---compute element unique id---*/
          size_t index = field +
                         vectors->num_field * (vector_local +
                                               vectors->num_vector_local_max *
                                                   ((size_t)env->proc_num));
          /*---randomize---*/
          index = randomize(index);
          /*---set to number between 0 and 1---*/
          Float_t value = index / (Float_t)randomize_max();
          Vectors_float_set(vectors, field, vector_local, value, env);
        } /*---field---*/
      }   /*---vector_local---*/
    } break;
    case DATA_TYPE_ID_BIT:
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    default:
      Assert(Bool_false ? "Invalid data type." : 0);
  }
}

/*===========================================================================*/
/*---Output the result metrics values---*/

void output_metrics(Metrics* metrics, Env* env) {
  switch (data_type_id_from_metric_type(env->metric_type, env)) {
    case DATA_TYPE_ID_FLOAT: {
      size_t index;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        printf("proc: %i, entry (", env->proc_num);
        int coord_num = 0;
        for (coord_num = 0; coord_num < env->num_way; ++coord_num) {
          if (coord_num > 0) {
            printf(",");
          }
          /*---Present to the user as 1-based---*/
          printf("%i",
                 1 + Metrics_coord_from_index(metrics, index, coord_num, env));
        }
        printf("): value: %e\n", ((Float_t*)(metrics->data))[index]);
      } /*---for index---*/
    } break;
    case DATA_TYPE_ID_BIT:
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    default:
      Assert(Bool_false ? "Invalid data type." : 0);
  }
}

/*===========================================================================*/
/*---Parse remaining unprocessed arguments---*/

void finish_parsing(int argc,
                    char** argv,
                    Env* env,
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
      Insist(env, i < argc ? "Missing value for num_field." : 0);
      *num_field = atoi(argv[i]);
      Insist(env, *num_field >= 0 ? "Invalid setting for num_field." : 0);

    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
      ++i;
      Insist(env, i < argc ? "Missing value for num_vector_local." : 0);
      *num_vector_local = atoi(argv[i]);
      Insist(env, *num_vector_local >= 0
                      ? "Invalid setting for num_vector_local."
                      : 0);
    } else if (strcmp(argv[i], "--verbosity") == 0) {
      ++i;
      Insist(env, i < argc ? "Missing value for verbosity." : 0);
      *verbosity = atoi(argv[i]);
      Insist(env, *verbosity >= 0 ? "Invalid setting for verbosity." : 0);
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      ++i; /*---processed elsewhere by Env---*/
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i; /*---processed elsewhere by Env---*/
    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i; /*---processed elsewhere by Env---*/
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i; /*---processed elsewhere by Env---*/
    } else {
      if (env->proc_num == 0) {
        fprintf(stderr, "Invalid argument \"%s\".", argv[i]);
      }
      Insist(env, Bool_false ? "Error: argument not recognized." : 0);
    } /*---if/else---*/

  } /*---for i---*/

  Insist(env, *num_field != uninitialized ? "Error: num_field not set." : 0);
  Insist(env, *num_vector_local != uninitialized
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

  Env env = Env_null();
  Env_create_from_args(&env, argc, argv);

  /*---Parse remaining unprocessed arguments---*/

  int num_field = 0;
  int num_vector_local = 0;
  int verbosity = 0;
  finish_parsing(argc, argv, &env, &num_field, &num_vector_local, &verbosity);

  /*---Initialize vectors---*/

  Vectors vectors = Vectors_null();
  Vectors_create(&vectors, data_type_id_from_metric_type(env.metric_type, &env),
                 num_field, num_vector_local, &env);

  input_vectors(&vectors, &env);

  /*---Set up metrics container for results---*/

  Metrics metrics = Metrics_null();
  Metrics_create(&metrics, data_type_id_from_metric_type(env.metric_type, &env),
                 num_vector_local, &env);

  /*---Calculate metrics---*/

  /*---Run once first, discard timing---*/
  if (env.compute_method == COMPUTE_METHOD_GPU) {
    compute_metrics(&metrics, &vectors, &env);
  }

  double time_begin = Env_get_synced_time(&env);

  compute_metrics(&metrics, &vectors, &env);

  double time_end = Env_get_synced_time(&env);

  /*---Output run information---*/

  double time_compute_metrics = time_end - time_begin;

  double checksum = Metrics_checksum(&metrics, &env);

  if (env.proc_num == 0 && verbosity > 0) {
    printf("metrics checksum %.17e compute time %.6f\n", checksum,
           time_compute_metrics);
  }

  /*---Output results---*/

  if (verbosity > 1) {
    output_metrics(&metrics, &env);
  }

  /*---Finalize---*/

  Metrics_destroy(&metrics, &env);
  Vectors_destroy(&vectors, &env);
  Env_destroy(&env);
  MPI_Finalize();
  return 0;
}

/*---------------------------------------------------------------------------*/
