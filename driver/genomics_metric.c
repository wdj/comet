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
  "        metric type to compute (1=Sorenson, 2=Czekanowski (default),\n"
  "        3=CCC)\n"
  "\n"
  "    --num_way <value>\n"
  "        dimension of metric to compute (2=2-way, 3=3-way\n"
  "\n"
  "    --global_all2all\n"
  "        whether to perform global all-to-all as rather than computing\n"
  "        on each processor separately (1=yes, 0=no (default))\n"
  "\n"
  "    --compute_method\n"
  "        manner of computing the result (0=CPU, 1=GPU)\n"
  "\n"
  );
  /* clang-format on */
}

/*===========================================================================*/
/*---Set the vectors to test values---*/

void input_vectors(Vectors* vectors, Env* env) {

  double num_vector_elts_total = ((double)vectors->num_field *
                                  (double)vectors->num_vector_local_max *
                                  (double)env->num_proc);
  int vector_local;
  for (vector_local = 0; vector_local < vectors->num_vector_local;
       ++vector_local) {
    int field;
    for (field = 0; field < vectors->num_field; ++field) {
      /*---compute element unique id---*/
      size_t index = field + vectors->num_field * (
                     vector_local + vectors->num_vector_local_max * (
                     (size_t)env->proc_num));
      /*---randomize---*/
      index = randomize ( index );
      /*---set---*/
      Float_t value = index / num_vector_elts_total;
/*
printf("%e\n",value);
*/
      Vectors_set(vectors, field, vector_local, value, env);
    } /*---field---*/
  } /*---vector_local---*/
}

/*===========================================================================*/
/*---Output the result metrics values---*/

void output_metrics(Metrics* metrics, Env* env) {
  switch (data_type_id_from_metric_type(env->metric_type, env)) {
    case DATA_TYPE_ID_FLOAT:
      {
        size_t index;
        for ( index = 0; index < metrics->num_elts_local; ++index ) {
          printf("proc: %i  entry: ", env->proc_num);
          int i = 0;
          for ( i = 0; i < env->num_way; ++i ) {
            printf("%i ", Metrics_coord_from_index(metrics, index, i, env));
          }
          printf(" value: %e\n", ((Float_t*)(metrics->data))[index]);
        } /*---for index---*/
      }
      break;
    case DATA_TYPE_ID_BIT:
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    default:
      Assert(Bool_false ? "Invalid data type." : 0);
  }
}

/*===========================================================================*/
/*---Parse remaining unprocessed arguments---*/

void finish_parsing(int argc, char** argv, Env* env,
                    int* num_field,
                    int* num_vector_local) {

  const int uninitialized = -1;
  *num_field = uninitialized;
  *num_vector_local = uninitialized;

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
      Insist(env, *num_vector_local >= 0 ?
             "Invalid setting for num_vector_local." : 0);
    } else if (strcmp(argv[i], "--metric_type") == 0) {
        ++i;
    } else if (strcmp(argv[i], "--num_way") == 0) {
        ++i;
    } else if (strcmp(argv[i], "--global_all2all") == 0) {
        ++i;
    } else if (strcmp(argv[i], "--compute_method") == 0) {
        ++i;
    } else {
      if ( env->proc_num==0 ) {
        fprintf(stderr, "Invalid argument \"%s\".", argv[i]);
      }
      Insist(env, Bool_false ? "Error: argument not recognized." : 0);
    } /*---if/else---*/

  }   /*---for i---*/

  Insist(env, *num_field != uninitialized ? "Error: num_field not set." : 0);
  Insist(env, *num_vector_local != uninitialized ?
                                      "Error: num_vector_local not set." : 0);
}

/*===========================================================================*/
/*---Main---*/

int main(int argc, char** argv) {

  /*---Initialize---*/

  MPI_Init(&argc, &argv);

  if ( argc == 1 ) {
    usage();
    MPI_Finalize();
    return 0;
  }

  Env env = Env_null();
  Env_create_from_args(&env, argc, argv);

  /*---Parse remaining unprocessed arguments---*/

  int num_field = 0;
  int num_vector_local = 0;
  finish_parsing(argc, argv, &env, &num_field, &num_vector_local);

  /*---Initialize vectors---*/

  Vectors vectors = Vectors_null();
  Vectors_create(&vectors, data_type_id_from_metric_type(env.metric_type, &env),
                 num_field, num_vector_local, &env);

  input_vectors(&vectors, &env);

  /*---Set up metrics---*/

  Metrics metrics = Metrics_null();
  Metrics_create(&metrics, data_type_id_from_metric_type(env.metric_type, &env),
                 num_vector_local, &env);

  /*---Calculate metrics---*/

  double time_begin = Env_get_synced_time(&env);

  compute_metrics(&metrics, &vectors, &env);

  double time_end = Env_get_synced_time(&env);

  double time_compute_metrics = time_end - time_begin;

  double checksum = Metrics_checksum ( &metrics, &env );

  if ( env.proc_num==0 ) {
    printf("metrics checksum %.17e compute time %.6f\n",
           checksum, time_compute_metrics);
  }

  /*---Output results---*/

  output_metrics(&metrics, &env);

  /*---Finalize---*/

  Vectors_destroy(&vectors, &env);
  Env_destroy(&env);
  MPI_Finalize();
  return 0;
}

/*---------------------------------------------------------------------------*/
