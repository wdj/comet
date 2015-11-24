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

#include "driver_utils.h"

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
  GMEnv_create_from_args(&env, argc, (const char**)argv);

  /*---Parse remaining unprocessed arguments---*/

  int num_field = 0;
  int num_vector_local = 0;
  int verbosity = 0;
  finish_parsing(argc, (const char**)argv, &env, &num_field,
                 &num_vector_local, &verbosity);

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

    GMChecksum checksum = GMMetrics_checksum(&metrics, &env);

    if (Env_proc_num(&env) == 0 && verbosity > 0) {
      printf("metrics checksum ");
      int i = 0;
      for (i = 0; i < GM_CHECKSUM_SIZE; ++i ) {
        printf("%s%li", i==0 ? "" : "-", checksum.data[GM_CHECKSUM_SIZE-1-i]);
      }
      printf(" compute time %.6f\n", time_compute_metrics);
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
