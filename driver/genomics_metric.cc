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

  /*---Perform preliminary run on GPU since sometimes first use is slower---*/

  const char* options1 = "--num_proc_vector 1 --num_field 1 --num_vector_local 2"
                                      " --compute_method GPU --verbosity 0";
  size_t len1 = strlen(options1);
  char* argstring1 = (char*)malloc((len1+1)*sizeof(char));
  const char* argv1[len1+1];
  int argc1 = 0;
  strcpy(argstring1, options1);
  create_args(argstring1, &argc1, argv1);

  perform_run(argc1, (const char**)argv1);

  free (argstring1);

  /*---Perform actual run---*/

  perform_run(argc, (const char**)argv);

  MPI_Finalize();
  return 0;
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
