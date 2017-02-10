/*---------------------------------------------------------------------------*/
/*!
 * \file   genomics_metric.cc
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

#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>

#include "mpi.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics.hh"

#include "driver_utils.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/* Stack tracing code */

#if 0
// http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes */

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}
#endif

// http://stackoverflow.com/questions/4636456/how-to-get-a-stack-trace-for-c-using-gcc-with-line-number-information
// http://stackoverflow.com/questions/3151779/how-its-better-to-invoke-gdb-from-program-to-print-its-stacktrace

void bt_sighandler(int sig, struct sigcontext ctx) {

  int proc_num = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

  char pid_buf[30];
  sprintf(pid_buf, "%d", getpid());
  char name_buf[512];
  name_buf[readlink("/proc/self/exe", name_buf, 511)] = 0;
  int child_pid = fork();
  if (!child_pid) {           
    dup2(2,1); // redirect output to stderr
    fprintf(stdout,"MPI rank %i: stack trace for %s pid=%s\n",
            proc_num, name_buf, pid_buf);

    char cmd[512];
    sprintf(cmd, "gdb --batch -n -ex thread -ex bt %s %s 2>&1"
                 " | grep '^#' | sed -e 's/^/MPI rank %i: /'",
            name_buf, pid_buf, proc_num);
    system(cmd);

    abort(); /* If gdb failed to start */
  } else {
    waitpid(child_pid,NULL,0);
  }

  exit(0);
}

/*---------------------------------------------------------------------------*/

void install_handler() {

  //signal(SIGSEGV, handler);

  struct sigaction sa;

  sa.sa_handler = (__sighandler_t)bt_sighandler;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART;

  // http://www.comptechdoc.org/os/linux/programming/linux_pgsignals.html

  sigaction(SIGSEGV, &sa, NULL);
  sigaction(SIGFPE, &sa, NULL);
  sigaction(SIGINT, &sa, NULL);
  sigaction(SIGILL, &sa, NULL);
  sigaction(SIGUSR1, &sa, NULL);
  sigaction(SIGINT, &sa, NULL);
}

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
  "    --num_field_local <value>\n"
  "        (Required) the number of elements in each vector on each processor\n"
  "\n"
  "    --num_vector <value>\n"
  "        (Required) the number of vectors to be processed"
  "\n"
  "    --num_vector_local <value>\n"
  "        (Required) the number of vectors to be processed on each processor\n"
  "\n"
  "    --metric_type <value>\n"
  "        metric type to compute (sorenson=Sorenson,\n"
  "        czekanowski=Czekanowski (default), ccc=CCC)\n"
  "\n"
  "    --num_way <value>\n"
  "        dimension of metric to compute (2=2-way (default), 3=3-way)\n"
  "\n"
  "    --all2all\n"
  "        whether to perform global all-to-all rather than computing\n"
  "        on each processor separately (yes=yes, no=no (default))\n"
  "\n"
  "    --compute_method\n"
  "        manner of computing the result (CPU=cpu, GPU=gpu (default),\n"
  "        REF=reference method)\n"
  "\n"
  "    --num_proc_vector\n"
  "        blocking factor to denote number of blocks used to decompose\n"
  "        the total number of vectors across processors \n"
  "        (default is the total number of procs requested)\n"
  "\n"
  "    --num_proc_field\n"
  "        blocking factor to denote number of blocks used to decompose\n"
  "        each vector across processors (default is 1)\n"
  "\n"
  "    --num_proc_repl\n"
  "        processor replication factor.  For each block along the vector\n"
  "        and field axes, this number of processors is applied to\n"
  "        computations for the block (default is 1)\n"
  "\n"
  "    --num_stage <value>\n"
  "        the number of stages the computation is divided into\n"
  "        (default is 1) (available for 3-way case only)\n"
  "\n"
  "    --stage_min <value>\n"
  "        the lowest stage number of the sequence of stages to be computed\n"
  "        for this run (default is 1)\n"
  "\n"
  "    --stage_max <value>\n"
  "        the highest stage number of the sequence of stages to be computed\n"
  "        for this run (default is num_stage)\n"
  "\n"
  "    --verbosity <value>\n"
  "      verbosity level of output (0=none, 1=some (default) 2=more)\n"
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

  //install_handler();

  /*---If using GPU---*/

  GMEnv env = GMEnv_null();
  GMEnv_create_from_args(&env, argc, (char**)argv, NULL);

  if (GMEnv_compute_method(&env) == GM_COMPUTE_METHOD_GPU) {
    /*---Perform preliminary run on GPU since sometimes first use is slower---*/

//FIX this to run on all nodes.
    const char* options1 =
        "--num_field 1 --num_vector_local 2 "
        "--compute_method GPU --verbosity 0";
    size_t len1 = strlen(options1);
    char* argstring1 = (char*)malloc((len1 + 1) * sizeof(char));
    GMAssertAlways(argstring1 != NULL);
    char* argv1[len1 + 1];
    int argc1 = 0;
    strcpy(argstring1, options1);
    gm_create_args(argstring1, &argc1, argv1);

    perform_run(argc1, (char**)argv1, NULL);

    free(argstring1);
  }

  GMEnv_destroy(&env);

  /*---Perform actual run---*/

  perform_run(argc, (char**)argv, NULL);

  MPI_Finalize();
  return 0;
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
