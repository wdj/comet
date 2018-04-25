//-----------------------------------------------------------------------------
/*!
 * \file   genomics_metric.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code for genomics metric calculation.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdio.h"
#include "stdlib.h"
#include "stddef.h"
#include "string.h"

#include "execinfo.h"
#include "signal.h"
#include "unistd.h"
#include "sys/wait.h"

#include "mpi.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics.hh"

#include "driver.hh"

//=============================================================================
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
  if (! child_pid) {           
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

//-----------------------------------------------------------------------------

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

//=============================================================================
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
  "        metric type to compute (czekanowski=Czekanowski (default),\n"
  "        ccc=CCC)\n"
  "\n"
  "    --ccc_multiplier <value>\n"
  "        fixed front multiplier value used to calculate the CCC metric\n"
  "        (default floating point value is 9/2).\n"
  "\n"
  "    --ccc_param <value>\n"
  "        fixed coefficient value used to calculate the CCC metric\n"
  "        (default floating point value is 2/3).\n"
  "\n"
  "    --sparse <value>\n"
  "        for the CCC metric, interpret vector entries of binary \"10\"\n"
  "        as empty or incomplete data (yes=yes, no=no (default))\n"
  "\n"
  "    --num_way <value>\n"
  "        dimension of metric to compute (2=2-way (default), 3=3-way)\n"
  "\n"
  "    --all2all <value>\n"
  "        whether to perform global all-to-all rather than computing\n"
  "        on each processor separately (yes=yes, no=no (default))\n"
  "\n"
  "    --compute_method <value>\n"
  "        manner of computing the result (CPU=cpu, GPU=gpu (default),\n"
  "        REF=reference method)\n"
  "\n"
  "    --num_proc_vector <value>\n"
  "        blocking factor to denote number of blocks used to decompose\n"
  "        the total number of vectors across processors \n"
  "        (default is the total number of procs requested)\n"
  "\n"
  "    --num_proc_field <value>\n"
  "        blocking factor to denote number of blocks used to decompose\n"
  "        each vector across processors (default is 1)\n"
  "\n"
  "    --num_proc_repl <value>\n"
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
  "        for this run (0-based, default is 0)\n"
  "\n"
  "    --stage_max <value>\n"
  "        the highest stage number of the sequence of stages to be computed\n"
  "        for this run (0-based, default is num_stage-1)\n"
  "\n"
  "    --num_phase <value>\n"
  "        the number of phases the computation is divided into\n"
  "        (default is 1) (available for 2-way case only)\n"
  "\n"
  "    --phase_min <value>\n"
  "        the lowest phase number of the sequence of phases to be computed\n"
  "        for this run (0-based, default is 0)\n"
  "\n"
  "    --phase_max <value>\n"
  "        the highest phase number of the sequence of phases to be computed\n"
  "        for this run (0-based, default is num_phase-1)\n"
  "\n"
  "    --input_file <value>\n"
  "        string denoting the filename or pathname file\n"
  "        used to store input vectors.  If this option not present,\n"
  "        a synthetic test case is run.\n"
  "\n"
  "    --problem_type <value>\n"
  "        the kind of synthetic test case to run. Allowed choices are\n"
  "        analytic (default) or random\n"
  "\n"
  "    --input_file <value>\n"
  "        string denoting the filename or pathname file\n"
  "        used to store input vectors.  If this option not present,\n"
  "        a synthetic test case is run.\n"
  "\n"
  "    --output_file_stub <value>\n"
  "        string denoting the filename or pathname stub of filenames\n"
  "        used to store result metrics.  Metric values are stored in files\n"
  "        whose names are formed by appending a unique identifier\n"
  "        (e.g., processor number) to the end of this string.  If this\n"
  "        option is absent, no output files are written.\n"
  "\n"
  "    --threshold <value>\n"
  "        output each result value only if its magnitude is greater than\n"
  "        this threshold.  If set negative, no thresholding is done\n"
  "        (default -1)\n"
  "\n"
  "    --checksum <value>\n"
  "        compute checksum of the metrics results (yes=yes (default), no=no)\n"
  "\n"
  "    --verbosity <value>\n"
  "       verbosity level of output (0=none, 1=some (default) 2,3=more)\n"
  "\n"
  );
  /* clang-format on */
}

//=============================================================================

void perform_run_preflight(int argc, char** argv) {

  GMEnv env_val = GMEnv_null(), *env = &env_val;;
  GMEnv_create(env, MPI_COMM_WORLD, argc, (char**)argv, NULL);

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {

    /*---Perform preliminary run on GPU since sometimes first use is slower---*/

    int num_proc = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    const char* options_template_1 =
        GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK ?
          "--num_field 1 --num_vector_local 2 "
          "--metric_type ccc "
          "--num_proc_vector %i --all2all no --num_way 2 "
          "--compute_method GPU --verbosity 0" :
          "--num_field 1 --num_vector_local 2 "
          "--metric_type czekanowski "
          "--num_proc_vector %i --all2all no --num_way 2 "
          "--compute_method GPU --verbosity 0";

    char options1[1024];
    sprintf(options1, options_template_1, num_proc);

    perform_run(options1);
  }

  GMEnv_destroy(env);
}

//=============================================================================
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

  /*---Perform preflight warmup---*/

  perform_run_preflight(argc, argv);

  /*---Perform actual run---*/

  perform_run(argc, (char**)argv, NULL);

  MPI_Finalize();
  return 0;
}

//-----------------------------------------------------------------------------
