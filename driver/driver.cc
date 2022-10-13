//-----------------------------------------------------------------------------
/*!
 * \file   driver.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#include "cstdio"
#include "cstdlib"
#include "cstddef"
#include "string.h"
#include "float.h"
#include "errno.h"

#include "unistd.h"

#include "env.hh"
#include "histograms.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "checksum.hh"
#include "compute_metrics.hh"

#include "test_problems.hh"
#include "vectors_io.hh"
#include "metrics_io.hh"
#include "driver.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Initialize driver.

Driver::Driver(CEnv& env)
  : options_(env)
  , counters_()
  , env_(env) {
}

Driver::Options::Options(CEnv& env)
  : num_field_local(0)
  , num_vector_local(0)
  , num_field_active(0)
  , num_vector_active(0)
  , is_inited_num_field_local(false)
  , is_inited_num_field_active(false)
  , is_inited_num_vector_local(false)
  , is_inited_num_vector_active(false)
  , verbosity(1)
  , stage_min(0)
  , stage_max(env.num_stage() - 1)
  , phase_min(0)
  , phase_max(env.num_phase() - 1)
  , input_file(NULL)
  , histograms_file(NULL)
  , output_file_stub(NULL)
  , problem_type(ProblemType::DEFAULT)
  , checksum(true)
  , env_(env) {
}

Driver::Counters::Counters()
  : num_incorrect(0)
  , max_incorrect_diff(0)
  , num_metric_items_computed(0)
  , num_metrics_active(0)
  , num_written(0)
  , vctime(0)
  , mctime(0)
  , cktime(0)
  , intime(0)
  , outtime(0)
  , tottime(0)
  , num_metric_items_local_computed(0)
  , num_metrics_active_local(0)
  , num_local_written(0) {
}

//-----------------------------------------------------------------------------
// Parse remaining unprocessed arguments.

void Driver::finish_parsing(int argc, char** argv) {

  errno = 0; // from std C.
  for (int i = 1; i < argc; ++i) {

    if (strcmp(argv[i], "--num_field") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for num_field.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0
                             && "Invalid setting for num_field.");
      options_.num_field_active = safe_cast<int>(value);
      options_.is_inited_num_field_active = true;
      options_.is_inited_num_field_local = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_field_local") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_,
                             i < argc && "Missing value for num_field_local.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0 &&
                    "Invalid setting for num_field_local.");
      options_.num_field_local = safe_cast<int>(value);
      options_.is_inited_num_field_local = true;
      options_.is_inited_num_field_active = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_vector") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_,
                             i < argc && "Missing value for num_vector.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0
                    && "Invalid setting for num_vector.");
      options_.num_vector_active = safe_cast<int>(value);
      options_.is_inited_num_vector_active = true;
      options_.is_inited_num_vector_local = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_vector_local") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_,
                             i < argc && "Missing value for num_vector_local.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0 &&
                    "Invalid setting for num_vector_local.");
      options_.num_vector_local = safe_cast<int>(value);
      options_.is_inited_num_vector_local = true;
      options_.is_inited_num_vector_active = false;

    //--------------------

    } else if (strcmp(argv[i], "--verbosity") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for verbosity.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0 &&
                    "Invalid setting for verbosity.");
      options_.verbosity = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--checksum") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for checksum.");
      if (strcmp(argv[i], "yes") == 0)
        options_.checksum = true;
      else if (strcmp(argv[i], "no") == 0)
        options_.checksum = false;
      else
        COMET_INSIST_INTERFACE(&env_, false && "Invalid setting for checksum.");

    //--------------------

    } else if (strcmp(argv[i], "--num_stage") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for num_stage.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 1
                    && "Invalid setting for num_stage.");
      env_.num_stage(safe_cast<int>(value));
      options_.stage_min = 0;
      options_.stage_max = env_.num_stage() - 1;

    //--------------------

    } else if (strcmp(argv[i], "--stage_min") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for stage_min.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0
                    && "Invalid setting for stage_min.");
      options_.stage_min = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--stage_max") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for stage_max.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value < env_.num_stage()
                    && "Invalid setting for stage_max.");
      options_.stage_max = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--num_phase") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for num_phase.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 1
                    && "Invalid setting for num_phase.");
      env_.num_phase(safe_cast<int>(value));
      options_.phase_min = 0;
      options_.phase_max = env_.num_phase() - 1;

    //--------------------

    } else if (strcmp(argv[i], "--phase_min") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for phase_min.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0
                    && "Invalid setting for phase_min.");
      options_.phase_min = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--phase_max") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for phase_max.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value < env_.num_phase()
                    && "Invalid setting for phase_max.");
      options_.phase_max = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--input_file") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_,
                             i < argc && "Missing value for input_file.");
      options_.input_file = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--histograms_file") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_,
                             i < argc && "Missing value for histograms_file.");
      options_.histograms_file = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--output_file_stub") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_,
                             i < argc && "Missing value for output_file_stub.");
      options_.output_file_stub = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--problem_type") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_,
                             i < argc && "Missing value for problem_type.");
      if (strcmp(argv[i], "random") == 0)
        options_.problem_type = ProblemType::RANDOM;
      else if (strcmp(argv[i], "analytic") == 0)
        options_.problem_type = ProblemType::ANALYTIC;
      else
        COMET_INSIST_INTERFACE(&env_,
                               false && "Invalid setting for problem_type.");

    } else if (strcmp(argv[i], "--threshold") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--ccc_multiplier") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--duo_multiplier") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--sparse") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--tc") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--num_tc_steps") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--metrics_shrink") == 0) {
      ++i; // processed by CEnv.
    } else if (strcmp(argv[i], "--fastnodes") == 0) {
      // processed by main.
    } else if (strcmp(argv[i], "--nopreflight") == 0) {
      // processed by main.
    } else {
    //----------
      if (env_.proc_num() == 0)
        fprintf(stderr, "Invalid argument \"%s\". ", argv[i]);
      COMET_INSIST_INTERFACE(&env_, false && "Error: argument not recognized.");
    //----------
    } // if/else

  } // for i

  // Checks.

  COMET_INSIST_INTERFACE(&env_, (options_.is_inited_num_field_local ||
                                 options_.is_inited_num_field_active)
                && "Error: must set either num_field_local or num_field.");
  COMET_INSIST_INTERFACE(&env_, (options_.is_inited_num_vector_local ||
                                 options_.is_inited_num_vector_active)
                && "Error: must set either num_vector_local or num_vector.");
}

//-----------------------------------------------------------------------------

void Driver::set_vectors(GMVectors& vectors) {
  if (options_.input_file)
    VectorsIO::read(vectors, options_.input_file, env_);
  else
    TestProblem::set_vectors_synthetic(&vectors, options_.problem_type,
                                       options_.verbosity, &env_);

  if (options_.verbosity > 2)
    VectorsIO::print(vectors, env_);
}

//-----------------------------------------------------------------------------
// Print a line of output to summarize result of run.

void Driver::print_output_sync(Checksum& cksum) {

  // Perform operations that may include allreduce.

  const double ops = env_.ops();
  const double ops_gemm = env_.ops_gemm();
  const double gemmtime_sum = env_.gemmtime_sum();
  const size_t cpu_mem_max = env_.cpu_mem_max();
  const size_t gpu_mem_max = env_.gpu_mem_max();

  if (!do_print())
    return;

  if (cksum.computing_checksum()) {
    printf("metrics checksum ");
    cksum.print(env_);
    printf(" ");
  }

  printf("ctime %.6f", env_.ctime());

  printf(" ops %e", ops);
  if (env_.ctime() > 0) {
    printf(" ops_rate %e", ops / env_.ctime());
    printf(" ops_rate/proc %e", ops / (env_.ctime() * env_.num_proc()) );
  }

  if (gemmtime_sum > 0)
    printf(" gemmrate/proc %e", ops_gemm / gemmtime_sum);

  printf(" vcmp %e", env_.vec_compares());
  printf(" vacmp %e", env_.vec_active_compares());
  if (options_.output_file_stub)
    printf(" vcmpout %e", (double)counters_.num_written);

  printf(" cmp %e", env_.entry_compares());
  printf(" acmp %e", env_.entry_active_compares());

  printf(" ecmp %e", env_.metric_compares());
  printf(" eacmp %e", env_.metric_active_compares());
  if (env_.ctime() > 0) {
    printf(" ecmp_rate %e", env_.metric_compares() / env_.ctime());
    printf(" ecmp_rate/proc %e", env_.metric_compares() /
      (env_.ctime() * env_.num_proc()) );
  }

  printf(" ment %e", (double)env_.metric_entries());
  printf(" mentc %e", (double)env_.metric_entries_computed());

  if (cksum.computing_checksum()) {
    printf(" me %.0f", cksum.num());
    printf(" mezero %.0f", cksum.num_zero());
    if (cksum.num() > 0)
      printf(" fracnonzero %.9f",
        (cksum.num()-cksum.num_zero()) / cksum.num());
  }

  if (env_.is_shrink())
    printf(" shrink %e", env_.shrink_achieved());

  printf(" vctime %.6f", counters_.vctime);
  printf(" mctime %.6f", counters_.mctime);
  if (cksum.computing_checksum())
    printf(" cktime %.6f", counters_.cktime);
  printf(" intime %.6f", counters_.intime);
  printf(" outtime %.6f", counters_.outtime);

  printf(" cpumem %e", (double)cpu_mem_max);
  printf(" gpumem %e", (double)gpu_mem_max);

  printf(" tottime %.6f", counters_.tottime);

  printf(" prec %s", env_.is_double_prec() ? "double" : "single");

  printf(" build %s", BuildHas::DEBUG ? "debug" : "release");

  if (env_.tc() != env_.tc_eff())
    printf(" tc_eff %i", env_.tc_eff());

  if (env_.is_shrink())
    printf(" is_shrink %s", "yes");

  printf("\n");
}

//=============================================================================
// Functions to perform a single metrics computation run.

void Driver::perform_run(int argc, char** argv, MPI_Comm base_comm) {
  COMET_INSIST(argc > 0 && argv);

  Checksum cksum;
  const char* const description = NULL;
  CEnv* env = NULL;

  perform_run_(cksum, argc, argv, description, base_comm, env);
}

//-----------------------------------------------------------------------------

void Driver::perform_run(const char* const options_str) {
  COMET_INSIST(options_str);

  Checksum cksum;
  MPI_Comm base_comm = MPI_COMM_WORLD;
  CEnv* env = NULL;

  perform_run(cksum, options_str, base_comm, env);
}

//-----------------------------------------------------------------------------

void Driver::perform_run(const char* const options_str, MPI_Comm base_comm,
                         CEnv& env) {
  COMET_INSIST(options_str);

  Checksum cksum;

  perform_run(cksum, options_str, base_comm, &env);
}

//-----------------------------------------------------------------------------

void Driver::perform_run(Checksum& cksum, const char* const options_str,
                         MPI_Comm base_comm, CEnv* env) {
  COMET_INSIST(options_str);

  // Convert options string to argc/argv.

  size_t len = strlen(options_str);
  char argstring[len+1];
  char* argv[len+1];
  int argc = 0;
  strcpy(argstring, options_str);
  CEnv::create_args(argstring, &argc, argv);

  const char* const description = options_str;

  return perform_run_(cksum, argc, argv, description, base_comm, env);
}

//-----------------------------------------------------------------------------

void Driver::perform_run_(Checksum& cksum_result, int argc, char** argv,
                          const char* const description, MPI_Comm base_comm,
                          CEnv* env_in) {

  // Initialize environment if needed.

  CEnv* env_local = NULL;

  if (!env_in) {
    env_local = new CEnv(base_comm, argc, argv, description);
    if (! env_local->is_proc_active()) {
      delete env_local;
      return;
    }
  }

  CEnv* const env = env_in ? env_in : env_local;

  perform_run_(cksum_result, argc, argv, base_comm, *env);

  if (env_local)
    delete env_local;
}

//-----------------------------------------------------------------------------
// Workhorse function to actually perform run.

void Driver::perform_run_(Checksum& cksum_result, int argc, char** argv,
                          MPI_Comm base_comm, CEnv& env) {
  COMET_INSIST(argc > 0 && argv);

  if (! env.is_proc_active())
    return;

  Driver::Timer timer_total(env);

  // Parse remaining unprocessed arguments.

  Driver driver(env);
  driver.finish_parsing(argc, argv);
  Driver::Options& options_ = driver.options_;
  Driver::Counters& counters_ = driver.counters_;

  // Set up parallel decomp for vectors, metrics.

  Driver::Timer timer(env);
  GMDecompMgr dm = GMDecompMgr_null();
  GMDecompMgr_create(&dm,
    options_.is_inited_num_field_local,
    options_.is_inited_num_vector_local,
    options_.is_inited_num_field_local ? options_.num_field_local
                                       : options_.num_field_active,
    options_.is_inited_num_vector_local ? options_.num_vector_local
                                        : options_.num_vector_active,
    env.data_type_vectors(), &env);
  timer.add_elapsed(counters_.vctime);

  // Allocate vectors.

  timer.start();
  GMVectors vectors;
  vectors.create(env.data_type_vectors(), dm, env);
  timer.add_elapsed(counters_.vctime);

  // Set vectors.

  timer.start();
  driver.set_vectors(vectors);
  timer.add_elapsed(counters_.intime);

  // More initializations.

  Checksum cksum(options_.checksum);
  Checksum cksum_local(options_.checksum);

  // Initialize output.

  timer.start();
  MetricsIO metrics_io(options_.output_file_stub, options_.verbosity, env);
  Histograms histograms(options_.histograms_file, env);
  dm.attach_histograms(&histograms);
  timer.add_elapsed(counters_.outtime);

  // Initialize metrics mem, compute metrics.

  timer.start();
  MetricsMem metrics_mem(&env);
  ComputeMetrics compute_metrics(dm, env);
  timer.add_elapsed(counters_.mctime);

  //--------------------
  // Begin loops over phases, stages.
  //--------------------

  for (int phase_num=options_.phase_min; phase_num<=options_.phase_max;
       ++phase_num) {
      env.phase_num(phase_num);

  for (int stage_num=options_.stage_min; stage_num<=options_.stage_max;
       ++stage_num) {
    env.stage_num(stage_num);

    // Set up metrics object to capture results.

    timer.start(); // 
    GMMetrics metrics = GMMetrics_null();
    GMMetrics_create(&metrics, env.data_type_metrics(), &dm, &metrics_mem,
                     &env);
    timer.add_elapsed(counters_.mctime);

    // Calculate metrics.

    compute_metrics.compute(metrics, vectors);
    counters_.num_metric_items_local_computed +=
      metrics.num_metric_items_local_computed;
    counters_.num_metrics_active_local += metrics.num_metrics_active_local;

    // Output results.

    timer.start();
    metrics_io.write(metrics);
    if (BuildHas::DEBUG)
      metrics_io.check_file(metrics);
    timer.add_elapsed(counters_.outtime);

    // Check correctness.

    timer.start();
    if (options_.checksum)
      TestProblem::check_metrics(&metrics, driver, &env);
    timer.add_elapsed(counters_.cktime);

    // Compute checksum.

    timer.start();
    if (options_.checksum)
      Checksum::compute(cksum, cksum_local, metrics, env);
    timer.add_elapsed(counters_.cktime);

    // Delete metrics object.

    timer.start();
    GMMetrics_destroy(&metrics, &env);
    timer.add_elapsed(counters_.mctime);

    // Do output.

    if (driver.do_print()) {
      if (env.num_phase() > 1 && env.num_stage() > 1)
        printf("Completed phase %i stage %i\n",
               env.phase_num(), env.stage_num());
      else if (env.num_phase() > 1)
        printf("Completed phase %i\n", env.phase_num());
      else if (env.num_stage() > 1)
        printf("Completed stage %i\n", env.stage_num());
    } // do_print
    driver.fflush_sync_();

  } // for stage_num

  } // for phase_num

  //--------------------
  // End loops over phases, stages.
  //--------------------

  // Finalize metrics mem.

  timer.start();
  metrics_mem.terminate();
  compute_metrics.terminate();
  timer.add_elapsed(counters_.mctime);

  // Finalize output.

  timer.start();
  counters_.num_local_written += metrics_io.num_written();
  metrics_io.terminate();
  histograms.finalize();
  histograms.output();
  if (options_.is_all_phase_all_stage())
    histograms.check(dm.num_vector_active);
  histograms.terminate();
  timer.add_elapsed(counters_.outtime);

  // Perform some checks.

  COMET_MPI_SAFE_CALL(MPI_Allreduce(
    &counters_.num_metric_items_local_computed,
    &counters_.num_metric_items_computed, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
    env.comm_repl_vector()));

  COMET_MPI_SAFE_CALL(MPI_Allreduce(
    &counters_.num_metrics_active_local,
    &counters_.num_metrics_active, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
    env.comm_repl_vector()));

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&counters_.num_local_written,
                      &counters_.num_written, 1, MPI_UNSIGNED_LONG_LONG,
                      MPI_SUM, env.comm_repl_vector()));

  if (env.all2all() && options_.is_all_phase_all_stage()) {
    const size_t num_metrics_expected =
      utils::nchoosek(dm.num_vector, env.num_way());

    const size_t num_metrics_active_expected =
      utils::nchoosek(dm.num_vector_active, env.num_way());

    COMET_INSIST(counters_.num_metrics_active == num_metrics_active_expected &&
                 "Inconsistent metrics count.");

    COMET_INSIST((env.num_metric_items_per_metric() * num_metrics_expected ==
                  counters_.num_metric_items_computed || env.is_shrink()) &&
                 "Inconsistent metrics count.");
  }

  // Deallocate vectors.

  timer.start();
  GMVectors_destroy(&vectors, &env);
  GMDecompMgr_destroy(&dm, &env);
  timer.add_elapsed(counters_.vctime);

  COMET_INSIST(env.cpu_mem_local() == 0 && "Memory leak detected.");
  COMET_INSIST(env.gpu_mem_local() == 0 && "Memory leak detected.");

  // Record the total time.

  counters_.tottime = timer_total.elapsed();

  // Output run information.

  driver.print_output_sync(cksum);
  driver.fflush_sync_();

  // Output a local checksum, if needed for testing purposes.

  const bool do_output_local_checksum = false;

  if (do_output_local_checksum) {
    if (options_.checksum && env.is_proc_active() && options_.verbosity > 0) {
      printf("local checksum: ");
      cksum_local.print(env);
      printf("\n");
    }
    driver.fflush_sync_();
  }

  // Validation: check for any wrong answers.

  if (counters_.num_incorrect) {
    const size_t hnlen = 256;
    char hn[hnlen];
    gethostname(hn, hnlen);
    int rank = 0;
    COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

    fprintf(stderr, "Error: incorrect results found.  num_incorrect  %zu  "
           "max_incorrect_diff  %e  hostname  %s  rank  %i\n",
           counters_.num_incorrect, counters_.max_incorrect_diff, hn, rank);
  }
  driver.fflush_sync_();

  COMET_INSIST(0 == counters_.num_incorrect && "Incorrect results found.");

  // Finalize.

  cksum_result.copy(cksum);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
