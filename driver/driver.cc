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
// Parse remaining unprocessed arguments.

void finish_parsing(int argc, char** argv, DriverOptions* do_, CEnv* env) {
  errno = 0;
  int i = 0;
  for (i = 1; i < argc; ++i) {
    //----------
    if (strcmp(argv[i], "--num_field") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_field.");
      const long num_field = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && num_field >= 0
                    && "Invalid setting for num_field.");
      do_->num_field_active = num_field;
      do_->num_field_active_initialized = true;
      do_->num_field_local_initialized = false;
    //----------
    } else if (strcmp(argv[i], "--num_field_local") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_field_local.");
      const long num_field_local = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && num_field_local >= 0 &&
                    (long)(int)num_field_local == num_field_local &&
                    "Invalid setting for num_field_local.");
      do_->num_field_local = num_field_local;
      do_->num_field_local_initialized = true;
      do_->num_field_active_initialized = false;
    //----------
    } else if (strcmp(argv[i], "--num_vector") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_vector.");
      const long num_vector = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && num_vector >= 0
                    && "Invalid setting for num_vector.");
      do_->num_vector_active = num_vector;
      do_->num_vector_active_initialized = true;
      do_->num_vector_local_initialized = false;
    //----------
    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_vector_local.");
      const long num_vector_local = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && num_vector_local >= 0 &&
                    (long)(int)num_vector_local == num_vector_local &&
                    "Invalid setting for num_vector_local.");
      do_->num_vector_local = num_vector_local;
      do_->num_vector_local_initialized = true;
      do_->num_vector_active_initialized = false;
    //----------
    } else if (strcmp(argv[i], "--verbosity") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for verbosity.");
      const long verbosity = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && verbosity >= 0 &&
                    "Invalid setting for verbosity.");
      do_->verbosity = verbosity;
      //--------------------
    } else if (strcmp(argv[i], "--checksum") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for checksum.");
      if (strcmp(argv[i], "yes") == 0) {
        do_->checksum = true;
      } else if (strcmp(argv[i], "no") == 0) {
        do_->checksum = false;
      } else {
        COMET_INSIST_INTERFACE(env, false && "Invalid setting for checksum.");
      }
    //----------
    } else if (strcmp(argv[i], "--num_stage") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_stage.");
      const long num_stage = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && num_stage >= 1
                    && (long)(int)num_stage == num_stage
                    && "Invalid setting for num_stage.");
      env->num_stage(num_stage);
      do_->stage_min_0based = 0;
      do_->stage_max_0based = env->num_stage() - 1;
    //----------
    } else if (strcmp(argv[i], "--stage_min") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for stage_min.");
      const long stage_min_0based = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && stage_min_0based >= 0
                    && (long)(int)stage_min_0based == stage_min_0based
                    && "Invalid setting for stage_min.");
      do_->stage_min_0based = stage_min_0based;
    //----------
    } else if (strcmp(argv[i], "--stage_max") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for stage_max.");
      const long stage_max_0based = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && stage_max_0based < env->num_stage()
                    && (long)(int)stage_max_0based == stage_max_0based
                    && "Invalid setting for stage_max.");
      do_->stage_max_0based = stage_max_0based;
    //----------
    } else if (strcmp(argv[i], "--num_phase") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_phase.");
      const long num_phase = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && num_phase >= 1
                    && (long)(int)num_phase == num_phase
                    && "Invalid setting for num_phase.");
      env->num_phase(num_phase);
      do_->phase_min_0based = 0;
      do_->phase_max_0based = env->num_phase() - 1;
    //----------
    } else if (strcmp(argv[i], "--phase_min") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for phase_min.");
      const long phase_min_0based = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && phase_min_0based >= 0
                    && (long)(int)phase_min_0based == phase_min_0based
                    && "Invalid setting for phase_min.");
      do_->phase_min_0based = phase_min_0based;
    //----------
    } else if (strcmp(argv[i], "--phase_max") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for phase_max.");
      const long phase_max_0based = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && phase_max_0based < env->num_phase()
                    && (long)(int)phase_max_0based == phase_max_0based
                    && "Invalid setting for phase_max.");
      do_->phase_max_0based = phase_max_0based;
    //----------
    } else if (strcmp(argv[i], "--input_file") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for input_file.");
      do_->input_file_path = argv[i];
    //----------
    } else if (strcmp(argv[i], "--output_file_stub") == 0) {
    //----------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for output_file_stub.");
      do_->metrics_file_path_stub = argv[i];
      //--------------------
    } else if (strcmp(argv[i], "--problem_type") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for problem_type.");
      if (strcmp(argv[i], "random") == 0) {
        do_->problem_type = GM_PROBLEM_TYPE_RANDOM;
      } else if (strcmp(argv[i], "analytic") == 0) {
        do_->problem_type = GM_PROBLEM_TYPE_ANALYTIC;
      } else {
        COMET_INSIST_INTERFACE(env, false && "Invalid setting for problem_type.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--detailed_output") == 0) {
      //--------------------
      ++i;
      if (strcmp(argv[i], "yes") == 0) {
        do_->detailed_output = true;
      } else if (strcmp(argv[i], "no") == 0) {
        do_->detailed_output = false;
      } else {
        COMET_INSIST_INTERFACE(env, false && "Invalid setting for detailed_output.");
      }
     //----------
    } else if (strcmp(argv[i], "--threshold") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--ccc_multiplier") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--duo_multiplier") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--sparse") == 0) {
      ++i; // processed elsewhere by CEnv.
    } else if (strcmp(argv[i], "--fastnodes") == 0) {
      // optionally processed by caller.
    } else if (strcmp(argv[i], "--nopreflight") == 0) {
      // optionally processed by caller.
    } else if (strcmp(argv[i], "--tc") == 0) {
      ++i; // optionally processed by caller.
    } else if (strcmp(argv[i], "--num_tc_steps") == 0) {
      ++i; // optionally processed by caller.
    } else if (strcmp(argv[i], "--num_kernel") == 0) {
      ++i; // optionally processed by caller.
    } else if (strcmp(argv[i], "--metrics_shrink") == 0) {
      ++i; // optionally processed by caller.
    } else if (strcmp(argv[i], "--print_details") == 0) {
      ++i; // optionally processed by caller.
    } else if (strcmp(argv[i], "--sync_time") == 0) {
      ++i; // optionally processed by caller.
    } else {
    //----------
      if (env->proc_num() == 0) {
        fprintf(stderr, "Invalid argument \"%s\". ", argv[i]);
      }
      COMET_INSIST_INTERFACE(env, false && "Error: argument not recognized.");
    //----------
    } // if/else

  } // for i

  COMET_INSIST_INTERFACE(env, (do_->num_field_local_initialized ||
                do_->num_field_active_initialized)
                && "Error: must set num_field_local or num_field.");
  COMET_INSIST_INTERFACE(env, (do_->num_vector_local_initialized ||
                do_->num_vector_active_initialized)
                && "Error: must set num_vector_local or num_vector.");
}

//-----------------------------------------------------------------------------

void set_vectors(GMVectors* vectors, DriverOptions* do_, CEnv* env) {
  COMET_INSIST(vectors && do_ && env);

  if (do_->input_file_path != NULL) {
    VectorsIO::read(*vectors, do_->input_file_path, *env);
  } else {
    set_vectors_synthetic(vectors, do_->problem_type, do_->verbosity, env);
  }
}

//-----------------------------------------------------------------------------

void print_output(bool do_print,
                  bool do_detailed,
                  Checksum& cksum,
                  CEnv& env,
                  char* metrics_file_path_stub,
                  size_t num_written,
                  double vctime,
                  double mctime,
                  double cmtime,
                  double cktime,
                  double intime,
                  double outtime,
                  double looptime,
                  double tottime) {

  const double ops = env.ops();
  const double simops = env.simops();
  const size_t cpu_mem_max = env.cpu_mem_max();
  const size_t gpu_mem_max = env.gpu_mem_max();

  if (!do_print)
    return;

  if(do_detailed) printf("\nOutput:\n");
  
  if (cksum.computing_checksum()) {
    printf("metrics checksum ");
    cksum.print(env);
    printf(" ");
  }

  printf("ctime %.6f", env.ctime());

  printf(" ops %e", ops);
  if (env.ctime() > 0) {
    printf(" ops_rate %e", ops / env.ctime());
    printf(" ops_rate/proc %e", ops / (env.ctime() * env.num_proc()) );
  }

  printf(" vcmp %e", env.vec_compares());
  if (metrics_file_path_stub) {
    printf(" vcmpout %e", (double)num_written);
  }

  printf(" cmp %e", env.entry_compares());

  printf(" ecmp %e", env.metric_compares());
  if (env.ctime() > 0) {
    printf(" ecmp_rate %e", env.metric_compares() / env.ctime());
    printf(" ecmp_rate/proc %e", env.metric_compares() /
      (env.ctime() * env.num_proc()) );
  }

  printf(" ment %e", (double)env.metric_entries());
  printf(" mentc %e", (double)env.metric_entries_computed());

  if (cksum.computing_checksum()) {
    printf(" me %.0f", cksum.num());
    printf(" mezero %.0f", cksum.num_zero());
    if (cksum.num() > 0) {
      printf(" fracnonzero %.9f",
        (cksum.num()-cksum.num_zero()) / cksum.num());
    }
  }

  if (env.is_shrink())
    printf(" shrink %e", env.shrink_achieved());

  printf(" vctime %.6f", vctime);
  printf(" mctime %.6f", mctime);
  if (cksum.computing_checksum()) {
    printf(" cktime %.6f", cktime);
  }
  printf(" intime %.6f", intime);
  printf(" outtime %.6f", outtime);

  printf(" cpumem %e", (double)cpu_mem_max);
  printf(" gpumem %e", (double)gpu_mem_max);

  printf(" tottime %.6f", tottime);

  printf(" prec %s", env.is_double_prec() ? "double" : "single");

  printf(" build %s", BuildHas::DEBUG ? "debug" : "release");

  if (env.tc() != env.tc_eff()) {
    printf(" tc_eff %i", env.tc_eff());
  }

  if (env.is_shrink()) {
    printf(" is_shrink %s", "yes");
  }

  if(do_detailed) printf("\n");

  if(do_detailed) {
    // More readable runtime output
    printf("\nDetailed Output:\n");
    printf("GEMM:\n"
           "Pre-GEMM time:        %.6f\n"
           "GEMM time:            %.6f\n"
           "Post-GEMM time:       %.6f\n",
            env.pregemmtime(), env.gemmtime(), env.postgemmtime());

    printf("\nComputeMetrics:\n"
           "Nums start time:      %.6f\n"
           "Nums wait time:       %.6f\n"
           "Combine time:         %.6f\n",
           env.numsstarttime(), env.numswaittime(), env.combinetime());

    printf("\nDriver:\n"
           "Input time:           %.6f\n"
           "Compute metric time:  %.6f\n"
           "Output time:          %.6f\n"
           "Compute metric total: %.6f\n"
           "Loop time:            %.6f\n"
           "Total time:           %.6f\n",
           intime, cmtime, outtime, env.ctime(), looptime, tottime);

    double tops = ops/(1024.0*1024.0*1024.0*1024.0);
    double gemmops = 0.0, cmops = 0.0, cops = 0.0;
    if(env.gemmtime()>0.0) gemmops = tops/env.gemmtime();
    if(cmtime>0.0)         cmops   = tops/cmtime;
    if(env.ctime()>0.0)    cops    = tops/env.ctime();

    printf("\nOps (Algorithm per Byte):\n"
           "Ops:                  %e\n"
           "TOps:                 %.2f\n"
           "GEMM TOps:            %.2f\n"
           "Metric TOps:          %.2f\n"
           "Metric total TOps:    %.2f\n\n",
           ops, tops, gemmops, cmops, cops);

    double tsimops = simops/(1024.0*1024.0*1024.0*1024.0);
    gemmops = 0.0; cmops = 0.0; cops = 0.0;
    if(env.gemmtime()>0.0) gemmops = tsimops/env.gemmtime();
    if(cmtime>0.0)         cmops   = tsimops/cmtime;
    if(env.ctime()>0.0)    cops    = tsimops/env.ctime();

    printf("\nSimulation Ops (Simulation per Bit):\n"
           "SimOps:               %e\n"
           "TSimOps:              %.2f\n"
           "GEMM TSimOps:         %.2f\n"
           "Metric TSimOps:       %.2f\n"
           "Metric total TSimOps: %.2f\n\n",
           simops, tsimops, gemmops, cmops, cops);
  }
}

//=============================================================================
// Perform a single metrics computation run.

void perform_run(const char* const options, MPI_Comm base_comm, CEnv* env) {
  COMET_INSIST(options);

  comet::Checksum cksum;

  perform_run(cksum, options, base_comm, env);
}

//-----------------------------------------------------------------------------

void perform_run(comet::Checksum& cksum, const char* const options,
                 MPI_Comm base_comm, CEnv* env) {
  COMET_INSIST(options);

  // Convert options string to args.

  size_t len = strlen(options);
  char argstring[len+1];
  char* argv[len+1];
  int argc = 0;
  strcpy(argstring, options);
  CEnv::create_args(argstring, &argc, argv);

  return perform_run(cksum, argc, argv, options, base_comm, env);
}

//-----------------------------------------------------------------------------

void perform_run(int argc, char** argv, const char* const description,
                            MPI_Comm base_comm, CEnv* env) {

  comet::Checksum cksum;

  perform_run(cksum, argc, argv, description, base_comm, env);
}

//-----------------------------------------------------------------------------

void perform_run(comet::Checksum& cksum_result, int argc, char** argv,
                 const char* const description,
                 MPI_Comm base_comm, CEnv* env_in) {

  // Initialize environment.

  CEnv* env_local = NULL;

  if (!env_in) {
    env_local = new CEnv(base_comm, argc, argv, description);
    if (! env_local->is_proc_active()) {
      delete env_local;
      return;
    }
  }

  CEnv* const env = env_in ? env_in : env_local;

  if(env->print_details()) printf("Starting CoMet run\n");
  double total_time_beg = env->synced_time();

  // Parse remaining unprocessed arguments.

  DriverOptions do_ = {0};
  do_.num_field_local_initialized = false;
  do_.num_field_active_initialized = false;
  do_.num_vector_local_initialized = false;
  do_.num_vector_active_initialized = false;
  do_.verbosity = 1;
  do_.stage_min_0based = 0;
  do_.stage_max_0based = env->num_stage() - 1;
  do_.phase_min_0based = 0;
  do_.phase_max_0based = env->num_phase() - 1;
  do_.input_file_path = NULL;
  do_.metrics_file_path_stub = NULL;
  do_.problem_type = problem_type_default();
  do_.checksum = true;
  do_.num_incorrect = 0;
  do_.max_incorrect_diff = 0.;

  finish_parsing(argc, argv, &do_, env);

  // Set up parallel deomp for vectors, metrics.

  double vctime = 0;
  double time_beg = env->synced_time();
  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm,
    do_.num_field_local_initialized,
    do_.num_vector_local_initialized,
    do_.num_field_local_initialized ? do_.num_field_local
                                    : do_.num_field_active,
    do_.num_vector_local_initialized ? do_.num_vector_local
                                     : do_.num_vector_active,
    env->data_type_vectors(), env);
  env->stream_compute(); //FIX
  double time_end = env->synced_time();
  vctime += time_end - time_beg;

  //TODO: possibly replace this with stuff from dm
  if (do_.num_vector_local_initialized) {
    do_.num_vector = do_.num_vector_local *
      (size_t)env->num_proc_vector();
    do_.num_vector_active = do_.num_vector;
  } else {
    // Pad up so that every proc has same number of vectors.
    do_.num_vector_local = gm_num_vector_local_required(
      utils::ceil(do_.num_vector_active, (size_t)env->num_proc_vector()), env);
    do_.num_vector = do_.num_vector_local *
      (size_t)env->num_proc_vector();
  }

  if (do_.num_field_local_initialized) {
    do_.num_field = do_.num_field_local * (size_t) env->num_proc_field();
    do_.num_field_active = do_.num_field;
  } else {
    // Pad up so that every proc has same number of fields.
    do_.num_field_local = utils::ceil(
        do_.num_field_active, (size_t)env->num_proc_field());
    do_.num_field = do_.num_field_local * (size_t) env->num_proc_field();
  }

  const bool do_print = env->is_proc_active() &&
     env->proc_num() == 0 && do_.verbosity > 0;
  const bool do_detailed = do_print && do_.detailed_output;

  // Allocate vectors.

  time_beg = env->synced_time();
  GMVectors vectors_value = GMVectors_null(), *vectors = &vectors_value;
  GMVectors_create(vectors, env->data_type_vectors(), dm, env);
  time_end = env->synced_time();
  vctime += time_end - time_beg;

  // Set vectors.

  double intime = 0;
  time_beg = env->synced_time();

  set_vectors(vectors, &do_, env);

  time_end = env->synced_time();
  intime += time_end - time_beg;

  // More initializations.

  comet::Checksum cksum(do_.checksum);
  comet::Checksum cksum_local(do_.checksum);

  double outtime = 0;
  double mctime = 0;
  double cmtime = 0;
  double cktime = 0;
  double looptime = 0;
  double time_beg_loop = 0, time_end_loop = 0;
  size_t num_metric_items_local_computed = 0;
  size_t num_local_written = 0;

  // Open output files.

  {
  time_beg = env->synced_time();
  MetricsIO metrics_io(do_.metrics_file_path_stub, do_.verbosity, *env);
  time_end = env->synced_time();
  outtime += time_end - time_beg;

  {
  MetricsMem metrics_mem(env);

  if(env->print_details()) printf("Setting up ComputeMetrics\n");
  ComputeMetrics compute_metrics(*dm, *env);
  if(env->print_details()) printf("Done setting up ComputeMetrics\n");

  // Loops over phases, stages.
  time_beg_loop = env->synced_time();  
  for (int phase_num=do_.phase_min_0based; phase_num<=do_.phase_max_0based;
       ++phase_num) {
      env->phase_num(phase_num);

    for (int stage_num=do_.stage_min_0based; stage_num<=do_.stage_max_0based;
         ++stage_num) {
      env->stage_num(stage_num);

      // Set up metrics container for results.
      if(env->print_details()) printf("Setting up metrics phase=%d/%d stage=%d/%d\n",phase_num,
        do_.phase_max_0based,stage_num,do_.stage_max_0based);
      time_beg = env->synced_time();
      GMMetrics metrics_value = GMMetrics_null(), *metrics = &metrics_value;
      GMMetrics_create(metrics, env->data_type_metrics(), dm,
                       &metrics_mem, env);
      time_end = env->synced_time();
      mctime += time_end - time_beg;

      // Calculate metrics.
      if(env->print_details()) printf("Computing metrics\n");
      time_beg = env->synced_time();
      compute_metrics.compute(*metrics, *vectors);
      time_end = env->synced_time();
      if(env->print_details()) printf("Done computing metrics\n");
      cmtime += time_end - time_beg;

      num_metric_items_local_computed += metrics->num_metric_items_local_computed;

      // Output results.
      if(env->print_details()) printf("Outputting results\n");
      time_beg = env->synced_time();
      metrics_io.write(*metrics);
      if (BuildHas::DEBUG) {
        metrics_io.check_file(*metrics);
      }
      time_end = env->synced_time();
      outtime += time_end - time_beg;

      // Check correctness.
      if (do_.checksum) {
        if(env->print_details()) printf("Checking correctness\n");
        time_beg = env->synced_time();
        check_metrics(metrics, &do_, env);
        time_end = env->synced_time();
        cktime += time_end - time_beg;
      }

      // Compute checksum.
      if (do_.checksum) {
        if(env->print_details()) printf("Computing checksum\n");
        time_beg = env->synced_time();
        comet::Checksum::compute(cksum, cksum_local, *metrics, *env);
        time_end = env->synced_time();
        cktime += time_end - time_beg;
      }
      time_beg = env->synced_time();
      GMMetrics_destroy(metrics, env);
      time_end = env->synced_time();
      mctime += time_end - time_beg;

      if (do_print) {
        if (env->num_phase() > 1 && env->num_stage() > 1) {
          if(env->print_details()) printf("Completed phase %i stage %i\n",
                 env->phase_num(), env->stage_num());
        } else if (env->num_phase() > 1) {
          if(env->print_details()) printf("Completed phase %i\n",
                 env->phase_num());
        } else if (env->num_stage() > 1) {
          if(env->print_details()) printf("Completed stage %i\n",
                 env->stage_num());
        }
      }

    } // for stage

  } // for phase

  // Finalize metrics mem.

  time_beg = env->synced_time();
  }
  time_end = env->synced_time();
  mctime += time_end - time_beg;
  // Close output files.

  num_local_written += metrics_io.num_written();
  time_beg = env->synced_time();
  }
  time_end = env->synced_time();
  outtime += time_end - time_beg;
  time_end_loop = env->synced_time();
  looptime = time_end_loop - time_beg_loop;
  // Deallocate vectors.

  time_beg = env->synced_time();
  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  time_end = env->synced_time();
  vctime += time_end - time_beg;

  // Perform some checks.

  COMET_INSIST(env->cpu_mem_local() == 0);
  COMET_INSIST(env->gpu_mem_local() == 0);

  size_t num_written = 0;
  if (env->is_proc_active()) {
    size_t num_metric_items_computed = 0;
    COMET_MPI_SAFE_CALL(MPI_Allreduce(&num_metric_items_local_computed,
      &num_metric_items_computed, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
      env->comm_repl_vector()));

    COMET_MPI_SAFE_CALL(MPI_Allreduce(&num_local_written, &num_written, 1,
      MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));

    if (env->num_way() == NumWay::_2 && env->all2all() && !env->is_shrink() &&
        do_.phase_min_0based==0 && do_.phase_max_0based==env->num_phase() - 1) {
      COMET_INSIST(num_metric_items_computed ==
                   env->num_metric_items_per_metric() * (
                   (do_.num_vector) * (size_t)
                   (do_.num_vector - 1) / 2 ) );
    }

    if (env->num_way() == NumWay::_3 && env->all2all() && !env->is_shrink() &&
        do_.phase_min_0based==0 && do_.phase_max_0based==env->num_phase() - 1 &&
        do_.stage_min_0based==0 && do_.stage_max_0based==env->num_stage() - 1) {
      COMET_INSIST(num_metric_items_computed ==
                   env->num_metric_items_per_metric() * (
                     (do_.num_vector) * (size_t)
                     (do_.num_vector - 1) * (size_t)
                     (do_.num_vector - 2) / 6 ) );
    }
  }

  if(env->print_details()) printf("Finishing CoMet run\n");

  double total_time_end = env->synced_time();
  double tottime = total_time_end - total_time_beg;

  // Output run information.

  print_output(do_print, do_detailed, cksum, *env, do_.metrics_file_path_stub, num_written,
    vctime, mctime, cmtime, cktime, intime, outtime, looptime, tottime);
    
  // Output a local checksum, for testing purposes.

  if (false) {
    // One more sync before checking num_correct, to allow flush of output.
    env->synced_time();
    if (do_.checksum && env->is_proc_active() && do_.verbosity > 0) {
      printf("local checksum: ");
      cksum_local.print(*env);
      printf("\n");
    }
  }
  env->synced_time();

  // Validation: check for any wrong answers.

  if (do_.num_incorrect) {
    const size_t hnlen = 256;
    char hn[hnlen];
    gethostname(hn, hnlen);
    int rank = 0;
    COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

    fprintf(stderr, "Error: incorrect results found.  num_incorrect  %zu  "
           "max_incorrect_diff  %e  hostname  %s  rank  %i\n",
           do_.num_incorrect, do_.max_incorrect_diff, hn, rank);
  }

  COMET_INSIST(do_.num_incorrect == 0);

  if(env->print_details()) printf("Finished CoMet run\n");

  // Finalize.

  if (env_local) {
    delete env_local;
  }

  cksum_result.copy(cksum);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
