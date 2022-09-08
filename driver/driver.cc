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
// Initialize driver options.

#if 0
DriverOptions::DriverOptions(CEnv& env)
  : num_field_local(0)
  , num_vector_local(0)
  , num_field(0)
  , num_vector(0)
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
  , num_incorrect(0)
  , max_incorrect_diff(0.)
  , checksum(true) {
}
#endif

Driver::Driver(CEnv& env)
  : num_field_local(0)
  , num_vector_local(0)
  , num_field(0)
  , num_vector(0)
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
  , num_incorrect(0)
  , max_incorrect_diff(0.)
  , checksum(true)
  , num_metric_items_computed(0)
  , num_metrics_active(0)
  , num_written(0)
  , vctime(0.)
  , mctime(0.)
  , cktime(0.)
  , intime(0.)
  , outtime(0.)
  , tottime(0.)
  , num_metric_items_local_computed(0)
  , num_metrics_active_local(0)
  , num_local_written(0)
  , env_(env) {
}

//-----------------------------------------------------------------------------
// Parse remaining unprocessed arguments.

#if 0
void DriverOptions::finish_parsing(int argc, char** argv, CEnv& env) {

  errno = 0; // from std C.
  for (int i = 1; i < argc; ++i) {

    if (strcmp(argv[i], "--num_field") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for num_field.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 0
                             && "Invalid setting for num_field.");
      this->num_field_active = safe_cast<int>(value);
      this->is_inited_num_field_active = true;
      this->is_inited_num_field_local = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_field_local") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for num_field_local.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 0 &&
                    "Invalid setting for num_field_local.");
      this->num_field_local = safe_cast<int>(value);
      this->is_inited_num_field_local = true;
      this->is_inited_num_field_active = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_vector") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for num_vector.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 0
                    && "Invalid setting for num_vector.");
      this->num_vector_active = safe_cast<int>(value);
      this->is_inited_num_vector_active = true;
      this->is_inited_num_vector_local = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_vector_local") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for num_vector_local.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 0 &&
                    "Invalid setting for num_vector_local.");
      this->num_vector_local = safe_cast<int>(value);
      this->is_inited_num_vector_local = true;
      this->is_inited_num_vector_active = false;

    //--------------------

    } else if (strcmp(argv[i], "--verbosity") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for verbosity.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 0 &&
                    "Invalid setting for verbosity.");
      this->verbosity = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--checksum") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for checksum.");
      if (strcmp(argv[i], "yes") == 0)
        this->checksum = true;
      else if (strcmp(argv[i], "no") == 0)
        this->checksum = false;
      else
        COMET_INSIST_INTERFACE(&env, false && "Invalid setting for checksum.");

    //--------------------

    } else if (strcmp(argv[i], "--num_stage") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for num_stage.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 1
                    && "Invalid setting for num_stage.");
      env.num_stage(safe_cast<int>(value));
      this->stage_min = 0;
      this->stage_max = env.num_stage() - 1;

    //--------------------

    } else if (strcmp(argv[i], "--stage_min") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for stage_min.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 0
                    && "Invalid setting for stage_min.");
      this->stage_min = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--stage_max") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for stage_max.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value < env.num_stage()
                    && "Invalid setting for stage_max.");
      this->stage_max = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--num_phase") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for num_phase.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 1
                    && "Invalid setting for num_phase.");
      env.num_phase(safe_cast<int>(value));
      this->phase_min = 0;
      this->phase_max = env.num_phase() - 1;

    //--------------------

    } else if (strcmp(argv[i], "--phase_min") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for phase_min.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value >= 0
                    && "Invalid setting for phase_min.");
      this->phase_min = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--phase_max") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for phase_max.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env, 0 == errno && value < env.num_phase()
                    && "Invalid setting for phase_max.");
      this->phase_max = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--input_file") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for input_file.");
      this->input_file = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--histograms_file") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for histograms_file.");
      this->histograms_file = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--output_file_stub") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for output_file_stub.");
      this->output_file_stub = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--problem_type") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env, i < argc && "Missing value for problem_type.");
      if (strcmp(argv[i], "random") == 0)
        this->problem_type = ProblemType::RANDOM;
      else if (strcmp(argv[i], "analytic") == 0)
        this->problem_type = ProblemType::ANALYTIC;
      else
        COMET_INSIST_INTERFACE(&env, false && "Invalid setting for problem_type.");

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
      if (env.proc_num() == 0) {
        fprintf(stderr, "Invalid argument \"%s\". ", argv[i]);
      }
      COMET_INSIST_INTERFACE(&env, false && "Error: argument not recognized.");
    //----------
    } // if/else

  } // for i

  COMET_INSIST_INTERFACE(&env, (this->is_inited_num_field_local ||
                this->is_inited_num_field_active)
                && "Error: must set either num_field_local or num_field.");
  COMET_INSIST_INTERFACE(&env, (this->is_inited_num_vector_local ||
                this->is_inited_num_vector_active)
                && "Error: must set either num_vector_local or num_vector.");
}
#endif







void Driver::finish_parsing(int argc, char** argv) {

  errno = 0; // from std C.
  for (int i = 1; i < argc; ++i) {

    if (strcmp(argv[i], "--num_field") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for num_field.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0
                             && "Invalid setting for num_field.");
      this->num_field_active = safe_cast<int>(value);
      this->is_inited_num_field_active = true;
      this->is_inited_num_field_local = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_field_local") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for num_field_local.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0 &&
                    "Invalid setting for num_field_local.");
      this->num_field_local = safe_cast<int>(value);
      this->is_inited_num_field_local = true;
      this->is_inited_num_field_active = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_vector") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for num_vector.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0
                    && "Invalid setting for num_vector.");
      this->num_vector_active = safe_cast<int>(value);
      this->is_inited_num_vector_active = true;
      this->is_inited_num_vector_local = false;

    //--------------------

    } else if (strcmp(argv[i], "--num_vector_local") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for num_vector_local.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0 &&
                    "Invalid setting for num_vector_local.");
      this->num_vector_local = safe_cast<int>(value);
      this->is_inited_num_vector_local = true;
      this->is_inited_num_vector_active = false;

    //--------------------

    } else if (strcmp(argv[i], "--verbosity") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for verbosity.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0 &&
                    "Invalid setting for verbosity.");
      this->verbosity = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--checksum") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for checksum.");
      if (strcmp(argv[i], "yes") == 0)
        this->checksum = true;
      else if (strcmp(argv[i], "no") == 0)
        this->checksum = false;
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
      this->stage_min = 0;
      this->stage_max = env_.num_stage() - 1;

    //--------------------

    } else if (strcmp(argv[i], "--stage_min") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for stage_min.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0
                    && "Invalid setting for stage_min.");
      this->stage_min = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--stage_max") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for stage_max.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value < env_.num_stage()
                    && "Invalid setting for stage_max.");
      this->stage_max = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--num_phase") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for num_phase.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 1
                    && "Invalid setting for num_phase.");
      env_.num_phase(safe_cast<int>(value));
      this->phase_min = 0;
      this->phase_max = env_.num_phase() - 1;

    //--------------------

    } else if (strcmp(argv[i], "--phase_min") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for phase_min.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value >= 0
                    && "Invalid setting for phase_min.");
      this->phase_min = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--phase_max") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for phase_max.");
      const auto value = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(&env_, 0 == errno && value < env_.num_phase()
                    && "Invalid setting for phase_max.");
      this->phase_max = safe_cast<int>(value);

    //--------------------

    } else if (strcmp(argv[i], "--input_file") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for input_file.");
      this->input_file = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--histograms_file") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for histograms_file.");
      this->histograms_file = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--output_file_stub") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for output_file_stub.");
      this->output_file_stub = argv[i];

    //--------------------

    } else if (strcmp(argv[i], "--problem_type") == 0) {

      ++i;
      COMET_INSIST_INTERFACE(&env_, i < argc && "Missing value for problem_type.");
      if (strcmp(argv[i], "random") == 0)
        this->problem_type = ProblemType::RANDOM;
      else if (strcmp(argv[i], "analytic") == 0)
        this->problem_type = ProblemType::ANALYTIC;
      else
        COMET_INSIST_INTERFACE(&env_, false && "Invalid setting for problem_type.");

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
      if (env_.proc_num() == 0) {
        fprintf(stderr, "Invalid argument \"%s\". ", argv[i]);
      }
      COMET_INSIST_INTERFACE(&env_, false && "Error: argument not recognized.");
    //----------
    } // if/else

  } // for i

  COMET_INSIST_INTERFACE(&env_, (this->is_inited_num_field_local ||
                this->is_inited_num_field_active)
                && "Error: must set either num_field_local or num_field.");
  COMET_INSIST_INTERFACE(&env_, (this->is_inited_num_vector_local ||
                this->is_inited_num_vector_active)
                && "Error: must set either num_vector_local or num_vector.");
}

//-----------------------------------------------------------------------------

#if 0
void set_vectors(GMVectors* vectors, DriverOptions* do_, CEnv* env) {
  COMET_INSIST(vectors && do_ && env);

  if (do_->input_file != NULL) {
    VectorsIO::read(*vectors, do_->input_file, *env);
  } else {
    set_vectors_synthetic(vectors, do_->problem_type, do_->verbosity, env);
  }
}
#endif



void Driver::set_vectors(GMVectors& vectors) {
  if (this->input_file)
    VectorsIO::read(vectors, this->input_file, env_);
  else
    set_vectors_synthetic(&vectors, this->problem_type, this->verbosity, &env_);
}

//-----------------------------------------------------------------------------







//-----------------------------------------------------------------------------
// Print a line of output to summarize result of run.

#if 0
void print_output(bool do_print,
                  Checksum& cksum,
                  CEnv& env,
                  char* output_file_stub,
                  size_t num_written,
                  double vctime,
                  double mctime,
                  double cktime,
                  double intime,
                  double outtime,
                  double tottime) {

  const double ops = env.ops();
  const double ops_gemm = env.ops_gemm();
  const double gemmtime_sum = env.gemmtime_sum();
  const size_t cpu_mem_max = env.cpu_mem_max();
  const size_t gpu_mem_max = env.gpu_mem_max();

  if (!do_print)
    return;

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

  if (gemmtime_sum > 0) {
    printf(" gemmrate/proc %e", ops_gemm / gemmtime_sum);
  }

  printf(" vcmp %e", env.vec_compares());
  printf(" vacmp %e", env.vec_active_compares());
  if (output_file_stub) {
    printf(" vcmpout %e", (double)num_written);
  }

  printf(" cmp %e", env.entry_compares());
  printf(" acmp %e", env.entry_active_compares());

  printf(" ecmp %e", env.metric_compares());
  printf(" eacmp %e", env.metric_active_compares());
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

  printf("\n");
}
#endif







void Driver::print_output(Checksum& cksum) {

  const double ops = env_.ops();
  const double ops_gemm = env_.ops_gemm();
  const double gemmtime_sum = env_.gemmtime_sum();
  const size_t cpu_mem_max = env_.cpu_mem_max();
  const size_t gpu_mem_max = env_.gpu_mem_max();

  const bool do_print = env_.is_proc_active() &&
     env_.proc_num() == 0 && this->verbosity > 0;

  if (!do_print)
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

  if (gemmtime_sum > 0) {
    printf(" gemmrate/proc %e", ops_gemm / gemmtime_sum);
  }

  printf(" vcmp %e", env_.vec_compares());
  printf(" vacmp %e", env_.vec_active_compares());
  if (output_file_stub) {
    printf(" vcmpout %e", (double)num_written);
  }

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
    if (cksum.num() > 0) {
      printf(" fracnonzero %.9f",
        (cksum.num()-cksum.num_zero()) / cksum.num());
    }
  }

  if (env_.is_shrink())
    printf(" shrink %e", env_.shrink_achieved());

  printf(" vctime %.6f", this->vctime);
  printf(" mctime %.6f", this->mctime);
  if (cksum.computing_checksum()) {
    printf(" cktime %.6f", this->cktime);
  }
  printf(" intime %.6f", this->intime);
  printf(" outtime %.6f", this->outtime);

  printf(" cpumem %e", (double)cpu_mem_max);
  printf(" gpumem %e", (double)gpu_mem_max);

  printf(" tottime %.6f", this->tottime);

  printf(" prec %s", env_.is_double_prec() ? "double" : "single");

  printf(" build %s", BuildHas::DEBUG ? "debug" : "release");

  if (env_.tc() != env_.tc_eff()) {
    printf(" tc_eff %i", env_.tc_eff());
  }

  if (env_.is_shrink()) {
    printf(" is_shrink %s", "yes");
  }

  printf("\n");
}



void Driver::print_output(Checksum& cksum, CEnv& env) {
  Driver driver(env);
  driver.print_output(cksum);
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

  //double total_time_beg = env->synced_time();
  Driver::Timer timer_total(*env);
  timer_total.start();

  // Parse remaining unprocessed arguments.

  //DriverOptions do_(*env);
  //do_.finish_parsing(argc, argv, *env);

  Driver driver(*env);
  driver.finish_parsing(argc, argv);

  // Set up parallel deomp for vectors, metrics.

  Driver::Timer timer(*env);
  //double time_beg = env->synced_time();
  timer.start();
  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm,
    driver.is_inited_num_field_local,
    driver.is_inited_num_vector_local,
    driver.is_inited_num_field_local ? driver.num_field_local
                                    : driver.num_field_active,
    driver.is_inited_num_vector_local ? driver.num_vector_local
                                     : driver.num_vector_active,
    env->data_type_vectors(), env);
  //driver.vctime += env->synced_time() - time_beg;
  driver.vctime += timer.elapsed();

//TODO: possibly replace this with stuff from dm
  if (driver.is_inited_num_vector_local) {
    driver.num_vector = driver.num_vector_local *
      (size_t)env->num_proc_vector();
    driver.num_vector_active = driver.num_vector;
  } else {
    // Pad up so that every proc has same number of vectors.
    driver.num_vector_local = gm_nvl_size_required(
      utils::ceil(driver.num_vector_active, (size_t)env->num_proc_vector()), *env);
    driver.num_vector = driver.num_vector_local *
      (size_t)env->num_proc_vector();
  }

  if (driver.is_inited_num_field_local) {
    driver.num_field = driver.num_field_local * (size_t) env->num_proc_field();
    driver.num_field_active = driver.num_field;
  } else {
    // Pad up so that every proc has same number of fields.
    driver.num_field_local = utils::ceil(
        driver.num_field_active, (size_t)env->num_proc_field());
    driver.num_field = driver.num_field_local * (size_t) env->num_proc_field();
  }

  const bool do_print = env->is_proc_active() &&
     env->proc_num() == 0 && driver.verbosity > 0;

  // Allocate vectors.

  //time_beg = env->synced_time();
  timer.start();
  GMVectors vectors_value = GMVectors_null(), *vectors = &vectors_value;
  GMVectors_create(vectors, env->data_type_vectors(), dm, env);
  //driver.vctime += env->synced_time() - time_beg;
  driver.vctime += timer.elapsed();

  // Set vectors.

  //time_beg = env->synced_time();
  timer.start();
  //set_vectors(vectors, &do_, env);
  driver.set_vectors(*vectors);
  //driver.intime += env->synced_time() - time_beg;
  driver.intime += timer.elapsed();

  // More initializations.

  comet::Checksum cksum(driver.checksum);
  comet::Checksum cksum_local(driver.checksum);

  {
    // Initialize output.

    //time_beg = env->synced_time();
    timer.start();
    MetricsIO metrics_io(driver.output_file_stub, driver.verbosity, *env);
    //driver.outtime += env->synced_time() - time_beg;
    driver.outtime += timer.elapsed();

    Histograms histograms(driver.histograms_file, *env);
    dm->attach_histograms(&histograms);

  {
    // Initialize metrics mem, compute metrics.

    MetricsMem metrics_mem(env);

    ComputeMetrics compute_metrics(*dm, *env);

  //--------------------
  // Begin loops over phases, stages.
  //--------------------

  for (int phase_num=driver.phase_min; phase_num<=driver.phase_max; ++phase_num) {
      env->phase_num(phase_num);

    for (int stage_num=driver.stage_min; stage_num<=driver.stage_max; ++stage_num) {
      env->stage_num(stage_num);

      // Set up metrics object to capture results.

      //time_beg = env->synced_time();
      timer.start();
      GMMetrics metrics_value = GMMetrics_null(), *metrics = &metrics_value;
      GMMetrics_create(metrics, env->data_type_metrics(), dm,
                       &metrics_mem, env);
      //driver.mctime += env->synced_time() - time_beg;
      driver.mctime += timer.elapsed();

      // Calculate metrics.

      compute_metrics.compute(*metrics, *vectors);
      driver.num_metric_items_local_computed += metrics->num_metric_items_local_computed;
      driver.num_metrics_active_local += metrics->num_metrics_active_local;

      // Output results.

      //time_beg = env->synced_time();
      timer.start();
      metrics_io.write(*metrics);
      if (BuildHas::DEBUG)
        metrics_io.check_file(*metrics);
      //driver.outtime += env->synced_time() - time_beg;
      driver.outtime += timer.elapsed();

      // Check correctness.

      timer.start();
      if (driver.checksum) {
        //time_beg = env->synced_time();
        check_metrics(metrics, driver, env);
        //driver.cktime += env->synced_time() - time_beg;
      }
      driver.cktime += timer.elapsed();

      // Compute checksum.

      timer.start();
      if (driver.checksum) {
        //time_beg = env->synced_time();
        comet::Checksum::compute(cksum, cksum_local, *metrics, *env);
        //driver.cktime += env->synced_time() - time_beg;
      }
      driver.cktime += timer.elapsed();

      //time_beg = env->synced_time();
      timer.start();
      GMMetrics_destroy(metrics, env);
      //driver.mctime += env->synced_time() - time_beg;
      driver.mctime += timer.elapsed();

      if (do_print) {
        if (env->num_phase() > 1 && env->num_stage() > 1) {
          printf("Completed phase %i stage %i\n",
                 env->phase_num(), env->stage_num());
        } else if (env->num_phase() > 1) {
          printf("Completed phase %i\n",
                 env->phase_num());
        } else if (env->num_stage() > 1) {
          printf("Completed stage %i\n",
                 env->stage_num());
        }
      } // do_print

    } // for stage_num

  } // for phase_num

  //--------------------
  // End loops over phases, stages.
  //--------------------

    // Finalize metrics mem.

    //time_beg = env->synced_time();
    timer.start();

  }
    //driver.mctime += env->synced_time() - time_beg;
    driver.mctime += timer.elapsed();

    // Finalize output.

    driver.num_local_written += metrics_io.num_written();
    //time_beg = env->synced_time();
    timer.start();

    histograms.finalize();
    histograms.output();
    if (driver.phase_min==0 && driver.phase_max==env->num_phase() - 1 &&
        driver.stage_min==0 && driver.stage_max==env->num_stage() - 1) {
      histograms.check(driver.num_vector_active);
    }
  }
  //driver.outtime += env->synced_time() - time_beg;
  driver.outtime += timer.elapsed();

  // Deallocate vectors.

  //time_beg = env->synced_time();
  timer.start();
  GMVectors_destroy(vectors, env);
  GMDecompMgr_destroy(dm, env);
  //driver.vctime += env->synced_time() - time_beg;
  driver.vctime += timer.elapsed();

  // Perform some checks.

  COMET_INSIST(env->cpu_mem_local() == 0);
  COMET_INSIST(env->gpu_mem_local() == 0);

  if (env->is_proc_active()) {

    COMET_MPI_SAFE_CALL(MPI_Allreduce(&driver.num_metric_items_local_computed,
      &driver.num_metric_items_computed, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
      env->comm_repl_vector()));

    COMET_MPI_SAFE_CALL(MPI_Allreduce(&driver.num_metrics_active_local,
      &driver.num_metrics_active, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
      env->comm_repl_vector()));

    COMET_MPI_SAFE_CALL(MPI_Allreduce(&driver.num_local_written, &driver.num_written, 1,
      MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));

    const size_t num_metrics_expected =
      utils::nchoosek(dm->num_vector, env->num_way());

    const size_t num_metrics_active_expected =
      utils::nchoosek(dm->num_vector_active, env->num_way());

    if (env->all2all() &&
        driver.phase_min==0 && driver.phase_max==env->num_phase() - 1 &&
        driver.stage_min==0 && driver.stage_max==env->num_stage() - 1) {

      COMET_INSIST(driver.num_metrics_active == num_metrics_active_expected);

      if (!env->is_shrink()) {
        COMET_INSIST(driver.num_metric_items_computed ==
                     env->num_metric_items_per_metric() * num_metrics_expected);
      }
    }

#if 0
    if (env->num_way() == NumWay::_2 && env->all2all() && !env->is_shrink() &&
        driver.phase_min==0 && driver.phase_max==env->num_phase() - 1) {

      //const size_t num_metrics_expected = ((driver.num_vector) * (size_t)
      //                                  (driver.num_vector - 1)) / 2;

      COMET_INSIST(driver.num_metric_items_computed ==
                   env->num_metric_items_per_metric() * num_metrics_expected);

      COMET_INSIST(num_metrics_active == num_metrics_active_expected);
    }

    if (env->num_way() == NumWay::_3 && env->all2all() && !env->is_shrink() &&
        driver.phase_min==0 && driver.phase_max==env->num_phase() - 1 &&
        driver.stage_min==0 && driver.stage_max==env->num_stage() - 1) {

      //const size_t num_metrics_expected = ((driver.num_vector) * (size_t)
      //                                  (driver.num_vector - 1) * (size_t)
      //                                  (driver.num_vector - 2)) / 6;

      COMET_INSIST(driver.num_metric_items_computed ==
                   env->num_metric_items_per_metric() * num_metrics_expected);

      COMET_INSIST(num_metrics_active == num_metrics_active_expected);
    }
#endif

  } // if is_proc_active

  //driver.tottime = env->synced_time() - total_time_beg;
  driver.tottime = timer_total.elapsed();

  // Output run information.

  //print_output(do_print, cksum, *env, driver.output_file_stub, driver.num_written,
  //  driver.vctime, driver.mctime, driver.cktime, driver.intime, driver.outtime, driver.tottime);
  driver.print_output(cksum);

  // Output a local checksum, for testing purposes.

  if (false) {
    // One more sync before checking num_correct, to allow flush of output.
    env->synced_time();
    if (driver.checksum && env->is_proc_active() && driver.verbosity > 0) {
      printf("local checksum: ");
      cksum_local.print(*env);
      printf("\n");
    }
  }
  env->synced_time();
  fflush(NULL);
  env->synced_time();

  // Validation: check for any wrong answers.

  if (driver.num_incorrect) {
    const size_t hnlen = 256;
    char hn[hnlen];
    gethostname(hn, hnlen);
    int rank = 0;
    COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

    fprintf(stderr, "Error: incorrect results found.  num_incorrect  %zu  "
           "max_incorrect_diff  %e  hostname  %s  rank  %i\n",
           driver.num_incorrect, driver.max_incorrect_diff, hn, rank);
  }
  env->synced_time();
  fflush(NULL);
  env->synced_time();

  COMET_INSIST(driver.num_incorrect == 0);

  // Finalize.

  if (env_local) {
    delete env_local;
  }

  cksum_result.copy(cksum);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
