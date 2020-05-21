//-----------------------------------------------------------------------------
/*!
 * \file   driver.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

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
    } else if (strcmp(argv[i], "--tc") == 0) {
      ++i; // optionally processed by caller.
    } else if (strcmp(argv[i], "--num_tc_steps") == 0) {
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
                  Checksum& cksum,
                  CEnv& env,
                  char* metrics_file_path_stub,
                  size_t num_written,
                  double vctime,
                  double mctime,
                  double cktime,
                  double intime,
                  double outtime,
                  double tottime) {

  const double ops = env.ops();
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

  printf(" vcmp %e", env.veccompares());
  if (metrics_file_path_stub) {
    printf(" vcmpout %e", (double)num_written);
  }

  printf(" cmp %e", env.compares());
  printf(" ecmp %e", env.eltcompares());
  if (env.ctime() > 0) {
    printf(" ecmp_rate %e", env.eltcompares() / env.ctime());
    printf(" ecmp_rate/proc %e", env.eltcompares() /
      (env.ctime() * env.num_proc()) );
  }
  if (cksum.computing_checksum()) {
    printf(" me %.0f", cksum.num());
    printf(" mezero %.0f", cksum.num_zero());
    if (cksum.num() > 0) {
      printf(" fracnonzero %.9f",
        (cksum.num()-cksum.num_zero()) / cksum.num());
    }
  }

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

  if (env.tc() != env.tc_eff()) {
    printf(" tc_eff %i", env.tc_eff());
  }

  printf("\n");
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
  double cktime = 0;

  size_t num_elts_local_computed = 0;
  size_t num_local_written = 0;

  // Open output files.

  {
  time_beg = env->synced_time();
  MetricsIO metrics_io(do_.metrics_file_path_stub, do_.verbosity, *env);
  time_end = env->synced_time();
  outtime += time_end - time_beg;

  {
  MetricsMem metrics_mem(env);

  ComputeMetrics compute_metrics(*dm, *env);

  // Loops over phases, stages.

  for (int phase_num=do_.phase_min_0based; phase_num<=do_.phase_max_0based;
       ++phase_num) {
      env->phase_num(phase_num);

    for (int stage_num=do_.stage_min_0based; stage_num<=do_.stage_max_0based;
         ++stage_num) {
      env->stage_num(stage_num);

      // Set up metrics container for results.

      time_beg = env->synced_time();
      GMMetrics metrics_value = GMMetrics_null(), *metrics = &metrics_value;
      GMMetrics_create(metrics, env->data_type_metrics(), dm,
                       &metrics_mem, env);
      time_end = env->synced_time();
      mctime += time_end - time_beg;

      // Calculate metrics.

      compute_metrics.compute(*metrics, *vectors);

      num_elts_local_computed += metrics->num_elts_local_computed;

      // Output results.

      time_beg = env->synced_time();
      metrics_io.write(*metrics);
      if (BuildHas::DEBUG) {
        metrics_io.check_file(*metrics);
      }
      time_end = env->synced_time();
      outtime += time_end - time_beg;

      // Check correctness.

      if (do_.checksum) {
        time_beg = env->synced_time();
        check_metrics(metrics, &do_, env);
        time_end = env->synced_time();
        cktime += time_end - time_beg;
      }

      // Compute checksum.

      if (do_.checksum) {
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
          printf("Completed phase %i stage %i\n",
                 env->phase_num(), env->stage_num());
        } else if (env->num_phase() > 1) {
          printf("Completed phase %i\n",
                 env->phase_num());
        } else if (env->num_stage() > 1) {
          printf("Completed stage %i\n",
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
    size_t num_elts_computed = 0;
    COMET_MPI_SAFE_CALL(MPI_Allreduce(&num_elts_local_computed,
      &num_elts_computed, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
      env->comm_repl_vector()));

    COMET_MPI_SAFE_CALL(MPI_Allreduce(&num_local_written, &num_written, 1,
      MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));

    if (env->num_way() == NUM_WAY::_2 && env->all2all() &&
        do_.phase_min_0based==0 && do_.phase_max_0based==env->num_phase() - 1) {
      COMET_INSIST(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) / 2);
    }

    if (env->num_way() == NUM_WAY::_3 && env->all2all() &&
        do_.phase_min_0based==0 && do_.phase_max_0based==env->num_phase() - 1 &&
        do_.stage_min_0based==0 && do_.stage_max_0based==env->num_stage() - 1) {
      COMET_INSIST(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) * (size_t)
                                          (do_.num_vector - 2) / 6);
    }
  }

  double total_time_end = env->synced_time();
  double tottime = total_time_end - total_time_beg;

  // Output run information.

  print_output(do_print, cksum, *env, do_.metrics_file_path_stub, num_written,
    vctime, mctime, cktime, intime, outtime, tottime);
    
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

  // Finalize.

  if (env_local) {
    delete env_local;
  }

  cksum_result.copy(cksum);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
