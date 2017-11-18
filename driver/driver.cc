//-----------------------------------------------------------------------------
/*!
 * \file   driver.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdio.h"
#include "stdlib.h"
#include "stddef.h"
#include "string.h"
#include "float.h"
#include "errno.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "checksums.hh"
#include "compute_metrics.hh"
#include "driver.hh"
#include "test_problems.hh"
#include "input_output.hh"

//=============================================================================
/*---Parse remaining unprocessed arguments---*/

void finish_parsing(int argc, char** argv, DriverOptions* do_, GMEnv* env) {
  errno = 0;
  int i = 0;
  for (i = 1; i < argc; ++i) {
    /*----------*/
    if (strcmp(argv[i], "--num_field") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_field.");
      const long num_field = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && num_field >= 0
                    && "Invalid setting for num_field.");
      do_->num_field_active = num_field;
      do_->num_field_active_initialized = true;
      do_->num_field_local_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_field_local") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_field_local.");
      const long num_field_local = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno && num_field_local >= 0 &&
                    (long)(int)num_field_local == num_field_local &&
                    "Invalid setting for num_field_local.");
      do_->num_field_local = num_field_local;
      do_->num_field_local_initialized = true;
      do_->num_field_active_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_vector.");
      const long num_vector = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && num_vector >= 0
                    && "Invalid setting for num_vector.");
      do_->num_vector_active = num_vector;
      do_->num_vector_active_initialized = true;
      do_->num_vector_local_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_vector_local.");
      const long num_vector_local = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno && num_vector_local >= 0 &&
                    (long)(int)num_vector_local == num_vector_local &&
                    "Invalid setting for num_vector_local.");
      do_->num_vector_local = num_vector_local;
      do_->num_vector_local_initialized = true;
      do_->num_vector_active_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--verbosity") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for verbosity.");
      const long verbosity = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno && verbosity >= 0 &&
                    "Invalid setting for verbosity.");
      do_->verbosity = verbosity;
      /*--------------------*/
    } else if (strcmp(argv[i], "--checksum") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for checksum.");
      if (strcmp(argv[i], "yes") == 0) {
        do_->checksum = true;
      } else if (strcmp(argv[i], "no") == 0) {
        do_->checksum = false;
      } else {
        GMInsistInterface(env, false && "Invalid setting for checksum.");
      }
    /*----------*/
    } else if (strcmp(argv[i], "--num_stage") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_stage.");
      const long num_stage = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && num_stage >= 1
                    && (long)(int)num_stage == num_stage
                    && "Invalid setting for num_stage.");
      env->num_stage = num_stage;
      do_->stage_min_1based = 1;
      do_->stage_max_1based = env->num_stage;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_min") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for stage_min.");
      const long stage_min_1based = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && stage_min_1based >= 1
                    && (long)(int)stage_min_1based == stage_min_1based
                    && "Invalid setting for stage_min.");
      do_->stage_min_1based = stage_min_1based;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_max") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for stage_max.");
      const long stage_max_1based = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && stage_max_1based <= env->num_stage
                    && (long)(int)stage_max_1based == stage_max_1based
                    && "Invalid setting for stage_max.");
      do_->stage_max_1based = stage_max_1based;
    /*----------*/
    } else if (strcmp(argv[i], "--num_phase") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_phase.");
      const long num_phase = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && num_phase >= 1
                    && (long)(int)num_phase == num_phase
                    && "Invalid setting for num_phase.");
      env->num_phase = num_phase;
      do_->phase_min_1based = 1;
      do_->phase_max_1based = env->num_phase;
    /*----------*/
    } else if (strcmp(argv[i], "--phase_min") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for phase_min.");
      const long phase_min_1based = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && phase_min_1based >= 1
                    && (long)(int)phase_min_1based == phase_min_1based
                    && "Invalid setting for phase_min.");
      do_->phase_min_1based = phase_min_1based;
    /*----------*/
    } else if (strcmp(argv[i], "--phase_max") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for phase_max.");
      const long phase_max_1based = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, errno == 0 && phase_max_1based <= env->num_phase
                    && (long)(int)phase_max_1based == phase_max_1based
                    && "Invalid setting for phase_max.");
      do_->phase_max_1based = phase_max_1based;
    /*----------*/
    } else if (strcmp(argv[i], "--input_file") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for input_file.");
      do_->input_file_path = argv[i];
    /*----------*/
    } else if (strcmp(argv[i], "--output_file_stub") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for output_file_stub.");
      do_->output_file_path_stub = argv[i];
      /*--------------------*/
    } else if (strcmp(argv[i], "--problem_type") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for problem_type.");
      if (strcmp(argv[i], "random") == 0) {
        GMEnv_set_compute_method(env, GM_PROBLEM_TYPE_RANDOM);
      } else if (strcmp(argv[i], "analytic") == 0) {
        GMEnv_set_compute_method(env, GM_PROBLEM_TYPE_ANALYTIC);
      } else {
        GMInsistInterface(env, false && "Invalid setting for problem_type.");
      }
    /*----------*/
    } else if (strcmp(argv[i], "--threshold") == 0) {
    /*----------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for threshold.");
      errno = 0;
      const double threshold = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && "Invalid setting for threshold.");
      do_->threshold = threshold;
     /*----------*/
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--ccc_multiplier") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--sparse") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else {
    /*----------*/
      if (GMEnv_proc_num(env) == 0) {
        fprintf(stderr, "Invalid argument \"%s\".", argv[i]);
      }
      GMInsistInterface(env, false && "Error: argument not recognized.");
    /*----------*/
    } /*---if/else---*/

  } /*---for i---*/

  GMInsistInterface(env, (do_->num_field_local_initialized ||
                do_->num_field_active_initialized)
                && "Error: must set num_field_local or num_field.");
  GMInsistInterface(env, (do_->num_vector_local_initialized ||
                do_->num_vector_active_initialized)
                && "Error: must set num_vector_local or num_vector.");
}

//-----------------------------------------------------------------------------

void set_vectors(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMInsist(vectors && do_ && env);

  if (do_->input_file_path != NULL) {
    set_vectors_from_file(vectors, do_, env);
  } else if (do_->problem_type == GM_PROBLEM_TYPE_RANDOM) {
    set_vectors_random(vectors, do_, env);
  } else {
    set_vectors_analytic(vectors, do_, env);
  }
}
//=============================================================================
/*---Perform a single metrics computation run---*/

GMChecksum perform_run(const char* const options) {
  GMInsist(options);

  /*---Convert options string to args---*/

  size_t len = strlen(options);
  char argstring[len+1];
  char* argv[len+1];
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);

  return perform_run(argc, argv, options);
}

//-----------------------------------------------------------------------------

GMChecksum perform_run(int argc, char** argv, const char* const description) {

  /*---Initialize environment---*/

  GMEnv env_value = GMEnv_null(), *env = &env_value;;

  GMEnv_create(env, MPI_COMM_WORLD, argc, argv, description);

  /*---Parse remaining unprocessed arguments---*/

  DriverOptions do_ = {0};
  do_.num_field_local_initialized = false;
  do_.num_field_active_initialized = false;
  do_.num_vector_local_initialized = false;
  do_.num_vector_active_initialized = false;
  do_.verbosity = 1;
  do_.stage_min_1based = 1;
  do_.stage_max_1based = env->num_stage;
  do_.phase_min_1based = 1;
  do_.phase_max_1based = env->num_phase;
  do_.input_file_path = NULL;
  do_.output_file_path_stub = NULL;
  //do_.problem_type = GM_PROBLEM_TYPE_RANDOM;
  do_.problem_type = GM_PROBLEM_TYPE_ANALYTIC;
  do_.threshold = -1.;
  do_.checksum = true;
  do_.num_incorrect = 0;

  finish_parsing(argc, argv, &do_, env);

  GMDecompMgr dm_value = GMDecompMgr_null(), *dm = &dm_value;
  GMDecompMgr_create(dm,
    do_.num_field_local_initialized,
    do_.num_vector_local_initialized,
    do_.num_field_local_initialized ? do_.num_field_local
                                    : do_.num_field_active,
    do_.num_vector_local_initialized ? do_.num_vector_local
                                     : do_.num_vector_active,
    GMEnv_data_type_vectors(env), env);





  if (do_.num_vector_local_initialized) {
    do_.num_vector = do_.num_vector_local *
      (size_t)GMEnv_num_proc_vector_i(env);
    do_.num_vector_active = do_.num_vector;
  } else {
    /*---Pad up so that every proc has same number of vectors---*/
    do_.num_vector_local = GMVectors_num_local_required(
        do_.num_vector_active, env);
    do_.num_vector = do_.num_vector_local *
      (size_t)GMEnv_num_proc_vector_i(env);
  }

  if (do_.num_field_local_initialized) {
    do_.num_field = do_.num_field_local * (size_t) GMEnv_num_proc_field(env);
    do_.num_field_active = do_.num_field;
  } else {
    /*---Pad up so that every proc has same number of fields---*/
    do_.num_field_local = gm_ceil_i8(
        do_.num_field_active, GMEnv_num_proc_field(env));
    do_.num_field = do_.num_field_local * (size_t) GMEnv_num_proc_field(env);
  }





  /*---Initialize vectors---*/

  double vctime = 0;
  double time_beg = GMEnv_get_synced_time(env);
  GMVectors vectors_value = GMVectors_null(), *vectors = &vectors_value;
  GMVectors_create(vectors, GMEnv_data_type_vectors(env), dm, env);
  double time_end = GMEnv_get_synced_time(env);
  vctime += time_end - time_beg;

  double intime = 0;
  time_beg = GMEnv_get_synced_time(env);
  set_vectors(vectors, &do_, env);
  time_end = GMEnv_get_synced_time(env);
  intime += time_end - time_beg;

  GMChecksum checksum_value = GMChecksum_null(), *checksum = &checksum_value;
  checksum->computing_checksum = do_.checksum;

  double outtime = 0;
  double mctime = 0;
  double cktime = 0;

  /*---Loops over phases, stages---*/

  size_t num_elts_local_computed = 0;

  for (env->stage_num=do_.stage_min_1based-1;
       env->stage_num<=do_.stage_max_1based-1; ++env->stage_num) {

    for (env->phase_num=do_.phase_min_1based-1;
         env->phase_num<=do_.phase_max_1based-1; ++env->phase_num) {

      /*---Set up metrics container for results---*/

      time_beg = GMEnv_get_synced_time(env);
      GMMetrics metrics_value = GMMetrics_null(), *metrics = &metrics_value;
      GMMetrics_create(metrics, GMEnv_data_type_metrics(env), dm, env);
      time_end = GMEnv_get_synced_time(env);
      mctime += time_end - time_beg;

      /*---Calculate metrics---*/

      gm_compute_metrics(metrics, vectors, env);

      num_elts_local_computed += metrics->num_elts_local_computed;

      /*---Output results---*/

      time_beg = GMEnv_get_synced_time(env);
      output_metrics(metrics, &do_, env);
      time_end = GMEnv_get_synced_time(env);
      outtime += time_end - time_beg;

      /*---Check correctness---*/

      check_metrics(metrics, &do_, env);

      /*---Compute checksum---*/

      if (do_.checksum) {
        time_beg = GMEnv_get_synced_time(env);
        GMChecksum_metrics(checksum, metrics, env);
        time_end = GMEnv_get_synced_time(env);
        cktime += time_end - time_beg;
      }

      GMMetrics_destroy(metrics, env);

    }

  } /*---End loops over phases, stages---*/

  GMVectors_destroy(vectors, env);

  /*---Perform some checks---*/

  GMInsist(env->cpu_mem == 0);
  GMInsist(env->gpu_mem == 0);

  if (GMEnv_is_proc_active(env)) {
  int mpi_code = 0;
    size_t num_elts_computed = 0;
    mpi_code = MPI_Allreduce(&num_elts_local_computed, &num_elts_computed, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                             GMEnv_mpi_comm_vector(env));
    GMInsist(mpi_code == MPI_SUCCESS);

    if (GMEnv_num_way(env) == GM_NUM_WAY_2 && GMEnv_all2all(env) &&
        do_.phase_min_1based==1 && do_.phase_max_1based==env->num_phase) {
      GMInsist(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) / 2);
    }

    if (GMEnv_num_way(env) == GM_NUM_WAY_3 && GMEnv_all2all(env) &&
        do_.stage_min_1based==1 && do_.stage_max_1based==env->num_stage) {
      GMInsist(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) * (size_t)
                                          (do_.num_vector - 2) / 6);
    }
  }

  /*---Output run information---*/

  if (GMEnv_is_proc_active(env) && GMEnv_proc_num(env) == 0 &&
      do_.verbosity > 0) {
    //-----
    if (do_.checksum) {
      printf("metrics checksum ");
      int i = 0;
      for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
        printf("%s%li", i == 0 ? "" : "-",
               checksum->data[GM_CHECKSUM_SIZE - 1 - i]);
      }
      if (checksum->is_overflowed) {
        printf("-OVFL");
        printf("-%e", checksum->value_max);
      }
      printf(" ");
    }
    //-----
    printf("time %.6f", env->time);
    //-----
    printf(" ops %e", env->ops);
    if (env->time > 0) {
      printf(" rate %e", env->ops / env->time);
      printf(" rate/proc %e", env->ops / (env->time*GMEnv_num_proc(env)) );
    }
    //-----
    printf(" cmp %e", env->compares);
    if (env->time > 0) {
      printf(" rate %e", env->compares / env->time);
      printf(" rate/proc %e", env->compares / (env->time*GMEnv_num_proc(env)) );
    }
    //-----
    printf(" vctime %.6f", vctime);
    printf(" mctime %.6f", mctime);
    if (do_.checksum) {
      printf(" cktime %.6f", cktime);
    }
    if (NULL != do_.input_file_path) {
      printf(" intime %.6f", intime);
    }
    if (NULL != do_.output_file_path_stub) {
      printf(" outtime %.6f", outtime);
    }
    //-----
    printf(" cpumem %e", (double)env->cpu_mem_max);
    printf(" gpumem %e", (double)env->gpu_mem_max);
    //-----
    printf("\n");
  }

  GMInsist(do_.num_incorrect == 0);

  /*---Finalize---*/

  GMDecompMgr_destroy(dm, env);
  GMEnv_destroy(env);

  return *checksum;
}

//=============================================================================
