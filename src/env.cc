/*---------------------------------------------------------------------------*/
/*!
 * \file   env.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Environment settings and general utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================


=============================================================================*/

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include <signal.h>

#ifdef TESTING
#include "gtest/gtest.h"
#endif

#include "mpi.h"
#include "cuda.h"

#include "env.hh"

/*---------------------------------------------------------------------------*/

static void gm_test_wrapper() {
#ifdef TESTING
  ASSERT_TRUE(0);
#endif
}

/*---------------------------------------------------------------------------*/

/*===========================================================================*/
/*---Assertions---*/

void gm_assert(const char* condition_string, const char* file, int line) {
  fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Assertion error",
          condition_string, file, line);
  gm_test_wrapper();
#ifdef GM_ASSERTIONS_ON
  raise(SIGUSR1);
#else
  exit(EXIT_FAILURE);
#endif
}

/*---------------------------------------------------------------------------*/

void gm_insist(const void* const env,
               const char* condition_string,
               const char* file,
               int line) {
  if (GMEnv_proc_num((const GMEnv* const)env) == 0) {
    fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Interface error",
            condition_string, file, line);
  }
  gm_test_wrapper();
  exit(EXIT_FAILURE);
}

/*===========================================================================*/
/*---Types---*/

GMMirroredPointer GMMirroredPointer_null(void) {
  GMMirroredPointer p;
  p.h = NULL;
  p.d = NULL;
  p.size = 0;
  p.dim0 = 0;
  p.dim1 = 0;
  return p;
}

/*===========================================================================*/
/*---Null object---*/

GMEnv GMEnv_null() {
  GMEnv result;
  memset((void*)&result, 0, sizeof(GMEnv));
  return result;
}

/*===========================================================================*/
/*---Utility to parse a string to construct arguments---*/

void gm_create_args(char* argstring, int* argc, char** argv) {
  size_t len = strlen(argstring);

  argv[0] = &argstring[0];
  *argc = 1;
  bool is_delim_prev = true;
  int i = 0;
  for (i = 0; i < (int)len; ++i) {
    const bool is_delim = argstring[i] == ' ' || argstring[i] == '\t';
    if (is_delim) {
      argstring[i] = 0;
    }
    if (is_delim_prev && !is_delim) {
      argv[*argc] = &(argstring[i]);
      (*argc)++;
    }
    is_delim_prev = is_delim;
  }
}

/*===========================================================================*/
/*---Initialize environment---*/

void GMEnv_create_impl_(GMEnv* const env, MPI_Comm comm, int argc,
                        char** argv, const char* const description,
                        bool make_comms, int num_proc, int proc_num) {
  GMAssertAlways(env != NULL);

  *env = GMEnv_null();

  /*---Set default values---*/
  env->metric_type_ = GM_METRIC_TYPE_CZEK;
  env->num_way_ = GM_NUM_WAY_2;
  env->all2all_ = false;
  env->are_cuda_streams_initialized_ = false;
  env->are_mpi_comms_initialized_ = false;
  GMEnv_set_compute_method(env, GM_COMPUTE_METHOD_GPU);
  env->num_stage = 1;
  env->stage_num = 0;
  env->num_phase = 1;
  env->phase_num = 0;
  env->sparse = false;
  env->ccc_param_ = ((double) 2) / ((double) 3);

  env->time = 0;
  env->compares = 0;
  env->ops_local = 0;
  env->ops = 0;
  env->cpu_mem = 0;
  env->cpu_mem_max = 0;
  env->gpu_mem = 0;
  env->gpu_mem_max = 0;
  env->description = description;

  env->mpi_comm_base_ = comm;
  env->make_comms_ = make_comms;

  if (env->make_comms_) {
    int mpi_code = MPI_Comm_size(env->mpi_comm_base_, &env->num_proc_base_);
    GMAssertAlways(mpi_code == MPI_SUCCESS);
    mpi_code = MPI_Comm_rank(env->mpi_comm_base_, &env->proc_num_base_);
    GMAssertAlways(mpi_code == MPI_SUCCESS);
  } else {
    env->num_proc_base_ = num_proc;
    env->proc_num_base_ = proc_num;
  }

  GMEnv_set_num_proc(env, env->num_proc_base_, 1, 1);

  /*---Modify based on user options---*/
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--metric_type") == 0) {
      /*--------------------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for metric_type." : 0);
      if (strcmp(argv[i], "czekanowski") == 0) {
        env->metric_type_ = GM_METRIC_TYPE_CZEK;
      } else if (strcmp(argv[i], "ccc") == 0) {
        env->metric_type_ = GM_METRIC_TYPE_CCC;
      } else {
        GMInsist(env, false ? "Invalid setting for metric_type." : 0);
      }
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_way") == 0) {
      /*--------------------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_way." : 0);
      errno = 0;
      const long num_way = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno && (num_way == GM_NUM_WAY_2 ||
                                   num_way == GM_NUM_WAY_3)
                               && "Invalid setting for num_way.");
      env->num_way_ = num_way;
      GMEnv_set_num_proc(env, env->num_proc_vector_i_, env->num_proc_repl_,
                       env->num_proc_field_);
      /*--------------------*/
    } else if (strcmp(argv[i], "--all2all") == 0) {
      /*--------------------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for all2all." : 0);
      if (strcmp(argv[i], "yes") == 0) {
        env->all2all_ = true;
      } else if (strcmp(argv[i], "no") == 0) {
        env->all2all_ = false;
      } else {
        GMInsist(env, false ? "Invalid setting for all2all." : 0);
      }
      /*--------------------*/
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      /*--------------------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for compute_method." : 0);
      if (strcmp(argv[i], "CPU") == 0) {
        GMEnv_set_compute_method(env, GM_COMPUTE_METHOD_CPU);
      } else if (strcmp(argv[i], "GPU") == 0) {
        GMEnv_set_compute_method(env, GM_COMPUTE_METHOD_GPU);
      } else if (strcmp(argv[i], "REF") == 0) {
        GMEnv_set_compute_method(env, GM_COMPUTE_METHOD_REF);
      } else {
        GMInsist(env, false ? "Invalid setting for compute_method." : 0);
      }
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      /*--------------------*/
      ++i;
      errno = 0;
      GMInsist(env, i < argc && "Missing value for num_proc_vector.");
      long num_proc_vector_i = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno
                    && (long)(int)num_proc_vector_i == num_proc_vector_i
                    && "Invalid setting for num_proc_vector.");
      GMEnv_set_num_proc(env, num_proc_vector_i, env->num_proc_repl_,
                         env->num_proc_field_);
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      /*--------------------*/
      ++i;
      errno = 0;
      GMInsist(env, i < argc && "Missing value for num_proc_field.");
      long num_proc_field = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno
                    && (long)(int)num_proc_field == num_proc_field
                    && "Invalid setting for num_proc_field.");
      GMEnv_set_num_proc(env, env->num_proc_vector_i_, env->num_proc_repl_,
                         num_proc_field);
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      /*--------------------*/
      ++i;
      errno = 0;
      GMInsist(env, i < argc && "Missing value for num_proc_repl.");
      long num_proc_repl = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno
                    && (long)(int)num_proc_repl == num_proc_repl
                    && "Invalid setting for num_proc_repl.");
      GMEnv_set_num_proc(env, env->num_proc_vector_i_, num_proc_repl,
                         env->num_proc_field_);
      /*--------------------*/
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      /*--------------------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for ccc_param." : 0);
      errno = 0;
      const double ccc_param = strtod(argv[i], NULL);
      GMInsist(env, 0 == errno && ccc_param >= 0
                               && "Invalid setting for ccc_param.");
      env->ccc_param_ = ccc_param;
      /*--------------------*/
    } else if (strcmp(argv[i], "--sparse") == 0) {
      /*--------------------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for sparse." : 0);
      if (strcmp(argv[i], "yes") == 0) {
        env->sparse = true;
      } else if (strcmp(argv[i], "no") == 0) {
        env->sparse = false;
      } else {
        GMInsist(env, false ? "Invalid setting for sparse." : 0);
      }
      /*--------------------*/
    } /*---if/else---*/
  }   /*---for i---*/

  /*---Helper variables---*/
  env->do_reduce = env->num_proc_field_ > 1;
  env->need_2way = env->metric_type_ == GM_METRIC_TYPE_CZEK;
}

/*---------------------------------------------------------------------------*/

void GMEnv_create(GMEnv* const env, MPI_Comm comm, int argc, char** argv,
                  const char* const description) {
  GMAssertAlways(env != NULL);

  GMEnv_create_impl_(env, comm, argc, argv, description, true, 0, 0);
}

/*---------------------------------------------------------------------------*/

void GMEnv_create(GMEnv* const env, MPI_Comm comm, const char* const options,
                  const char* const description) {
  GMAssertAlways(env != NULL);

  /*---Convert options string to args---*/

  size_t len = strlen(options);
  char argstring[len+1];
  char* argv[len+1];
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);

  GMEnv_create_impl_(env, comm, argc, argv, description, true, 0, 0);
}

/*---------------------------------------------------------------------------*/

void GMEnv_create_no_comms(GMEnv* const env, int argc, char** argv,
                           const char* const description,
                           int num_proc, int proc_num) {
  GMAssertAlways(env != NULL);

  GMEnv_create_impl_(env, MPI_COMM_WORLD, argc, argv, description,
                     false, num_proc, proc_num);
}

/*---------------------------------------------------------------------------*/

void GMEnv_create_no_comms(GMEnv* const env, const char* const options,
                           const char* const description,
                           int num_proc, int proc_num) {
  GMAssertAlways(env != NULL);

  /*---Convert options string to args---*/

  size_t len = strlen(options);
  char argstring[len+1];
  char* argv[len+1];
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);

  GMEnv_create_impl_(env, MPI_COMM_WORLD, argc, argv, description,
                     false, num_proc, proc_num);
}

/*===========================================================================*/
/*---Manage cuda streams---*/

void GMEnv_initialize_streams(GMEnv* const env) {
  GMAssertAlways(env != NULL);

  /*---NOTE: this is used for lazy initialization---*/

  if (env->are_cuda_streams_initialized_) {
    return;
  }

  if (env->compute_method_ != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  cudaStreamCreate(&env->stream_compute_);
  GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));

  cudaStreamCreate(&env->stream_togpu_);
  GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));

  cudaStreamCreate(&env->stream_fromgpu_);
  GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));

  env->are_cuda_streams_initialized_ = true;
}

/*---------------------------------------------------------------------------*/

void GMEnv_terminate_streams(GMEnv* const env) {
  GMAssertAlways(env != NULL);

  if (!env->are_cuda_streams_initialized_) {
    return;
  }

  cudaStreamDestroy(env->stream_compute_);
  GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));

  cudaStreamDestroy(env->stream_togpu_);
  GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));

  cudaStreamDestroy(env->stream_fromgpu_);
  GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));

  env->are_cuda_streams_initialized_ = false;
}

/*===========================================================================*/
/*---Manage MPI comms---*/

void GMEnv_initialize_comms(GMEnv* const env) {
  GMAssertAlways(env != NULL);

  if (env->are_mpi_comms_initialized_) {
    return;
  }

  if (!env->make_comms_) {
    return;
  }

  int mpi_code = MPI_Comm_split(env->mpi_comm_base_, env->is_proc_active_,
                            env->proc_num_, &env->mpi_comm_);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  mpi_code = MPI_Comm_split(env->mpi_comm_base_,
      env->is_proc_active_ ? env->proc_num_field_ : env->num_proc_,
      env->proc_num_, &env->mpi_comm_vector_);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  mpi_code = MPI_Comm_split(env->mpi_comm_base_,
      env->is_proc_active_ ? env->proc_num_vector_ : env->num_proc_,
      env->proc_num_, &env->mpi_comm_field_);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  env->are_mpi_comms_initialized_ = true;
}

/*---------------------------------------------------------------------------*/

void GMEnv_terminate_comms(GMEnv* const env) {
  GMAssertAlways(env != NULL);

  if (!env->are_mpi_comms_initialized_) {
    return;
  }

  /*---Destroy any nontrivial communicators---*/

  int mpi_code = MPI_Comm_free(&(env->mpi_comm_));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  mpi_code = MPI_Comm_free(&(env->mpi_comm_vector_));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  mpi_code = MPI_Comm_free(&(env->mpi_comm_field_));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  env->are_mpi_comms_initialized_ = false;
}

/*===========================================================================*/
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* const env) {
  GMAssertAlways(env != NULL);

  GMEnv_terminate_comms(env);

  GMEnv_terminate_streams(env);

  *env = GMEnv_null();
}

/*===========================================================================*/
/*---Accessors---*/

void GMEnv_set_compute_method(GMEnv* const env, int compute_method) {
  GMAssertAlways(env != NULL);
  GMAssertAlways(compute_method >= 0);
  GMAssertAlways(compute_method < GM_NUM_COMPUTE_METHOD);

  env->compute_method_ = compute_method;
}

/*---------------------------------------------------------------------------*/

int GMEnv_data_type_vectors(const GMEnv* const env) {
  GMAssertAlways(env != NULL);

  switch (env->metric_type_) {
    case GM_METRIC_TYPE_CZEK:
      return GM_DATA_TYPE_FLOAT;
    case GM_METRIC_TYPE_CCC:
      return GM_DATA_TYPE_BITS2;
  }
  GMAssertAlways(false ? "Invalid metric type." : 0);
  return 0;
}

/*---------------------------------------------------------------------------*/

int GMEnv_data_type_metrics(const GMEnv* const env) {
  GMAssertAlways(env != NULL);

  switch (env->metric_type_) {
    case GM_METRIC_TYPE_CZEK:
      return GM_DATA_TYPE_FLOAT;
    case GM_METRIC_TYPE_CCC:
      return env->num_way_ == GM_NUM_WAY_2 ? GM_DATA_TYPE_TALLY2X2
                                           : GM_DATA_TYPE_TALLY4X2;
  }
  GMAssertAlways(false ? "Invalid metric type." : 0);
  return 0;
}

/*---------------------------------------------------------------------------*/

void GMEnv_set_num_proc(GMEnv* const env, int num_proc_vector_i,
                      int num_proc_repl, int num_proc_field) {
  GMAssertAlways(env != NULL);
  GMAssertAlways(num_proc_vector_i > 0);
  GMAssertAlways(num_proc_repl > 0);
  GMAssertAlways(num_proc_field > 0);

  GMAssertAlways(env->num_proc_base_ != 0);
  //GMAssertAlways(env->proc_num_base_ is initialized);

  /*---Set proc counts---*/

  env->num_proc_vector_i_ = num_proc_vector_i;
  env->num_proc_repl_ = num_proc_repl;
  env->num_proc_field_ = num_proc_field;

  env->num_proc_vector_total_ = env->num_proc_vector_i_ * env->num_proc_repl_;

  env->num_proc_ = env->num_proc_vector_total_ * num_proc_field;
  GMAssertAlways(env->num_proc_ <= env->num_proc_base_);

  /*---Set proc nums---*/

  env->proc_num_ = env->proc_num_base_;

  env->is_proc_active_ = env->proc_num_ < env->num_proc_;

  int itmp = env->proc_num_;

  env->proc_num_repl_ = itmp % env->num_proc_repl_;
  itmp /= env->num_proc_repl_;

  env->proc_num_vector_i_ = itmp % env->num_proc_vector_i_;
  itmp /= env->num_proc_vector_i_;

  env->proc_num_field_ = itmp % env->num_proc_field_;
  env->proc_num_vector_ = env->proc_num_ % env->num_proc_vector_total_;

  /*---Destroy old communicators if necessary---*/

  GMEnv_terminate_comms(env);

  /*---Make new communicators---*/

  GMEnv_initialize_comms(env);
}

/*---------------------------------------------------------------------------*/

cudaStream_t GMEnv_stream_compute(GMEnv* const env) {
  GMAssertAlways(env != NULL);
  // GMAssertAlways(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_compute_;
}

/*---------------------------------------------------------------------------*/

cudaStream_t GMEnv_stream_togpu(GMEnv* const env) {
  GMAssertAlways(env != NULL);
  // GMAssertAlways(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_togpu_;
}

/*---------------------------------------------------------------------------*/

cudaStream_t GMEnv_stream_fromgpu(GMEnv* const env) {
  GMAssertAlways(env != NULL);
  // GMAssertAlways(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_fromgpu_;
}

/*===========================================================================*/
/*---Timer functions---*/

double GMEnv_get_time(const GMEnv* const env) {
  GMAssertAlways(env);

  struct timeval tv;
  gettimeofday(&tv, NULL);
  double result = ((double)tv.tv_sec + (double)tv.tv_usec * 1.e-6);

  return result;
}

/*---------------------------------------------------------------------------*/

double GMEnv_get_synced_time(const GMEnv* const env) {
  GMAssertAlways(env != NULL);

  /*
  cudaThreadSynchronize();
  */

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    cudaDeviceSynchronize();
    GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));
  }

  const int mpi_code = MPI_Barrier(GMEnv_mpi_comm(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);
  return GMEnv_get_time(env);
}

/*===========================================================================*/
/*---Misc.---*/

bool GMEnv_cuda_last_call_succeeded(const GMEnv* const env) {
  GMAssertAlways(env);

  bool result = true;

  /*---NOTE: this read of the last error is a destructive read---*/
  cudaError_t error = cudaGetLastError();

  if (error != cudaSuccess) {
    result = false;
    printf("CUDA error detected: %s\n", cudaGetErrorString(error));
  }

  return result;
}

/*---------------------------------------------------------------------------*/

void* gm_malloc(size_t n, GMEnv* env) {
  void* p = malloc(n);
  GMAssertAlways(p != NULL);
  env->cpu_mem += n;
  env->cpu_mem_max = gm_max_i8(env->cpu_mem_max, env->cpu_mem);
  return p;
}

/*---------------------------------------------------------------------------*/

void gm_free(void* p, size_t n, GMEnv* env) {
  free(p);
  env->cpu_mem -= n;
}

/*---------------------------------------------------------------------------*/

GMFloat* GMFloat_malloc(size_t n, GMEnv* env) {
  GMFloat* p = (GMFloat*)gm_malloc(n * sizeof(GMFloat), env);
  GMAssertAlways(p != NULL);
  GMFloat_fill_nan(p, n);
  return p;
}

/*---------------------------------------------------------------------------*/

void GMFloat_free(GMFloat* p, size_t n, GMEnv* env) {
  gm_free(p, n * sizeof(GMFloat), env);
}

/*---------------------------------------------------------------------------*/

void GMFloat_fill_nan(GMFloat* const a, size_t n) {
  GMAssertAlways(a != NULL);
  GMAssertAlways(n+1 >= 1);
#ifdef GM_ASSERTIONS_ON
  GMFloat value = sqrt(-1);
  size_t i = 0;
  for (i=0; i<n; ++i) {
    a[i] = value;
  }
#endif
}

/*---------------------------------------------------------------------------*/

void GMFloat_check(GMFloat* const a, size_t n) {
  GMAssertAlways(a != NULL);
  GMAssertAlways(n+1 >= 1);
#ifdef GM_ASSERTIONS_ON
  bool no_nans_found = true;
  size_t i = 0;
  for (i=0; i<n; ++i) {
    if (a[i] != a[i]) {
      no_nans_found = false;
    }
  }
  GMAssertAlways(no_nans_found);
#endif
}

/*---------------------------------------------------------------------------*/

int gm_mpi_type(const GMEnv* const env) {
  GMAssertAlways(env != NULL);

  /* clang-format off */
  const int mpi_type = GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK ?
                         GM_MPI_FLOAT :
                       GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC ?
                         MPI_DOUBLE_COMPLEX :
                       0;
  /* clang-format on */

  return mpi_type;
}

/*===========================================================================*/

//#ifdef __cplusplus
//extern "C" {
//#endif
//
//#ifdef __cplusplus
//} /*---extern "C"---*/
//#endif

/*---------------------------------------------------------------------------*/
