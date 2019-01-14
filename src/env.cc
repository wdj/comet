//-----------------------------------------------------------------------------
/*!
 * \file   env.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Environment settings and general utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "sys/time.h"
#include "stdio.h"
#include "stdlib.h"
#include "stddef.h"
#include "string.h"
#include "math.h"
#include "errno.h"

#include "signal.h"

#include "mpi.h"
#include "cuda.h"

#include "env.hh"

//=============================================================================
/*---Null object---*/

GMEnv GMEnv_null() {
  GMEnv result;
  memset((void*)&result, 0, sizeof(GMEnv));
  return result;
}

//=============================================================================
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
    if (is_delim_prev && ! is_delim) {
      argv[*argc] = &(argstring[i]);
      (*argc)++;
    }
    is_delim_prev = is_delim;
  }
}

//=============================================================================
/*---Initialize environment---*/

void GMEnv_create_impl_(GMEnv* const env, MPI_Comm base_comm, int argc,
                        char** argv, const char* const description,
                        bool make_comms, int num_proc, int proc_num) {
  GMInsist(env);

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
  GMEnv_ccc_param_set(GMEnv_ccc_param_default(), env);
  GMEnv_ccc_multiplier_set(GMEnv_ccc_multiplier_default(), env);

  env->time = 0;
  env->compares = 0;
  env->eltcompares = 0;
  env->veccompares = 0;
  env->ops_local = 0;
  env->ops = 0;
  env->cpu_mem = 0;
  env->cpu_mem_max = 0;
  env->gpu_mem = 0;
  env->gpu_mem_max = 0;
  env->description = description;
  env->tc = 0;
  env->num_tc_steps = 1;

  env->mpi_comm_base_ = base_comm;
  env->make_comms_ = make_comms;

  if (env->make_comms_) {
    int mpi_code = MPI_Comm_size(env->mpi_comm_base_, &env->num_proc_base_);
    GMInsist(mpi_code == MPI_SUCCESS);
    mpi_code = MPI_Comm_rank(env->mpi_comm_base_, &env->proc_num_base_);
    GMInsist(mpi_code == MPI_SUCCESS);
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
      GMInsistInterface(env, i < argc && "Missing value for metric_type.");
      if (strcmp(argv[i], "czekanowski") == 0) {
        env->metric_type_ = GM_METRIC_TYPE_CZEK;
      } else if (strcmp(argv[i], "ccc") == 0) {
        env->metric_type_ = GM_METRIC_TYPE_CCC;
      } else {
        GMInsistInterface(env, false && "Invalid setting for metric_type.");
      }
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_way") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_way.");
      errno = 0;
      const long num_way = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno && (num_way == GM_NUM_WAY_2 ||
                                   num_way == GM_NUM_WAY_3)
                               && "Invalid setting for num_way.");
      env->num_way_ = num_way;
      GMEnv_set_num_proc(env, env->num_proc_vector_i_, env->num_proc_repl_,
                       env->num_proc_field_);
      /*--------------------*/
    } else if (strcmp(argv[i], "--all2all") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for all2all.");
      if (strcmp(argv[i], "yes") == 0) {
        env->all2all_ = true;
      } else if (strcmp(argv[i], "no") == 0) {
        env->all2all_ = false;
      } else {
        GMInsistInterface(env, false && "Invalid setting for all2all.");
      }
      /*--------------------*/
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for compute_method.");
      if (strcmp(argv[i], "CPU") == 0) {
        GMEnv_set_compute_method(env, GM_COMPUTE_METHOD_CPU);
      } else if (strcmp(argv[i], "GPU") == 0) {
        GMEnv_set_compute_method(env, GM_COMPUTE_METHOD_GPU);
      } else if (strcmp(argv[i], "REF") == 0) {
        GMEnv_set_compute_method(env, GM_COMPUTE_METHOD_REF);
      } else {
        GMInsistInterface(env, false && "Invalid setting for compute_method.");
      }
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      /*--------------------*/
      ++i;
      errno = 0;
      GMInsistInterface(env, i < argc && "Missing value for num_proc_vector.");
      long num_proc_vector_i = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_proc_vector_i == num_proc_vector_i
                    && "Invalid setting for num_proc_vector.");
      GMEnv_set_num_proc(env, num_proc_vector_i, env->num_proc_repl_,
                         env->num_proc_field_);
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      /*--------------------*/
      ++i;
      errno = 0;
      GMInsistInterface(env, i < argc && "Missing value for num_proc_field.");
      long num_proc_field = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_proc_field == num_proc_field
                    && "Invalid setting for num_proc_field.");
      GMEnv_set_num_proc(env, env->num_proc_vector_i_, env->num_proc_repl_,
                         num_proc_field);
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      /*--------------------*/
      ++i;
      errno = 0;
      GMInsistInterface(env, i < argc && "Missing value for num_proc_repl.");
      long num_proc_repl = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_proc_repl == num_proc_repl
                    && "Invalid setting for num_proc_repl.");
      GMEnv_set_num_proc(env, env->num_proc_vector_i_, num_proc_repl,
                         env->num_proc_field_);
      /*--------------------*/
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for ccc_param.");
      errno = 0;
      const double ccc_param = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && ccc_param >= 0
                               && "Invalid setting for ccc_param.");
      GMEnv_ccc_param_set(ccc_param, env);
      /*--------------------*/
    } else if (strcmp(argv[i], "--ccc_multiplier") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for ccc_multiplier.");
      errno = 0;
      const double ccc_multiplier = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && ccc_multiplier >= 0
                               && "Invalid setting for ccc_multiplier.");
      GMEnv_ccc_multiplier_set(ccc_multiplier, env);
      /*--------------------*/
    } else if (strcmp(argv[i], "--sparse") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for sparse.");
      if (strcmp(argv[i], "yes") == 0) {
        env->sparse = true;
      } else if (strcmp(argv[i], "no") == 0) {
        env->sparse = false;
      } else {
        GMInsistInterface(env, false && "Invalid setting for sparse.");
      }
      /*--------------------*/
    } else if (strcmp(argv[i], "--tc") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for tc.");
      errno = 0;
      const long tc = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)tc == tc
                    && tc >= 0
                    && tc < GM_NUM_TC_GEMM_SOURCE_TYPE
                    && "Invalid setting for tc.");
      env->tc = tc;
      /*--------------------*/
    } else if (strcmp(argv[i], "--num_tc_steps") == 0) {
      /*--------------------*/
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_tc_steps.");
      errno = 0;
      const long num_tc_steps = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_tc_steps == num_tc_steps
                    && num_tc_steps >= 1
                    && "Invalid setting for tc.");
      env->num_tc_steps = num_tc_steps;
      /*--------------------*/
    } /*---if/else---*/
  }   /*---for i---*/

  /*---Helper variables---*/
  env->do_reduce = env->num_proc_field_ > 1;
  env->need_2way = env->metric_type_ == GM_METRIC_TYPE_CZEK;
}

//-----------------------------------------------------------------------------

void GMEnv_create(GMEnv* const env, MPI_Comm base_comm, int argc, char** argv,
                  const char* const description) {
  GMInsist(env);

  GMEnv_create_impl_(env, base_comm, argc, argv, description, true, 0, 0);
}

//-----------------------------------------------------------------------------

void GMEnv_create(GMEnv* const env, MPI_Comm base_comm,
                  const char* const options,
                  const char* const description) {
  GMInsist(env);

  /*---Convert options string to args---*/

  size_t len = strlen(options);
  char argstring[len+1];
  char* argv[len+1];
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);

  GMEnv_create_impl_(env, base_comm, argc, argv, description, true, 0, 0);
}

//-----------------------------------------------------------------------------

void GMEnv_create_no_comms(GMEnv* const env, int argc, char** argv,
                           const char* const description,
                           int num_proc, int proc_num) {
  GMInsist(env);

  GMEnv_create_impl_(env, MPI_COMM_WORLD, argc, argv, description,
                     false, num_proc, proc_num);
}

//-----------------------------------------------------------------------------

void GMEnv_create_no_comms(GMEnv* const env, const char* const options,
                           const char* const description,
                           int num_proc, int proc_num) {
  GMInsist(env);

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

//=============================================================================
/*---Manage cuda streams---*/

void GMEnv_initialize_streams(GMEnv* const env) {
  GMInsist(env);

  /*---NOTE: this is used for lazy initialization---*/

  if (env->are_cuda_streams_initialized_) {
    return;
  }

  if (env->compute_method_ != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  cudaStreamCreate(&env->stream_compute_);
  GMInsist(GMEnv_cuda_last_call_succeeded(env));

  cudaStreamCreate(&env->stream_togpu_);
  GMInsist(GMEnv_cuda_last_call_succeeded(env));

  cudaStreamCreate(&env->stream_fromgpu_);
  GMInsist(GMEnv_cuda_last_call_succeeded(env));

  env->are_cuda_streams_initialized_ = true;
}

//-----------------------------------------------------------------------------

void GMEnv_terminate_streams(GMEnv* const env) {
  GMInsist(env);

  if (! env->are_cuda_streams_initialized_) {
    return;
  }

  cudaStreamDestroy(env->stream_compute_);
  GMInsist(GMEnv_cuda_last_call_succeeded(env));

  cudaStreamDestroy(env->stream_togpu_);
  GMInsist(GMEnv_cuda_last_call_succeeded(env));

  cudaStreamDestroy(env->stream_fromgpu_);
  GMInsist(GMEnv_cuda_last_call_succeeded(env));

  env->are_cuda_streams_initialized_ = false;
}

//=============================================================================
/*---Manage MPI comms---*/

void GMEnv_initialize_comms(GMEnv* const env) {
  GMInsist(env);

  if (env->are_mpi_comms_initialized_) {
    return;
  }

  if (! env->make_comms_) {
    return;
  }

  int mpi_code = MPI_Comm_split(env->mpi_comm_base_, env->is_proc_active_,
                            env->proc_num_, &env->mpi_comm_);
  GMInsist(mpi_code == MPI_SUCCESS);

  // Communicator along repl / vector axis.

  mpi_code = MPI_Comm_split(env->mpi_comm_base_,
      env->is_proc_active_ ? env->proc_num_field_ : env->num_proc_,
      //env->proc_num_,
      env->is_proc_active_ ? env->proc_num_repl_vector_ : env->proc_num_,
      &env->mpi_comm_repl_vector_);
  GMInsist(mpi_code == MPI_SUCCESS);

  // Communicator along field axis.

  mpi_code = MPI_Comm_split(env->mpi_comm_base_,
      env->is_proc_active_ ? env->proc_num_repl_vector_ : env->num_proc_,
      //env->proc_num_,
      env->is_proc_active_ ? env->proc_num_field_ : env->proc_num_,
      &env->mpi_comm_field_);
  GMInsist(mpi_code == MPI_SUCCESS);

  env->are_mpi_comms_initialized_ = true;
}

//-----------------------------------------------------------------------------

void GMEnv_terminate_comms(GMEnv* const env) {
  GMInsist(env);

  if (! env->are_mpi_comms_initialized_) {
    return;
  }

  /*---Destroy any nontrivial communicators---*/

  int mpi_code = MPI_Comm_free(&(env->mpi_comm_));
  GMInsist(mpi_code == MPI_SUCCESS);

  mpi_code = MPI_Comm_free(&(env->mpi_comm_repl_vector_));
  GMInsist(mpi_code == MPI_SUCCESS);

  mpi_code = MPI_Comm_free(&(env->mpi_comm_field_));
  GMInsist(mpi_code == MPI_SUCCESS);

  env->are_mpi_comms_initialized_ = false;
}

//=============================================================================
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* const env) {
  GMInsist(env);

  GMEnv_terminate_comms(env);

  GMEnv_terminate_streams(env);

  *env = GMEnv_null();
}

//=============================================================================
/*---Accessors---*/

void GMEnv_set_compute_method(GMEnv* const env, int compute_method) {
  GMInsist(env);
  GMInsist(compute_method >= 0);
  GMInsist(compute_method < GM_NUM_COMPUTE_METHOD);

  env->compute_method_ = compute_method;
}

//-----------------------------------------------------------------------------

int GMEnv_data_type_vectors(const GMEnv* const env) {
  GMInsist(env);

  switch (env->metric_type_) {
    case GM_METRIC_TYPE_CZEK:
      return GM_DATA_TYPE_FLOAT;
    case GM_METRIC_TYPE_CCC:
      return GM_DATA_TYPE_BITS2;
  }
  GMInsist(false && "Invalid metric type.");
  return 0;
}

//-----------------------------------------------------------------------------

int GMEnv_data_type_metrics(const GMEnv* const env) {
  GMInsist(env);

  switch (env->metric_type_) {
    case GM_METRIC_TYPE_CZEK:
      return GM_DATA_TYPE_FLOAT;
    case GM_METRIC_TYPE_CCC:
      return env->num_way_ == GM_NUM_WAY_2 ? GM_DATA_TYPE_TALLY2X2
                                           : GM_DATA_TYPE_TALLY4X2;
  }
  GMInsist(false && "Invalid metric type.");
  return 0;
}

//-----------------------------------------------------------------------------

void GMEnv_set_num_proc(GMEnv* const env, int num_proc_vector_i,
                      int num_proc_repl, int num_proc_field) {
  GMInsist(env);
  GMInsist(num_proc_vector_i > 0);
  GMInsist(num_proc_repl > 0);
  GMInsist(num_proc_field > 0);

#ifdef NOMPI
  GMInsist(num_proc_vector_i == 1);
  GMInsist(num_proc_repl == 1);
  GMInsist(num_proc_field == 1);
#endif

  GMInsist(env->num_proc_base_ != 0);
  //GMInsist(env->proc_num_base_ is initialized);

  /*---Set proc counts---*/

  env->num_proc_vector_i_ = num_proc_vector_i;
  env->num_proc_repl_ = num_proc_repl;
  env->num_proc_field_ = num_proc_field;

  env->num_proc_repl_vector_ = env->num_proc_vector_i_ * env->num_proc_repl_;

  env->num_proc_ = env->num_proc_repl_vector_ * num_proc_field;
  GMInsist(env->num_proc_ <= env->num_proc_base_);

  /*---Set proc nums---*/

  env->proc_num_ = env->proc_num_base_;

  env->is_proc_active_ = env->proc_num_ < env->num_proc_;

  enum {ORDER_FRV = 0,
        ORDER_RVF = 1,
        ORDER_FVR = 2};

  //const int order = ORDER_FRV;
  const int order = ORDER_FVR;

  if (order == ORDER_FRV) {
    env->proc_num_field_ = env->proc_num_ % env->num_proc_field_;
    env->proc_num_repl_ = (env->proc_num_ / env->num_proc_field_)
                                          % env->num_proc_repl_;
    env->proc_num_vector_i_ = (env->proc_num_ / env->num_proc_field_)
                                              / env->num_proc_repl_;
  }

  if (order == ORDER_RVF) {
    env->proc_num_repl_ = env->proc_num_ % env->num_proc_repl_;
    env->proc_num_vector_i_ = (env->proc_num_ / env->num_proc_repl_)
                                              % env->num_proc_vector_i_;
    env->proc_num_field_ = (env->proc_num_ / env->num_proc_repl_)
                                           / env->num_proc_vector_i_;
  }

  if (order == ORDER_FVR) {
    env->proc_num_field_ = env->proc_num_ % env->num_proc_field_;
    env->proc_num_vector_i_ = (env->proc_num_ / env->num_proc_field_)
                                              % env->num_proc_vector_i_;
    env->proc_num_repl_ = (env->proc_num_ / env->num_proc_field_)
                                          / env->num_proc_vector_i_;
  }

  env->proc_num_repl_vector_ = env->proc_num_repl_ + env->num_proc_repl_ *
                               env->proc_num_vector_i_;

  /*---Destroy old communicators if necessary---*/

  GMEnv_terminate_comms(env);

  /*---Make new communicators---*/

  GMEnv_initialize_comms(env);
}

//-----------------------------------------------------------------------------

cudaStream_t GMEnv_stream_compute(GMEnv* const env) {
  GMInsist(env);
  // GMInsist(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_compute_;
}

//-----------------------------------------------------------------------------

cudaStream_t GMEnv_stream_togpu(GMEnv* const env) {
  GMInsist(env);
  // GMInsist(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_togpu_;
}

//-----------------------------------------------------------------------------

cudaStream_t GMEnv_stream_fromgpu(GMEnv* const env) {
  GMInsist(env);
  // GMInsist(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_fromgpu_;
}

//=============================================================================
/*---Timer functions---*/

double GMEnv_get_time(const GMEnv* const env) {

  struct timeval tv;
  gettimeofday(&tv, NULL);
  double result = ((double)tv.tv_sec + (double)tv.tv_usec * 1.e-6);

  return result;
}

//-----------------------------------------------------------------------------

double GMEnv_get_synced_time(const GMEnv* const env) {
  GMInsist(env);

  if (! GMEnv_is_proc_active(env)) {
    return 0;
  }

  /*
  cudaThreadSynchronize();
  */

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    cudaDeviceSynchronize();
    GMInsist(GMEnv_cuda_last_call_succeeded(env));
  }

  const int mpi_code = MPI_Barrier(GMEnv_mpi_comm(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  return GMEnv_get_time(env);
}

//=============================================================================
/*---Misc.---*/

bool GMEnv_cuda_last_call_succeeded(const GMEnv* const env) {
  GMInsist(env);

  bool result = true;

  /*---NOTE: this read of the last error is a destructive read---*/
  cudaError_t error = cudaGetLastError();

  if (error != cudaSuccess) {
    result = false;
    printf("CUDA error detected: %s\n", cudaGetErrorString(error));
  }

  return result;
}

//-----------------------------------------------------------------------------

int gm_gpu_compute_capability() {
  cudaDeviceProp deviceProp;

  // Assume only one GPU visible per rank.

  cudaGetDeviceProperties(&deviceProp, 0);

  return deviceProp.major * 100 + deviceProp.minor;
}

//-----------------------------------------------------------------------------

void* gm_malloc(size_t n, GMEnv* env) {
  GMInsist(env);
  void* p = malloc(n);
  GMInsist(p && "Invalid pointer from malloc,"
                " possibly due to insufficient memory.");
  env->cpu_mem += n;
  env->cpu_mem_max = gm_max_i8(env->cpu_mem_max, env->cpu_mem);
  return p;
}

//-----------------------------------------------------------------------------

void gm_free(void* p, size_t n, GMEnv* env) {
  GMInsist(p && env);
  free(p);
  env->cpu_mem -= n;
}

//-----------------------------------------------------------------------------

bool GMEnv_is_ppc64() {
#ifdef __powerpc64__
  return true;
#else
  return false;
#endif
   //return strcmp("__PPC64__", "__" "PPC64" "__") != 0;
}

//-----------------------------------------------------------------------------

GMFloat* GMFloat_malloc(size_t n, GMEnv* env) {
  GMInsist(env);
  GMFloat* p = (GMFloat*)gm_malloc(n * sizeof(GMFloat), env);
  GMInsist(p);
  GMFloat_fill_nan(p, n);
  return p;
}

//-----------------------------------------------------------------------------

void GMFloat_free(GMFloat* p, size_t n, GMEnv* env) {
  GMInsist(p && env);
  gm_free(p, n * sizeof(GMFloat), env);
}

//-----------------------------------------------------------------------------

void GMFloat_fill_nan(GMFloat* const a, size_t n) {
  GMInsist(a);
  GMInsist(n+1 >= 1);
#ifdef GM_ASSERTIONS_ON
  GMFloat value = sqrt(-1);
  size_t i = 0;
  for (i=0; i<n; ++i) {
    a[i] = value;
  }
#endif
}

//-----------------------------------------------------------------------------

void GMFloat_check(GMFloat* const a, size_t n) {
  GMInsist(a);
  GMInsist(n+1 >= 1);
#ifdef GM_ASSERTIONS_ON
  bool no_nans_found = true;
  size_t i = 0;
  for (i=0; i<n; ++i) {
    if (a[i] != a[i]) {
      no_nans_found = false;
    }
  }
  GMInsist(no_nans_found);
#endif
}

//-----------------------------------------------------------------------------

//int GMFloat_mant_dig() {
//  GMInsist(FLT_RADIX == 2);
//  return sizeof(GMFloat) == 8 ? DBL_MANT_DIG : FLT_MANT_DIG;
//}

template<> int gm_mant_dig<float>() {
  GMInsist(FLT_RADIX == 2);
  return FLT_MANT_DIG;
}

template<> int gm_mant_dig<double>() {
  GMInsist(FLT_RADIX == 2);
  return DBL_MANT_DIG;
}

//-----------------------------------------------------------------------------

MPI_Datatype gm_mpi_type(const GMEnv* const env) {
  GMInsist(env);

  /* clang-format off */
  const MPI_Datatype mpi_type =
    GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK ? GM_MPI_FLOAT :
    GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC ? MPI_DOUBLE_COMPLEX :
                                                   0; // should never get here
  /* clang-format on */

  return mpi_type;
}

//-----------------------------------------------------------------------------

size_t gm_num_vector_local_required(size_t num_vector_active_local,
                                    GMEnv* const env) {
  GMInsist(env);
  // NOTE: this function should receive the same num_vector_active_local
  // and give the same result independent of MPI rank.

  const bool need_divisible_by_6 = GMEnv_num_way(env) == GM_NUM_WAY_3 &&
                                   GMEnv_all2all(env) &&
                                   GMEnv_num_proc_vector_i(env) > 2;

  const bool need_divisible_by_4 = env->tc;

  const int round_factor = (need_divisible_by_4 && need_divisible_by_6) ? 12 :
                            need_divisible_by_4 ? 4 :
                            need_divisible_by_6 ? 6 : 1;

  const size_t num_vector_local = gm_ceil_i8(num_vector_active_local,
                                             round_factor)*round_factor;

  return num_vector_local;
}

//-----------------------------------------------------------------------------

size_t gm_gemm_size_required(size_t size_requested, GMEnv* const env) {
  GMInsist(env);

  const bool need_divisible_by_4 = env->tc;

  return need_divisible_by_4 ? gm_ceil_i8(size_requested, 4)*4 : size_requested;
}

//-----------------------------------------------------------------------------

size_t gm_array_cksum(unsigned char* a, size_t n) {
  GMInsist(a);

  size_t result = 0;

  const size_t mask = (((size_t)1) << 32) - 1;

#pragma omp parallel for schedule(dynamic,1000) reduction(+:result)
  for (size_t i=0; i<n; ++i) {
    result += (a[i] * i) & mask;
  }

  return result;
}

//-----------------------------------------------------------------------------
