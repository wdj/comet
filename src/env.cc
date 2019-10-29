//-----------------------------------------------------------------------------
/*!
 * \file   env.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Basic environment - settings, MPI communicators, etc.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdio"
#include "cstdlib"
#include "cstddef"
#include "string.h"
#include "math.h"

#include "errno.h"
#include "sys/time.h"
#include "signal.h"

#include "mpi.h"

#if defined USE_CUDA
#include "cuda.h"
#endif

#include "env.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Null object

GMEnv GMEnv_null() {
  GMEnv result;
  memset((void*)&result, 0, sizeof(GMEnv));
  return result;
}

//-----------------------------------------------------------------------------

void Env::set_defaults() {

  *this = GMEnv_null();

  // CoMet Settings
  metric_type_ = GM_METRIC_TYPE_CZEK;
  num_way_ = GM_NUM_WAY_2;
  all2all_ = false;
  GMEnv_set_compute_method(this, GM_COMPUTE_METHOD_GPU);
  num_stage = 1;
  stage_num = 0;
  num_phase = 1;
  phase_num = 0;
  GMEnv_ccc_param_set(GMEnv_ccc_param_default(), this);
  GMEnv_ccc_multiplier_set(GMEnv_ccc_multiplier_default(), this);
  GMEnv_duo_multiplier_set(GMEnv_duo_multiplier_default(), this);
  sparse = false;
  tc_ = Env::TC_NONE;
  num_tc_steps = 1;
  // Counters
  //env->ctime_ = 0;
  compares = 0;
  eltcompares = 0;
  veccompares = 0;
  ops_local = 0;
  ops = 0;
  cpu_mem_local = 0;
  cpu_mem_max_local = 0;
  cpu_mem_max = 0;
  gpu_mem_local = 0;
  gpu_mem_max_local = 0;
  gpu_mem_max = 0;
  // MPI
  are_comms_initialized_ = false;
}

//-----------------------------------------------------------------------------

void Env::parse_args(int argc, char** argv) {

  Env* env = this;

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--metric_type") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for metric_type.");
      if (strcmp(argv[i], "czekanowski") == 0) {
        metric_type_ = GM_METRIC_TYPE_CZEK;
      } else if (strcmp(argv[i], "ccc") == 0) {
        metric_type_ = GM_METRIC_TYPE_CCC;
      } else if (strcmp(argv[i], "duo") == 0) {
        metric_type_ = GM_METRIC_TYPE_DUO;
      } else {
        GMInsistInterface(env, false && "Invalid setting for metric_type.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--num_way") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_way.");
      errno = 0;
      const long num_way = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno && (num_way == GM_NUM_WAY_2 ||
                                            num_way == GM_NUM_WAY_3)
                               && "Invalid setting for num_way.");
      num_way_ = num_way;
      GMEnv_set_num_proc(env, num_proc_vector_i_, num_proc_repl_,
                       num_proc_field_);
      //--------------------
    } else if (strcmp(argv[i], "--all2all") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for all2all.");
      if (strcmp(argv[i], "yes") == 0) {
        all2all_ = true;
      } else if (strcmp(argv[i], "no") == 0) {
        all2all_ = false;
      } else {
        GMInsistInterface(env, false && "Invalid setting for all2all.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      //--------------------
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
      //--------------------
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      //--------------------
      ++i;
      errno = 0;
      GMInsistInterface(env, i < argc && "Missing value for num_proc_vector.");
      long num_proc_vector_i = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_proc_vector_i == num_proc_vector_i
                    && "Invalid setting for num_proc_vector.");
      GMEnv_set_num_proc(env, num_proc_vector_i, num_proc_repl_,
                         num_proc_field_);
      //--------------------
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      //--------------------
      ++i;
      errno = 0;
      GMInsistInterface(env, i < argc && "Missing value for num_proc_field.");
      long num_proc_field = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_proc_field == num_proc_field
                    && "Invalid setting for num_proc_field.");
      GMEnv_set_num_proc(env, num_proc_vector_i_, num_proc_repl_,
                         num_proc_field);
      //--------------------
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      //--------------------
      ++i;
      errno = 0;
      GMInsistInterface(env, i < argc && "Missing value for num_proc_repl.");
      long num_proc_repl = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_proc_repl == num_proc_repl
                    && "Invalid setting for num_proc_repl.");
      GMEnv_set_num_proc(env, num_proc_vector_i_, num_proc_repl,
                         num_proc_field_);
      //--------------------
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for ccc_param.");
      errno = 0;
      const double ccc_param = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && ccc_param >= 0
                               && "Invalid setting for ccc_param.");
      GMEnv_ccc_param_set(ccc_param, env);
      //--------------------
    } else if (strcmp(argv[i], "--ccc_multiplier") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for ccc_multiplier.");
      errno = 0;
      const double ccc_multiplier = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && ccc_multiplier >= 0
                               && "Invalid setting for ccc_multiplier.");
      GMEnv_ccc_multiplier_set(ccc_multiplier, env);
      //--------------------
    } else if (strcmp(argv[i], "--duo_multiplier") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for duo_multiplier.");
      errno = 0;
      const double duo_multiplier = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && duo_multiplier >= 0
                               && "Invalid setting for duo_multiplier.");
      GMEnv_duo_multiplier_set(duo_multiplier, env);
      //--------------------
    } else if (strcmp(argv[i], "--sparse") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for sparse.");
      if (strcmp(argv[i], "yes") == 0) {
        sparse = true;
      } else if (strcmp(argv[i], "no") == 0) {
        sparse = false;
      } else {
        GMInsistInterface(env, false && "Invalid setting for sparse.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--tc") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for tc.");
      errno = 0;
      const long tc = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)tc == tc
                    && tc >=0 && tc < Env::NUM_TC
                    && "Invalid setting for tc.");
      tc_ = tc;
      //--------------------
    } else if (strcmp(argv[i], "--num_tc_steps") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for num_tc_steps.");
      errno = 0;
      const long num_tc_steps = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_tc_steps == num_tc_steps
                    && num_tc_steps >= 1
                    && "Invalid setting for tc.");
      num_tc_steps = num_tc_steps;
      //--------------------
    } // if/else
  }   // for i
}

//=============================================================================
// Utility to parse a string to construct arguments

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
// Initialize environment

void GMEnv_create_impl_(GMEnv* const env, MPI_Comm base_comm, int argc,
                        char** argv, const char* const description,
                        bool make_comms, int num_proc, int proc_num) {
  GMInsist(env);

  env->set_defaults();

  env->comm_base_ = base_comm;
  env->make_comms_ = make_comms;
  // MPI proc counts, proc numbers
  if (env->make_comms_) {
    COMET_MPI_SAFE_CALL(MPI_Comm_size(env->comm_base_, &env->num_proc_base_));
    COMET_MPI_SAFE_CALL(MPI_Comm_rank(env->comm_base_, &env->proc_num_base_));
  } else {
    env->num_proc_base_ = num_proc;
    env->proc_num_base_ = proc_num;
  }
  GMEnv_set_num_proc(env, env->num_proc_base_, 1, 1);
  // accelerator streams
  //env->are_streams_initialized_ = false;
  //env->are_comms_initialized_ = false;
  // OTHER
  env->description = description;

  env->parse_args(argc, argv);

  if (make_comms) {
    GMInsistInterface(env, env->can_run() &&
                      "Invalid problem for this system and build.");
  }

  // Helper variables
  env->do_reduce = env->num_proc_field_ > 1;
  env->need_2way = env->metric_type_ == GM_METRIC_TYPE_CZEK;
}

//-----------------------------------------------------------------------------

void GMEnv_create(GMEnv* const env, MPI_Comm base_comm,
                  int argc, char** argv,
                  const char* const description) {
  GMInsist(env);

  const int num_proc = 0;
  const int proc_num = 0;
  const bool make_comms = true;

  GMEnv_create_impl_(env, base_comm, argc, argv, description, make_comms,
    num_proc, proc_num);
}

//-----------------------------------------------------------------------------

void GMEnv_create(GMEnv* const env, MPI_Comm base_comm,
                  const char* const options,
                  const char* const description) {
  GMInsist(env);

  // Convert options string to args

  size_t len = strlen(options);

  char* argstring = (char*)malloc((len+1)*sizeof(char));
  char ** argv = (char**)malloc((len+1)*sizeof(char*));
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);

  const int num_proc = 0;
  const int proc_num = 0;
  const bool make_comms = true;

  GMEnv_create_impl_(env, base_comm, argc, argv, description, make_comms, 
    num_proc, proc_num);

  free(argstring);
  free(argv);
}

//-----------------------------------------------------------------------------

void GMEnv_create(GMEnv* const env,
                  const char* const options,
                  int num_proc, int proc_num) {
  GMInsist(env);

  // Convert options string to args

  size_t len = strlen(options);
  char* argstring = (char*)malloc((len+1)*sizeof(char));
  char ** argv = (char**)malloc((len+1)*sizeof(char*));
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);

  const char* const description = NULL;
  const bool make_comms = false;

  GMEnv_create_impl_(env, MPI_COMM_WORLD, argc, argv, description,
                     make_comms, num_proc, proc_num);

  free(argstring);
  free(argv);
}

//-----------------------------------------------------------------------------

Env::Env(const char* const options, int num_proc, int proc_num) {

  // Convert options string to args

  size_t len = strlen(options);
  char* argstring = (char*)malloc((len+1)*sizeof(char));
  char ** argv = (char**)malloc((len+1)*sizeof(char*));
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);
  
  GMEnv_create_impl_(this, MPI_COMM_WORLD, argc, argv, NULL,
                     false, num_proc, proc_num);
  
  free(argstring);
  free(argv);
}

//-----------------------------------------------------------------------------

Env::~Env() {
  comms_terminate_();
  streams_terminate_();
}

//=============================================================================
// Manage accelerator streams

void Env::streams_initialize_() {

  // NOTE: this is used for lazy initialization

  if (are_streams_initialized_)
    return;

  if (compute_method_ != GM_COMPUTE_METHOD_GPU)
    return;

#if defined USE_CUDA
  cudaStreamCreate(&stream_compute_);
#elif defined USE_HIP
  hipStreamCreate(&stream_compute_);
#endif
  GMInsist(GMEnv_accel_last_call_succeeded(this) &&
           "Failure in call to stream create.");

#if defined USE_CUDA
  cudaStreamCreate(&stream_togpu_);
#elif defined USE_HIP
  hipStreamCreate(&stream_togpu_);
#endif
  GMInsist(GMEnv_accel_last_call_succeeded(this) &&
           "Failure in call to stream create.");

#if defined USE_CUDA
  cudaStreamCreate(&stream_fromgpu_);
#endif
#ifdef USE_HIP
  hipStreamCreate(&stream_fromgpu_);
#endif
  GMInsist(GMEnv_accel_last_call_succeeded(this) &&
           "Failure in call to stream create.");

  are_streams_initialized_ = true;
}

//-----------------------------------------------------------------------------

void Env::streams_terminate_() {

  if (! are_streams_initialized_)
    return;

#if defined USE_CUDA
  cudaStreamDestroy(stream_compute_);
#elif defined USE_HIP
  hipStreamDestroy(stream_compute_);
#endif
  GMInsist(GMEnv_accel_last_call_succeeded(this) &&
           "Failure in call to stream destroy.");

#if defined USE_CUDA
  cudaStreamDestroy(stream_togpu_);
#elif defined USE_HIP
  hipStreamDestroy(stream_togpu_);
#endif
  GMInsist(GMEnv_accel_last_call_succeeded(this) &&
           "Failure in call to stream destroy.");

#if defined USE_CUDA
  cudaStreamDestroy(stream_fromgpu_);
#elif defined USE_HIP
  hipStreamDestroy(stream_fromgpu_);
#endif
  GMInsist(GMEnv_accel_last_call_succeeded(this) &&
           "Failure in call to stream destroy.");

  are_streams_initialized_ = false;
}

//-----------------------------------------------------------------------------

void Env::stream_synchronize(Stream_t stream) const {

  if (GMEnv_compute_method(this) != GM_COMPUTE_METHOD_GPU)
    return;

#if defined USE_CUDA
  cudaStreamSynchronize(stream);
#elif defined USE_HIP
  hipStreamSynchronize(stream);
#endif
  GMInsist(GMEnv_accel_last_call_succeeded(this) &&
           "Failure in call to stream synchronize.");
}

//=============================================================================
// Manage MPI comms

void GMEnv_initialize_comms(GMEnv* const env) {
  GMInsist(env);

  if (env->are_comms_initialized_)
    return;

  if (! env->make_comms_)
    return;

  int mpi_code = MPI_Comm_split(env->comm_base_, env->is_proc_active_,
                            env->proc_num_, &env->comm_);
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Comm_split.");

  // Communicator along repl / vector axis.

  mpi_code = MPI_Comm_split(env->comm_base_,
      env->is_proc_active_ ? env->proc_num_field_ : env->num_proc_,
      //env->proc_num_,
      env->is_proc_active_ ? env->proc_num_repl_vector_ : env->proc_num_,
      &env->comm_repl_vector_);
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Comm_split.");

  // Communicator along field axis.

  mpi_code = MPI_Comm_split(env->comm_base_,
      env->is_proc_active_ ? env->proc_num_repl_vector_ : env->num_proc_,
      //env->proc_num_,
      env->is_proc_active_ ? env->proc_num_field_ : env->proc_num_,
      &env->comm_field_);
  GMInsist(mpi_code == MPI_SUCCESS && "Failure in call to MPI_Comm_split.");

  env->are_comms_initialized_ = true;
}

//-----------------------------------------------------------------------------

void Env::comms_terminate_() {

  if (! are_comms_initialized_)
    return;

  // Destroy any nontrivial communicators

  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_));
  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_repl_vector_));
  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_field_));

  are_comms_initialized_ = false;
}

//=============================================================================
// Finalize environment

void GMEnv_destroy(GMEnv* const env) {
  GMInsist(env);

  env->comms_terminate_();
  env->streams_terminate_();
  *env = GMEnv_null();
}

//=============================================================================
// Accessors

void GMEnv_set_compute_method(GMEnv* const env, int compute_method) {
  GMInsist(env);
  GMInsist(compute_method >= 0 && compute_method < GM_NUM_COMPUTE_METHOD);

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
    case GM_METRIC_TYPE_DUO:
      return GM_DATA_TYPE_BITS2;
  }
  GMInsist(false && "Invalid metric_type.");
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
    case GM_METRIC_TYPE_DUO:
      return GM_DATA_TYPE_TALLY2X2; // 2-way only for now
  }
  GMInsist(false && "Invalid metric_type.");
  return 0;
}

//-----------------------------------------------------------------------------

void GMEnv_set_num_proc(GMEnv* const env, int num_proc_vector_i,
                      int num_proc_repl, int num_proc_field) {
  GMInsist(env);
  GMInsist(num_proc_vector_i > 0);
  GMInsist(num_proc_repl > 0);
  GMInsist(num_proc_field > 0);

#ifndef USE_MPI
  GMInsist(num_proc_vector_i == 1);
  GMInsist(num_proc_repl == 1);
  GMInsist(num_proc_field == 1);
#endif

  GMInsist(env->num_proc_base_ != 0);
  //GMInsist(env->proc_num_base_ is initialized);

  // Set proc counts

  env->num_proc_vector_i_ = num_proc_vector_i;
  env->num_proc_repl_ = num_proc_repl;
  env->num_proc_field_ = num_proc_field;

  env->num_proc_repl_vector_ = env->num_proc_vector_i_ * env->num_proc_repl_;

  env->num_proc_ = env->num_proc_repl_vector_ * num_proc_field;
  GMInsist(env->num_proc_ <= env->num_proc_base_ &&
           "Number of procs requested exceeds number available.");

  // Set proc nums

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

  // Destroy old communicators if necessary

  env->comms_terminate_();

  // Make new communicators

  GMEnv_initialize_comms(env);
}

//-----------------------------------------------------------------------------

Stream_t Env::stream_compute() {
  streams_initialize_();
  return stream_compute_;
}

Stream_t Env::stream_togpu() {
  streams_initialize_();
  return stream_togpu_;
}

Stream_t Env::stream_fromgpu() {
  streams_initialize_();
  return stream_fromgpu_;
}

//=============================================================================
// Timer functions

double Env::time() {

  struct timeval tv;
  gettimeofday(&tv, NULL);
  double result = ((double)tv.tv_sec + (double)tv.tv_usec * 1.e-6);

  return result;
}

//-----------------------------------------------------------------------------

void GMEnv_accel_sync(const GMEnv* const env) {
  GMInsist(env);

  if (! env->is_proc_active()) {
    return;
  }

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

#if defined USE_CUDA
  cudaDeviceSynchronize();
#elif defined USE_HIP
  hipDeviceSynchronize();
#endif
  GMInsist(GMEnv_accel_last_call_succeeded(env) &&
           "Failure in call to device synchronize.");
}

//-----------------------------------------------------------------------------

double Env::synced_time() {

  if (! is_proc_active())
    return 0;

  GMEnv_accel_sync(this);

  COMET_MPI_SAFE_CALL(MPI_Barrier(comm_));

  return Env::time();
}

//=============================================================================
// Memory, arrays and floating point

void gm_cpu_mem_inc(size_t n, GMEnv* env) {
  GMInsist(env);
  env->cpu_mem_local += n;
  env->cpu_mem_max_local = gm_max_i8(env->cpu_mem_max_local,
                                     env->cpu_mem_local);
}

//----------

void gm_cpu_mem_dec(size_t n, GMEnv* env) {
  GMInsist(env);
  env->cpu_mem_local -= n;
}

//----------

void gm_gpu_mem_inc(size_t n, GMEnv* env) {
  GMInsist(env);
  env->gpu_mem_local += n;
  env->gpu_mem_max_local = gm_max_i8(env->gpu_mem_max_local,
                                     env->gpu_mem_local);
}

//----------

void gm_gpu_mem_dec(size_t n, GMEnv* env) {
  GMInsist(env);
  env->gpu_mem_local -= n;
}

//-----------------------------------------------------------------------------

void* gm_malloc(size_t n, GMEnv* env) {
  GMInsist(env);
  void* p = malloc(n);
  GMInsist(p &&
           "Invalid pointer from malloc, possibly due to insufficient memory.");
  gm_cpu_mem_inc(n, env);
  return p;
}

//-----------------------------------------------------------------------------

void gm_free(void* p, size_t n, GMEnv* env) {
  GMInsist(p && env);
  free(p);
  gm_cpu_mem_dec(n, env);
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
    GMEnv_metric_type(env) == GM_METRIC_TYPE_DUO ? MPI_DOUBLE_COMPLEX :
                                                   0; // should never get here
  /* clang-format on */

  return mpi_type;
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

//=============================================================================
// Misc.

bool GMEnv_accel_last_call_succeeded(const GMEnv* const env) {
  GMInsist(env);

#if defined USE_CUDA
  // NOTE: this read of the last error is a destructive read.
  cudaError_t error = cudaGetLastError();
  const bool result = error == cudaSuccess;

  if (!result) {
    printf("CUDA error detected: %s\n", cudaGetErrorString(error));
  }

  return result;
#elif defined USE_HIP
  // NOTE: this read of the last error is (apparently) a destructive read.
  hipError_t error = hipGetLastError();
  const bool result = error == hipSuccess;

  if (!result) {
    printf("HIP error detected: %s\n", hipGetErrorString(error));
  }

  return result;
#endif

  return true;
}

//-----------------------------------------------------------------------------

int Env::num_proc() {
  int num_proc = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_size(MPI_COMM_WORLD, &num_proc));
  return num_proc;
}

//-----------------------------------------------------------------------------

int Env::compute_capability() {
#if defined USE_CUDA
  cudaDeviceProp deviceProp;
  // Assume only one GPU per rank.
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, 0);
  const int compute_capability = error != cudaSuccess ? 0 :
    deviceProp.major * 100 + deviceProp.minor;
#elif defined USE_HIP
  hipDeviceProp_t deviceProp;
  hipGetDeviceProperties(&deviceProp, 0); // Assume only one GPU per rank.
//FIX this
  const int compute_capability = deviceProp.major * 100 + deviceProp.minor;
#else
  const int compute_capability = 0;
#endif
  return compute_capability;
}

//-----------------------------------------------------------------------------
/// \brief Determine whether requested run possible on this hardware and build.

bool Env::can_run() const {

  bool result = true;

  result = result && num_proc_ <= Env::num_proc();

  if (num_proc_ > 1) {
    result = result && Env::HAVE_MPI;
  }

  if (is_metric_type_bitwise() && compute_method_ == Env::COMPUTE_METHOD_CPU) {
    result = result && (Env::TC_NONE == tc_ ||
                        (Env::TC_FLOAT32 == tc_ && Env::HAVE_CPUBLAS));
  }

  if (compute_method_ == Env::COMPUTE_METHOD_GPU) {
    result = result && Env::HAVE_ACCEL && Env::compute_capability() > 0;
  }

  if (is_metric_type_bitwise() && compute_method_ == Env::COMPUTE_METHOD_GPU &&
      Env::TC_FLOAT16 == tc_) {
    result = result && Env::HAVE_CUDA && Env::compute_capability() >= 700;
  }

  if (is_metric_type_bitwise() && compute_method_ == Env::COMPUTE_METHOD_GPU &&
      Env::TC_INT8 == tc_) {
    result = result && Env::HAVE_CUDA && Env::compute_capability() >= 750;
  }

  if (is_metric_type_bitwise() && compute_method_ == Env::COMPUTE_METHOD_GPU &&
      Env::TC_FLOAT32 == tc_) {
    result = result && (Env::HAVE_CUDA || Env::HAVE_HIP);
  }

  const bool tc_on = Env::TC_FLOAT16 == tc_ || Env::TC_INT8 == tc_ ||
                     Env::TC_FLOAT32 == tc_;

  if (compute_method_ == Env::COMPUTE_METHOD_GPU && (!is_metric_type_bitwise()
      || (is_metric_type_bitwise() && tc_on))) {
    result = result && Env::HAVE_MAGMA;
  }

  return result;
}

//-----------------------------------------------------------------------------

int Env::tc_eff() const {

  return tc_; // FIX
}


#if 0
//-----------------------------------------------------------------------------

bool Env::is_tc_valid(int tc) {
  if (tc == Env::TC_NONE)
    return true;
  if (!(tc == TC_FLOAT16 || tc == TC_INT8 || TC_FLOAT32))
    return false;

#if defined USE_CUDA
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0); // Assume only one GPU per rank.
  const int compute_capability = deviceProp.major * 100 + deviceProp.minor;
  if (tc == TC_FLOAT16 && compute_capability >= 700) return true;
  if (tc == TC_INT8 && compute_capability >= 700) return true;
  //if (tc == GM_TC_METHOD_INT4 && compute_capability < 750) return false;
  //if (tc == GM_TC_METHOD_INT1 && compute_capability < 750) return false;
#elif defined USE_HIP
  hipDeviceProp_t deviceProp;
  hipGetDeviceProperties(&deviceProp, 0); // Assume only one GPU per rank.
//FIX this
//  const int compute_capability = deviceProp.major * 100 + deviceProp.minor;
//  if (tc == GM_TC_METHOD_FLOAT16 && compute_capability >= 700) return true;
//  if (tc == GM_TC_METHOD_INT8 && compute_capability >= 700) return true;
#endif

  return false;
}
#endif

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
