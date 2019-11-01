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
#include "cstring"
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

void Env::set_defaults() {

  memset((void*)this, 0, sizeof(*this));

  // CoMet Settings
  metric_type_ = MetricType::CZEK;
  num_way_ = NUM_WAY::_2;
  all2all_ = false;
  compute_method_ = ComputeMethod::GPU;
  num_stage = 1;
  stage_num = 0;
  num_phase = 1;
  phase_num = 0;
  ccc_param_set_(Env::ccc_param_default());
  ccc_multiplier_set_(Env::ccc_multiplier_default());
  duo_multiplier_set_(Env::duo_multiplier_default());
  sparse_ = false;
  tc_ = TC::NONE;
  num_tc_steps_ = 1;
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
    if (false) {
      //--------------------
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for metric_type.");
      if (strcmp(argv[i], "czekanowski") == 0) {
        metric_type_ = MetricType::CZEK;
      } else if (strcmp(argv[i], "ccc") == 0) {
        metric_type_ = MetricType::CCC;
      } else if (strcmp(argv[i], "duo") == 0) {
        metric_type_ = MetricType::DUO;
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
      GMInsistInterface(env, 0 == errno && (num_way == NUM_WAY::_2 ||
                                            num_way == NUM_WAY::_3)
                               && "Invalid setting for num_way.");
      num_way_ = num_way;
      set_num_proc_(num_proc_vector_, num_proc_repl_, num_proc_field_);
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
        compute_method_ =  ComputeMethod::CPU;
      } else if (strcmp(argv[i], "GPU") == 0) {
        compute_method_ =  ComputeMethod::GPU;
      } else if (strcmp(argv[i], "REF") == 0) {
        compute_method_ =  ComputeMethod::REF;
      } else {
        GMInsistInterface(env, false && "Invalid setting for compute_method.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      //--------------------
      ++i;
      errno = 0;
      GMInsistInterface(env, i < argc && "Missing value for num_proc_vector.");
      long num_proc_vector = strtol(argv[i], NULL, 10);
      GMInsistInterface(env, 0 == errno
                    && (long)(int)num_proc_vector == num_proc_vector
                    && "Invalid setting for num_proc_vector.");
      set_num_proc_(num_proc_vector, num_proc_repl_, num_proc_field_);
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
      set_num_proc_(num_proc_vector_, num_proc_repl_, num_proc_field);
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
      set_num_proc_(num_proc_vector_, num_proc_repl, num_proc_field_);
      //--------------------
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for ccc_param.");
      errno = 0;
      const double ccc_param = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && ccc_param >= 0
                               && "Invalid setting for ccc_param.");
      env->ccc_param_set_(ccc_param);
      //--------------------
    } else if (strcmp(argv[i], "--ccc_multiplier") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for ccc_multiplier.");
      errno = 0;
      const double ccc_multiplier = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && ccc_multiplier >= 0
                               && "Invalid setting for ccc_multiplier.");
      ccc_multiplier_set_(ccc_multiplier);
      //--------------------
    } else if (strcmp(argv[i], "--duo_multiplier") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for duo_multiplier.");
      errno = 0;
      const double duo_multiplier = strtod(argv[i], NULL);
      GMInsistInterface(env, 0 == errno && duo_multiplier >= 0
                               && "Invalid setting for duo_multiplier.");
      duo_multiplier_set_(duo_multiplier);
      //--------------------
    } else if (strcmp(argv[i], "--sparse") == 0) {
      //--------------------
      ++i;
      GMInsistInterface(env, i < argc && "Missing value for sparse.");
      if (strcmp(argv[i], "yes") == 0) {
        sparse_ = true;
      } else if (strcmp(argv[i], "no") == 0) {
        sparse_ = false;
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
                    && TC::is_valid(tc)
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
      num_tc_steps_ = num_tc_steps;
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

void Env::create_impl_(MPI_Comm base_comm, int argc,
                       char** argv, const char* const description,
                       bool make_comms, int num_proc, int proc_num) {
  set_defaults();

  comm_base_ = base_comm;
  make_comms_ = make_comms;
  // MPI proc counts, proc numbers
  if (make_comms_) {
    COMET_MPI_SAFE_CALL(MPI_Comm_size(comm_base_, &num_proc_base_));
    COMET_MPI_SAFE_CALL(MPI_Comm_rank(comm_base_, &proc_num_base_));
  } else {
    num_proc_base_ = num_proc;
    proc_num_base_ = proc_num;
  }
  set_num_proc_(num_proc_base_, 1, 1);
  // OTHER
  this->description_ = description; //FIX

  parse_args(argc, argv);

  if (make_comms) {
    GMInsistInterface(this, can_run() &&
                      "Invalid problem for this system and build.");
  }
}

//-----------------------------------------------------------------------------

Env::Env(MPI_Comm base_comm,
        int argc,
        char** argv,
        const char* const description) {

  const int num_proc = 0;
  const int proc_num = 0;
  const bool make_comms = true;

  create_impl_(base_comm, argc, argv, description, make_comms,
    num_proc, proc_num);
}

//-----------------------------------------------------------------------------

Env::Env(MPI_Comm base_comm,
         const char* const options,
         const char* const description) {

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

  create_impl_(base_comm, argc, argv, description, make_comms, 
    num_proc, proc_num);

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
  
  create_impl_(MPI_COMM_WORLD, argc, argv, NULL, false, num_proc, proc_num);
  
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

  if (compute_method_ != ComputeMethod::GPU)
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

  if (compute_method() != ComputeMethod::GPU)
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

void Env::comms_initialize_() {

  if (are_comms_initialized_)
    return;

  if (! make_comms_)
    return;

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_, is_proc_active_,
    proc_num_, &comm_));

  // Communicator along repl / vector axis.

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_,
      is_proc_active_ ? proc_num_field_ : num_proc_,
      //proc_num_,
      is_proc_active_ ? proc_num_repl_vector_ : proc_num_,
      &comm_repl_vector_));

  // Communicator along field axis.

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_,
      is_proc_active_ ? proc_num_repl_vector_ : num_proc_,
      //proc_num_,
      is_proc_active_ ? proc_num_field_ : proc_num_,
      &comm_field_));

  are_comms_initialized_ = true;
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
// Accessors

int GMEnv_data_type_vectors(const GMEnv* const env) {
  GMInsist(env);

  switch (env->metric_type()) {
    case MetricType::CZEK:
      return GM_DATA_TYPE_FLOAT;
    case MetricType::CCC:
      return GM_DATA_TYPE_BITS2;
    case MetricType::DUO:
      return GM_DATA_TYPE_BITS2;
  }
  GMInsist(false && "Invalid metric_type.");
  return 0;
}

//-----------------------------------------------------------------------------

int GMEnv_data_type_metrics(const GMEnv* const env) {
  GMInsist(env);

  switch (env->metric_type()) {
    case MetricType::CZEK:
      return GM_DATA_TYPE_FLOAT;
    case MetricType::CCC:
      return env->num_way() == NUM_WAY::_2 ? GM_DATA_TYPE_TALLY2X2
                                           : GM_DATA_TYPE_TALLY4X2;
    case MetricType::DUO:
      return GM_DATA_TYPE_TALLY2X2; // 2-way only for now
  }
  GMInsist(false && "Invalid metric_type.");
  return 0;
}

//-----------------------------------------------------------------------------

void Env::set_num_proc_(int num_proc_vector,
                      int num_proc_repl, int num_proc_field) {
  GMInsist(num_proc_vector > 0);
  GMInsist(num_proc_repl > 0);
  GMInsist(num_proc_field > 0);

  if (!BuildHas::MPI) {
    GMInsist(num_proc_vector == 1);
    GMInsist(num_proc_repl == 1);
    GMInsist(num_proc_field == 1);
  }

  GMInsist(num_proc_base_ != 0);
  //GMInsist(proc_num_base_ is initialized);

  // Set proc counts

  num_proc_vector_ = num_proc_vector;
  num_proc_repl_ = num_proc_repl;
  num_proc_field_ = num_proc_field;

  num_proc_repl_vector_ = num_proc_vector_ * num_proc_repl_;

  num_proc_ = num_proc_repl_vector_ * num_proc_field;
  GMInsist(num_proc_ <= num_proc_base_ &&
           "Number of procs requested exceeds number available.");

  // Set proc nums

  proc_num_ = proc_num_base_;

  is_proc_active_ = proc_num_ < num_proc_;

  enum {ORDER_FRV = 0,
        ORDER_RVF = 1,
        ORDER_FVR = 2};

  //const int order = ORDER_FRV;
  const int order = ORDER_FVR;

  if (order == ORDER_FRV) {
    proc_num_field_ = proc_num_ % num_proc_field_;
    proc_num_repl_ = (proc_num_ / num_proc_field_)
                                          % num_proc_repl_;
    proc_num_vector_ = (proc_num_ / num_proc_field_)
                                              / num_proc_repl_;
  }

  if (order == ORDER_RVF) {
    proc_num_repl_ = proc_num_ % num_proc_repl_;
    proc_num_vector_ = (proc_num_ / num_proc_repl_)
                                              % num_proc_vector_;
    proc_num_field_ = (proc_num_ / num_proc_repl_)
                                           / num_proc_vector_;
  }

  if (order == ORDER_FVR) {
    proc_num_field_ = proc_num_ % num_proc_field_;
    proc_num_vector_ = (proc_num_ / num_proc_field_)
                                              % num_proc_vector_;
    proc_num_repl_ = (proc_num_ / num_proc_field_)
                                          / num_proc_vector_;
  }

  proc_num_repl_vector_ = proc_num_repl_ + num_proc_repl_ *
                               proc_num_vector_;

  // Destroy old communicators if necessary

  comms_terminate_();

  // Make new communicators

  comms_initialize_();
}

//-----------------------------------------------------------------------------

Env::Stream_t Env::stream_compute() {
  streams_initialize_();
  return stream_compute_;
}

Env::Stream_t Env::stream_togpu() {
  streams_initialize_();
  return stream_togpu_;
}

Env::Stream_t Env::stream_fromgpu() {
  streams_initialize_();
  return stream_fromgpu_;
}

//=============================================================================
// Timer functions

double System::time() {

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

  if (env->compute_method() != ComputeMethod::GPU) {
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

  return System::time();
}

//=============================================================================
// Memory, arrays and floating point

void* gm_malloc(size_t n, GMEnv* env) {
  GMInsist(env);
  void* p = malloc(n);
  GMInsist(p &&
           "Invalid pointer from malloc, possibly due to insufficient memory.");
  env->cpu_mem_inc(n);
  return p;
}

//-----------------------------------------------------------------------------

void gm_free(void* p, size_t n, GMEnv* env) {
  GMInsist(p && env);
  free(p);
  env->cpu_mem_dec(n);
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
    env->metric_type() == MetricType::CZEK ? GM_MPI_FLOAT :
    env->metric_type() == MetricType::CCC ? MPI_DOUBLE_COMPLEX :
    env->metric_type() == MetricType::DUO ? MPI_DOUBLE_COMPLEX :
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

int System::num_proc() {
  int num_proc = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_size(MPI_COMM_WORLD, &num_proc));
  return num_proc;
}

//-----------------------------------------------------------------------------

int System::proc_num() {
  int proc_num = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &proc_num));
  return proc_num;
}

//-----------------------------------------------------------------------------

int System::compute_capability() {
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

  result = result && num_proc_ <= System::num_proc();

  if (num_proc_ > 1) {
    result = result && BuildHas::MPI;
  }

  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::CPU) {
    result = result && (TC::NONE == tc_ ||
                        (TC::FP32 == tc_ && BuildHas::CPUBLAS));
  }

  if (compute_method_ == ComputeMethod::GPU) {
    result = result && BuildHas::ACCEL && System::compute_capability() > 0;
  }

  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::GPU &&
      TC::FP16 == tc_) {
    result = result && BuildHas::CUDA && System::compute_capability() >= 700;
  }

  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::GPU &&
      TC::INT8 == tc_) {
    result = result && BuildHas::CUDA && System::compute_capability() >= 750;
  }

  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::GPU &&
      TC::FP32 == tc_) {
    result = result && (BuildHas::CUDA || BuildHas::HIP);
  }

  const bool tc_on = TC::FP16 == tc_ || TC::INT8 == tc_ ||
                     TC::FP32 == tc_;

  if (compute_method_ == ComputeMethod::GPU && (!is_metric_type_bitwise()
      || (is_metric_type_bitwise() && tc_on))) {
    result = result && BuildHas::MAGMA;
  }

  return result;
}

//-----------------------------------------------------------------------------

int Env::tc_eff() const {

  return tc_; // FIX
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
