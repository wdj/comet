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

//=============================================================================
// System-related operations.

//-----------------------------------------------------------------------------
/// \brief System (wallclock) timer.

double System::time() {

  struct timeval tv;
  gettimeofday(&tv, NULL);
  double result = ((double)tv.tv_sec + (double)tv.tv_usec * 1.e-6);

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Number of processors (MPI ranks) available.

int System::num_proc() {
  int num_proc = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_size(MPI_COMM_WORLD, &num_proc));
  return num_proc;
}

//-----------------------------------------------------------------------------
/// \brief MPI rank in comm world.

int System::proc_num() {
  int proc_num = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &proc_num));
  return proc_num;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator compute capability.

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
/// \brief Accelerator did most recent call succeed.

bool System::accel_last_call_succeeded() {

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

//=============================================================================
// Constructor/destructor

void Env::create_impl_(MPI_Comm base_comm, int argc,
                       char** argv, const char* const description,
                       bool make_comms, int num_proc, int proc_num) {
  set_defaults_();

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
  // Other
  this->description_ = description; //FIX

  parse_args_(argc, argv);

  if (make_comms) {
    COMET_INSIST_INTERFACE(this, can_run() &&
                      "Invalid problem for this system and build.");
  }
}

//-----------------------------------------------------------------------------
/// \brief Env constructor using argc/argv.

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
/// \brief Env constructor using options c_string.

Env::Env(MPI_Comm base_comm,
         const char* const options,
         const char* const description) {

  // Convert options string to args

  size_t len = strlen(options);
  char* argstring = (char*)malloc((len+1)*sizeof(char));
  char ** argv = (char**)malloc((len+1)*sizeof(char*));
  int argc = 0;
  strcpy(argstring, options);
  Env::create_args(argstring, &argc, argv);

  const int num_proc = 0;
  const int proc_num = 0;
  const bool make_comms = true;

  create_impl_(base_comm, argc, argv, description, make_comms, 
    num_proc, proc_num);

  free(argstring);
  free(argv);
}

//-----------------------------------------------------------------------------
/// \brief Env constructor using options c_string, don't alloc comms.

Env::Env(const char* const options, int num_proc, int proc_num) {

  // Convert options string to args

  size_t len = strlen(options);
  char* argstring = (char*)malloc((len+1)*sizeof(char));
  char ** argv = (char**)malloc((len+1)*sizeof(char*));
  int argc = 0;
  strcpy(argstring, options);
  Env::create_args(argstring, &argc, argv);
  
  create_impl_(MPI_COMM_WORLD, argc, argv, NULL, false, num_proc, proc_num);
  
  free(argstring);
  free(argv);
}

//-----------------------------------------------------------------------------
/// \brief Env destructor.

Env::~Env() {
  comms_terminate_();
  streams_terminate_();
}

//-----------------------------------------------------------------------------
/// \brief Set scalar entries of Env to default values.

void Env::set_defaults_() {

  // Set all to zero to start.
  memset((void*)this, 0, sizeof(*this));

  // CoMet Settings
  metric_type_ = MetricType::CZEK;
  num_way_ = NUM_WAY::_2;
  all2all_ = false;
  compute_method_ = ComputeMethod::GPU;
  num_stage_ = 1;
  stage_num_ = 0;
  num_phase_ = 1;
  phase_num_ = 0;
  ccc_param_set_(Env::ccc_param_default());
  ccc_multiplier_set_(Env::ccc_multiplier_default());
  duo_multiplier_set_(Env::duo_multiplier_default());
  sparse_ = false;
  tc_ = TC::NO;
  num_tc_steps_ = 1;
}

//-----------------------------------------------------------------------------
/// \brief Parse command line arguments, set Env accordingly.

void Env::parse_args_(int argc, char** argv) {

  Env* env = this;

  for (int i = 1; i < argc; ++i) {
    if (false) {
      //--------------------
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for metric_type.");
      if (strcmp(argv[i], "czekanowski") == 0) {
        metric_type_ = MetricType::CZEK;
      } else if (strcmp(argv[i], "ccc") == 0) {
        metric_type_ = MetricType::CCC;
      } else if (strcmp(argv[i], "duo") == 0) {
        metric_type_ = MetricType::DUO;
      } else {
        COMET_INSIST_INTERFACE(env, false && "Invalid setting for metric_type.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--num_way") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_way.");
      errno = 0;
      const long num_way = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno && (num_way == NUM_WAY::_2 ||
                                            num_way == NUM_WAY::_3)
                               && "Invalid setting for num_way.");
      num_way_ = num_way;
      set_num_proc_(num_proc_vector_, num_proc_repl_, num_proc_field_);
      //--------------------
    } else if (strcmp(argv[i], "--all2all") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for all2all.");
      if (strcmp(argv[i], "yes") == 0) {
        all2all_ = true;
      } else if (strcmp(argv[i], "no") == 0) {
        all2all_ = false;
      } else {
        COMET_INSIST_INTERFACE(env, false && "Invalid setting for all2all.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for compute_method.");
      if (strcmp(argv[i], "CPU") == 0) {
        compute_method_ =  ComputeMethod::CPU;
      } else if (strcmp(argv[i], "GPU") == 0) {
        compute_method_ =  ComputeMethod::GPU;
      } else if (strcmp(argv[i], "REF") == 0) {
        compute_method_ =  ComputeMethod::REF;
      } else {
        COMET_INSIST_INTERFACE(env, false && "Invalid setting for compute_method.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      //--------------------
      ++i;
      errno = 0;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_proc_vector.");
      long num_proc_vector = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno
                    && (long)(int)num_proc_vector == num_proc_vector
                    && "Invalid setting for num_proc_vector.");
      set_num_proc_(num_proc_vector, num_proc_repl_, num_proc_field_);
      //--------------------
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      //--------------------
      ++i;
      errno = 0;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_proc_field.");
      long num_proc_field = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno
                    && (long)(int)num_proc_field == num_proc_field
                    && "Invalid setting for num_proc_field.");
      set_num_proc_(num_proc_vector_, num_proc_repl_, num_proc_field);
      //--------------------
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      //--------------------
      ++i;
      errno = 0;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_proc_repl.");
      long num_proc_repl = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno
                    && (long)(int)num_proc_repl == num_proc_repl
                    && "Invalid setting for num_proc_repl.");
      set_num_proc_(num_proc_vector_, num_proc_repl, num_proc_field_);
      //--------------------
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for ccc_param.");
      errno = 0;
      const double ccc_param = strtod(argv[i], NULL);
      COMET_INSIST_INTERFACE(env, 0 == errno && ccc_param >= 0
                               && "Invalid setting for ccc_param.");
      env->ccc_param_set_(ccc_param);
      //--------------------
    } else if (strcmp(argv[i], "--ccc_multiplier") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for ccc_multiplier.");
      errno = 0;
      const double ccc_multiplier = strtod(argv[i], NULL);
      COMET_INSIST_INTERFACE(env, 0 == errno && ccc_multiplier >= 0
                               && "Invalid setting for ccc_multiplier.");
      ccc_multiplier_set_(ccc_multiplier);
      //--------------------
    } else if (strcmp(argv[i], "--duo_multiplier") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for duo_multiplier.");
      errno = 0;
      const double duo_multiplier = strtod(argv[i], NULL);
      COMET_INSIST_INTERFACE(env, 0 == errno && duo_multiplier >= 0
                               && "Invalid setting for duo_multiplier.");
      duo_multiplier_set_(duo_multiplier);
      //--------------------
    } else if (strcmp(argv[i], "--sparse") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for sparse.");
      if (strcmp(argv[i], "yes") == 0) {
        sparse_ = true;
      } else if (strcmp(argv[i], "no") == 0) {
        sparse_ = false;
      } else {
        COMET_INSIST_INTERFACE(env, false && "Invalid setting for sparse.");
      }
      //--------------------
    } else if (strcmp(argv[i], "--tc") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for tc.");
      errno = 0;
      const long tc = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno
                    && (long)(int)tc == tc
                    && TC::is_valid(tc)
                    && "Invalid setting for tc.");
      tc_ = tc;
      //--------------------
    } else if (strcmp(argv[i], "--num_tc_steps") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for num_tc_steps.");
      errno = 0;
      const long num_tc_steps = strtol(argv[i], NULL, 10);
      COMET_INSIST_INTERFACE(env, 0 == errno
                    && (long)(int)num_tc_steps == num_tc_steps
                    && num_tc_steps >= 1
                    && "Invalid setting for tc.");
      num_tc_steps_ = num_tc_steps;
      //--------------------
    } // if/else
  }   // for i
}

//-----------------------------------------------------------------------------
/// \brief Utility to parse a string to construct argc/argv.

void Env::create_args(char* argstring, int* argc, char** argv) {
  const size_t len = strlen(argstring);

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
// CoMet settings

//-----------------------------------------------------------------------------
/// \brief Indicate the scalar type to use for vectors entries.

int Env::data_type_vectors() const {

  switch (metric_type()) {
    case MetricType::CZEK:
      return GM_DATA_TYPE_FLOAT;
    case MetricType::CCC:
      return GM_DATA_TYPE_BITS2;
    case MetricType::DUO:
      return GM_DATA_TYPE_BITS2;
  }
  COMET_INSIST(false && "Invalid metric_type.");
  return 0;
}

//-----------------------------------------------------------------------------
/// \brief Indicate the scalar type to use for metrics entries.

int Env::data_type_metrics() const {

  switch (metric_type()) {
    case MetricType::CZEK:
      return GM_DATA_TYPE_FLOAT;
    case MetricType::CCC:
      return num_way() == NUM_WAY::_2 ? GM_DATA_TYPE_TALLY2X2
                                      : GM_DATA_TYPE_TALLY4X2;
    case MetricType::DUO:
      return GM_DATA_TYPE_TALLY2X2; // 2-way only for now
  }
  COMET_INSIST(false && "Invalid metric_type.");
  return 0;
}

//-----------------------------------------------------------------------------
/// \brief Determine whether requested run possible on this hardware and build.

bool Env::can_run(int tc) const {

  COMET_INSIST(TC::AUTO != tc);

  bool result = true;

  result = result && num_proc_ <= System::num_proc();

  if (num_proc_ > 1) {
    result = result && BuildHas::MPI;
  }

  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::CPU) {
    result = result && (TC::NO == tc ||
                        (TC::FP32 == tc && BuildHas::CPUBLAS));
  }

  if (compute_method_ == ComputeMethod::GPU) {
    result = result && BuildHas::ACCEL && System::compute_capability() > 0;
  }

//TODO: adjust this for HIP case.
  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::GPU &&
      TC::FP16 == tc) {
    result = result && BuildHas::CUDA && System::compute_capability() >= 700;
  }

//TODO: adjust this for HIP case.
  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::GPU &&
      TC::INT8 == tc) {
    result = result && BuildHas::CUDA && System::compute_capability() >= 750;
  }

  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::GPU &&
      TC::FP32 == tc) {
    result = result && (BuildHas::CUDA || BuildHas::HIP);
  }

  if (compute_method_ == ComputeMethod::GPU && (!is_metric_type_bitwise()
      || (is_metric_type_bitwise() && TC::NO == tc))) {
    result = result && BuildHas::MAGMA;
  }

  return result;
}

//-----------------------------------------------------------------------------
// \brief Select best tc value if AUTO specified.

int Env::tc_eff() const {

  if (TC::AUTO != tc_)
    return tc_;

  // NOTE: order is important here: fastest first.
  for (auto tc : {TC::INT8, TC::FP16, TC::FP32}) {
    if (can_run(tc))
      return tc;
  }

  COMET_INSIST(false && "Suitable tc setting not found for this platform / build.");
  return 0;
}

//-----------------------------------------------------------------------------
/// \brief MPI type to be used for metrics.

MPI_Datatype Env::metrics_mpi_type() const {

  if (metric_type() == MetricType::CZEK) {
    return GM_MPI_FLOAT;
  } else if (metric_type() == MetricType::CCC) {
    return MPI_DOUBLE_COMPLEX;
  } else if (metric_type() == MetricType::DUO) {
    return MPI_DOUBLE_COMPLEX;
  }

  COMET_INSIST(false && "Invalid metric_type.");
  return MPI_DOUBLE_COMPLEX; // Should never get here.
}

//=============================================================================
// Counters

//-----------------------------------------------------------------------------
/// \brief Synchronize CPU to all GPU activity.

void Env::accel_sync_() const {

  if (! is_proc_active())
    return;

  if (compute_method() != ComputeMethod::GPU)
    return;

# if defined USE_CUDA
    cudaDeviceSynchronize();
# elif defined USE_HIP
    hipDeviceSynchronize();
#endif
  COMET_INSIST(System::accel_last_call_succeeded() &&
           "Failure in call to device synchronize.");
}

//-----------------------------------------------------------------------------
/// \brief Get system time, synced across all active processes and CPUs/GPUs.

double Env::synced_time() {

  if (! is_proc_active())
    return 0;

  accel_sync_();

  COMET_MPI_SAFE_CALL(MPI_Barrier(comm_));

  return System::time();
}

//-----------------------------------------------------------------------------
/// \brief Increment byte count of per-rank CPU memory used.

void Env::cpu_mem_local_inc(size_t n) {
 cpu_mem_local_ += n;
 cpu_mem_max_local_ = utils::max(cpu_mem_max_local_, cpu_mem_local_);
}

//-----------------------------------------------------------------------------
/// \brief Increment byte count of per-rank GPU memory used.

void Env::gpu_mem_local_inc(size_t n) {
 gpu_mem_local_ += n;
 gpu_mem_max_local_ = utils::max(gpu_mem_max_local_, gpu_mem_local_);
}

//-----------------------------------------------------------------------------
/// \brief Decrement byte count of per-rank CPU memory used.

void Env::cpu_mem_local_dec(size_t n) {
  COMET_INSIST(n <= cpu_mem_local_);
  cpu_mem_local_ -= n;
}

//-----------------------------------------------------------------------------
/// \brief Decrement byte count of per-rank GPU memory used.

void Env::gpu_mem_local_dec(size_t n) {
  COMET_INSIST(n <= gpu_mem_local_);
  gpu_mem_local_ -= n;
}

//-----------------------------------------------------------------------------
/// \brief Compute and return per-rank CPU memory (global) high water mark.

size_t Env::cpu_mem_max() const {
  size_t result = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&cpu_mem_max_local_, &result, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm()));
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Compute and return per-rank GPU memory (global) high water mark.

size_t Env::gpu_mem_max() const {
  size_t result = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&gpu_mem_max_local_, &result, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm()));
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Compute and return (global) number of operations performed.

double Env::ops() const {
  double result = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&ops_local_, &result, 1, MPI_DOUBLE,
    MPI_SUM, comm()));
  return result;
}

//=============================================================================
// MPI comms

//-----------------------------------------------------------------------------
/// \brief Allocate MPI communicators needed for computation.

void Env::comms_initialize_() {

  if (are_comms_initialized_)
    return;

  if (! make_comms_)
    return;

  // Communicator for active procs to use for calculation.

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_, is_proc_active_, proc_num_,
    &comm_));

  // Communicator along repl / vector axis.

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_,
      is_proc_active_ ? proc_num_field_ : num_proc_,
      is_proc_active_ ? proc_num_repl_vector_ : proc_num_,
      &comm_repl_vector_));

  // Communicator along field axis.

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_,
      is_proc_active_ ? proc_num_repl_vector_ : num_proc_,
      is_proc_active_ ? proc_num_field_ : proc_num_,
      &comm_field_));

  are_comms_initialized_ = true;
}

//-----------------------------------------------------------------------------
/// \brief Dellocate previously allocated MPI communicators.

void Env::comms_terminate_() {

  if (! are_comms_initialized_)
    return;

  // Destroy communicators.

  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_));
  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_repl_vector_));
  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_field_));

  are_comms_initialized_ = false;
}

//=============================================================================
// Accelerator streams

//-----------------------------------------------------------------------------
/// \brief Allocated accelerator streams needed for computation.

void Env::streams_initialize_() {

  if (are_streams_initialized_)
    return;

  if (compute_method() != ComputeMethod::GPU)
    return;

  for (Stream_t* const stream : {&stream_compute_, &stream_togpu_,
                                 &stream_fromgpu_}) {
#   if defined USE_CUDA
      cudaStreamCreate(stream);
#   elif defined USE_HIP
      hipStreamCreate(stream);
#   else
      if (stream) {}
#   endif
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to stream create.");
  }

  are_streams_initialized_ = true;
}

//-----------------------------------------------------------------------------
/// \brief Dellocate previously allocated accelerator streams.

void Env::streams_terminate_() {

  if (! are_streams_initialized_)
    return;

  for (const Stream_t stream : {stream_compute_, stream_togpu_,
                                stream_fromgpu_}) {
#   if defined USE_CUDA
      cudaStreamDestroy(stream);
#   elif defined USE_HIP
      hipStreamDestroy(stream);
#   else
      if (stream) {}
#   endif
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to stream destroy.");
  }

  are_streams_initialized_ = false;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator stream for kernel launches on accelerator.

Env::Stream_t Env::stream_compute() {
  streams_initialize_(); // Lazy initialization.
  return stream_compute_;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator stream for transfers from CPU to GPU.

Env::Stream_t Env::stream_togpu() {
  streams_initialize_(); // Lazy initialization.
  return stream_togpu_;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator stream for transfers to CPU from GPU.

Env::Stream_t Env::stream_fromgpu() {
  streams_initialize_(); // Lazy initialization.
  return stream_fromgpu_;
}

//-----------------------------------------------------------------------------
/// \brief CPU wait for accelerator stream to complete queued work.

void Env::stream_synchronize(Stream_t stream) const {

  if (compute_method() != ComputeMethod::GPU)
    return;

  COMET_INSIST(are_streams_initialized_);

# if defined USE_CUDA
    cudaStreamSynchronize(stream);
# elif defined USE_HIP
    hipStreamSynchronize(stream);
# endif
  COMET_INSIST(System::accel_last_call_succeeded() &&
           "Failure in call to stream synchronize.");
}

//=============================================================================
// MPI proc counts.

//-----------------------------------------------------------------------------
/// \brief Set up proc counts and proc numbers.

void Env::set_num_proc_(int num_proc_vector,
                      int num_proc_repl, int num_proc_field) {
  COMET_INSIST(num_proc_vector > 0);
  COMET_INSIST(num_proc_repl > 0);
  COMET_INSIST(num_proc_field > 0);

  if (!BuildHas::MPI) {
    COMET_INSIST(num_proc_vector == 1);
    COMET_INSIST(num_proc_repl == 1);
    COMET_INSIST(num_proc_field == 1);
  }

  COMET_INSIST(num_proc_base_);
  //COMET_INSIST(proc_num_base_ is initialized);

  // Set proc counts

  num_proc_vector_ = num_proc_vector;
  num_proc_repl_ = num_proc_repl;
  num_proc_field_ = num_proc_field;

  num_proc_repl_vector_ = num_proc_repl_ * num_proc_vector_;

  num_proc_ = num_proc_repl_vector_ * num_proc_field;
  COMET_INSIST(num_proc_ <= num_proc_base_ &&
           "Number of procs requested exceeds number available.");

  // Set proc nums

  proc_num_ = proc_num_base_;

  is_proc_active_ = proc_num_ < num_proc_;

  // Choose axis ordering for proc axes.

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

//=============================================================================
// Memory, arrays and floating point

void* gm_malloc(size_t n, GMEnv* env) {
  COMET_INSIST(env);
  void* p = malloc(n);
  COMET_INSIST(p &&
           "Invalid pointer from malloc, possibly due to insufficient memory.");
  env->cpu_mem_local_inc(n);
  return p;
}

//-----------------------------------------------------------------------------

void gm_free(void* p, size_t n, GMEnv* env) {
  COMET_INSIST(p && env);
  free(p);
  env->cpu_mem_local_dec(n);
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
  COMET_INSIST(env);
  GMFloat* p = (GMFloat*)gm_malloc(n * sizeof(GMFloat), env);
  GMFloat_fill_nan(p, n);
  return p;
}

//-----------------------------------------------------------------------------

void GMFloat_free(GMFloat* p, size_t n, GMEnv* env) {
  COMET_INSIST(p && env);
  gm_free(p, n * sizeof(GMFloat), env);
}

//-----------------------------------------------------------------------------

void GMFloat_fill_nan(GMFloat* const a, size_t n) {
  COMET_INSIST(a);
  COMET_INSIST(n+1 >= 1);
#ifdef COMET_ASSERTIONS_ON
  GMFloat value = sqrt(-1);
  size_t i = 0;
  for (i=0; i<n; ++i) {
    a[i] = value;
  }
#endif
}

//-----------------------------------------------------------------------------

void GMFloat_check(GMFloat* const a, size_t n) {
  COMET_INSIST(a);
  COMET_INSIST(n+1 >= 1);
#ifdef COMET_ASSERTIONS_ON
  bool no_nans_found = true;
  size_t i = 0;
  for (i=0; i<n; ++i) {
    if (a[i] != a[i]) {
      no_nans_found = false;
    }
  }
  COMET_INSIST(no_nans_found);
#endif
}

//-----------------------------------------------------------------------------

size_t gm_array_cksum(unsigned char* a, size_t n) {
  COMET_INSIST(a);

  size_t result = 0;

  const size_t mask = (((size_t)1) << 32) - 1;

#pragma omp parallel for schedule(dynamic,1000) reduction(+:result)
  for (size_t i=0; i<n; ++i) {
    result += (a[i] * i) & mask;
  }

  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
