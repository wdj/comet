//-----------------------------------------------------------------------------
/*!
 * \file   env.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Basic environment - settings, MPI communicators, etc.
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
#include "cstring"
#include "math.h"
#include "limits"

#include "errno.h"
#include "sys/time.h"
#include "signal.h"

#include "mpi.h"

#if defined COMET_USE_CUDA
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
#if defined COMET_USE_CUDA
  cudaDeviceProp deviceProp;
  // Assume only one GPU per rank.
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, 0);
  const int compute_capability = error != cudaSuccess ? 0 :
    deviceProp.major * 100 + deviceProp.minor;
#elif defined COMET_USE_HIP
  hipDeviceProp_t deviceProp;
  hipGetDeviceProperties(&deviceProp, 0); // Assume only one GPU per rank.
  //const int compute_capability = deviceProp.major * 100 + deviceProp.minor;
  // This seems more stable than major/minor.
  const int compute_capability = deviceProp.gcnArch;
#else
  const int compute_capability = 0;
#endif
  COMET_INSIST(System::accel_last_call_succeeded());
  return compute_capability;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator did most recent call succeed.

bool System::accel_last_call_succeeded() {

#if defined COMET_USE_CUDA
  // NOTE: this read of the last error is a destructive read.
  cudaError_t error = cudaGetLastError();
  const bool result = error == cudaSuccess;

  if (!result) {
    fprintf(stderr, "CUDA error detected: %s\n", cudaGetErrorString(error));
  }

  return result;
#elif defined COMET_USE_HIP
  // NOTE: this read of the last error is (apparently) a destructive read.
  hipError_t error = hipGetLastError();
  const bool result = error == hipSuccess;

  if (!result) {
    fprintf(stderr, "HIP error detected: %s\n", hipGetErrorString(error));
  }

  return result;
#endif

  return true;
}

//=============================================================================
// Constructor/destructor

void CEnv::create_impl_(MPI_Comm base_comm, int argc,
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
/// \brief CEnv constructor using argc/argv.

CEnv::CEnv(MPI_Comm base_comm,
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
/// \brief CEnv constructor using options c_string.

CEnv::CEnv(MPI_Comm base_comm,
         const char* const options,
         const char* const description) {

  // Convert options string to args

  size_t len = strlen(options);
  char* argstring = (char*)malloc((len+1)*sizeof(char));
  char ** argv = (char**)malloc((len+1)*sizeof(char*));
  int argc = 0;
  strcpy(argstring, options);
  CEnv::create_args(argstring, &argc, argv);

  const int num_proc = 0;
  const int proc_num = 0;
  const bool make_comms = true;

  create_impl_(base_comm, argc, argv, description, make_comms, 
    num_proc, proc_num);

  free(argstring);
  free(argv);
}

//-----------------------------------------------------------------------------
/// \brief CEnv constructor using options c_string, don't alloc comms.

CEnv::CEnv(const char* const options, int num_proc, int proc_num) {

  // Convert options string to args

  size_t len = strlen(options);
  char* argstring = (char*)malloc((len+1)*sizeof(char));
  char ** argv = (char**)malloc((len+1)*sizeof(char*));
  int argc = 0;
  strcpy(argstring, options);
  CEnv::create_args(argstring, &argc, argv);
  
  create_impl_(MPI_COMM_WORLD, argc, argv, NULL, false, num_proc, proc_num);
  
  free(argstring);
  free(argv);
}

//-----------------------------------------------------------------------------
/// \brief CEnv destructor.

CEnv::~CEnv() {
  comms_terminate_();
//  handles_terminate_();
  streams_terminate_();
}

//-----------------------------------------------------------------------------
/// \brief Set scalar members of CEnv to default values.

void CEnv::set_defaults_() {

  // Set all to zero to start.
  memset((void*)this, 0, sizeof(*this));

  // CoMet Settings
  metric_type_ = MetricType::CZEK;
  num_way_ = NumWay::_2;
  all2all_ = false;
  compute_method_ = ComputeMethod::GPU;
  num_stage_ = 1;
  stage_num_ = 0;
  num_phase_ = 1;
  phase_num_ = 0;
  ccc_param_set_(CEnv::ccc_param_default());
  ccc_multiplier_set_(CEnv::ccc_multiplier_default());
  duo_multiplier_set_(CEnv::duo_multiplier_default());
  sparse_ = false;
  tc_ = TC::NO;
  tc_eff_ = tc_eff_compute_();
  num_tc_steps_ = 1;
  threshold_ = CEnv::threshold_eff(-1);
  threshold_eff_cache_ = threshold_;
  metrics_shrink_ = 1;
  coords_type_cache_ = 0;

  ctime_ = 0;
  gemmtime_ = 0;
  ops_local_ = 0;
  ops_gemm_local_ = 0;
  cpu_mem_local_ = 0;
  gpu_mem_local_ = 0;
  cpu_mem_max_local_ = 0;
  gpu_mem_max_local_ = 0;
  entry_compares_ = 0;
  metric_compares_ = 0;
  vec_compares_ = 0;
  metric_entries_ = 0;
  metric_entries_computed_ = 0;
  shrink_achieved_ = DBL_MAX;

  compute_capability_cache_ = System::compute_capability();

//FIX
  queue_compute_.is_initialized = false;
  queue_fromgpu_.is_initialized = false;
  queue_togpu_.is_initialized = false;
  queue_compute_.magma_queue = NULL;
  queue_fromgpu_.magma_queue = NULL;
  queue_togpu_.magma_queue = NULL;
}

//-----------------------------------------------------------------------------
/// \brief Parse command line arguments, set CEnv accordingly.

void CEnv::parse_args_(int argc, char** argv) {

  CEnv* env = this;

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
      COMET_INSIST_INTERFACE(env, 0 == errno && (num_way == NumWay::_2 ||
                                                 num_way == NumWay::_3)
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
      tc_eff_ = tc_eff_compute_();
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
    } else if (strcmp(argv[i], "--threshold") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for threshold.");
      errno = 0;
      const double threshold = strtod(argv[i], NULL);
      COMET_INSIST_INTERFACE(env, 0 == errno && "Invalid setting for threshold.");
      threshold_ = threshold;
      //--------------------
    } else if (strcmp(argv[i], "--metrics_shrink") == 0) {
      //--------------------
      ++i;
      COMET_INSIST_INTERFACE(env, i < argc && "Missing value for metrics_shrink.");
      errno = 0;
      const double metrics_shrink = strtod(argv[i], NULL);
      COMET_INSIST_INTERFACE(env, 0 == errno
                    && metrics_shrink > 0
                    && "Invalid setting for metrics_shrink.");
      metrics_shrink_ = metrics_shrink;
      //--------------------
    } // if/else
  }   // for i

  // Finalize settings.

  threshold_eff_cache_ = CEnv::threshold_eff(threshold_);
  coords_type_cache_ = coords_type_compute_();
}

//-----------------------------------------------------------------------------
/// \brief Utility to parse a string to construct argc/argv.

void CEnv::create_args(char* argstring, int* argc, char** argv) {
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
/// \brief Indicate the scalar type to use for vectors elements.

int CEnv::data_type_vectors() const {

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
/// \brief Indicate the scalar type to use for metrics elements.

int CEnv::data_type_metrics() const {

  switch (metric_type()) {
    case MetricType::CZEK:
      return GM_DATA_TYPE_FLOAT;
    case MetricType::CCC:
      return num_way() == NumWay::_2 ? GM_DATA_TYPE_TALLY2X2
                                     : GM_DATA_TYPE_TALLY4X2;
    case MetricType::DUO:
      return num_way() == NumWay::_2 ? GM_DATA_TYPE_TALLY2X2
                                     : GM_DATA_TYPE_TALLY4X2;
  }
  COMET_INSIST(false && "Invalid metric_type.");
  return 0;
}

//-----------------------------------------------------------------------------
/// \brief Determine whether requested run possible on this hardware and build.

bool CEnv::can_run(int tc_try) const {

  // Using "auto" here is inadequate to determine whether can run.
  COMET_INSIST(TC::AUTO != tc_try);

  bool result = true;

  // Reference methods never use linalg or tc packages.
  if (compute_method_ == ComputeMethod::REF)
    result = result && !is_try_tc_(tc_try) && !is_using_linalg();

  // MPI requirements.
  if (make_comms_) {
    result = result && num_proc_ <= System::num_proc();
    if (num_proc_ > 1)
      result = result && BuildHas::MPI;
  }

  // Can use tc package only for bitwise methods.
  result = result && !(!is_metric_type_bitwise() && is_try_tc_(tc_try));

  // CPU can use tc package only if have CPUBLAS and doing FP32 GEMM.
  if (is_metric_type_bitwise() && compute_method_ == ComputeMethod::CPU) {
    result = result && (!is_try_tc_(tc_try) ||
                        (TC::FP32 == tc_try && BuildHas::CPUBLAS));
  }

  // Use of GPU requires presence of GPU.
  if (is_compute_method_gpu()) {
    result = result && BuildHas::ACCEL && System::compute_capability() > 0;
  }

  // 3-way DUO requires either tc package or non-linalg code path.
  if (can_use_linalg_(tc_try) && !is_try_tc_(tc_try) &&
      num_way() == NumWay::_3 && metric_type() == MetricType::DUO) {
    result = false; // currently unimplemented
  }

  // If using GPU, must have MAGMA use if non-bitwise method OR
  // if bitwise method and not using tc package.
  if (is_compute_method_gpu() && (!is_metric_type_bitwise()
      || (is_metric_type_bitwise() && !is_try_tc_(tc_try)))) {
    result = result && BuildHas::MAGMA;
  }

  // Checks based on GPU compute capability.

  // NOTE: choices for AMD/HIP:
  // /opt/rocm/hip/bin/hipcc:@knownTargets = ('gfx701', 'gfx801', 'gfx802',
  // 'gfx803', 'gfx900', 'gfx906', 'gfx908', 'gfx1010', 'gfx1011', 'gfx1012');
  // NOTE: MI60 is 906
  const int cc_mi60 = 906; // 900;
  const int cc_mi100 = 908; // 900;
  const int cc_minone = 1000;

  if (is_metric_type_bitwise() && is_compute_method_gpu() &&
      TC::FP32 == tc_try) {
    // ISSUE: may need to adjust CUDA compute capability here.
    result = result && ((BuildHas::CUDA && System::compute_capability() >= 400)
                     || (BuildHas::HIP && System::compute_capability() >= cc_mi60));
  }

  if (is_metric_type_bitwise() && is_compute_method_gpu() &&
      TC::FP16 == tc_try) {
    // ISSUE: may need to adjust HIP compute capability here.
    result = result &&((BuildHas::CUDA && System::compute_capability() >= 700)
                  || (BuildHas::HIP && System::compute_capability() >= cc_minone));
                  //|| (BuildHas::HIP && System::compute_capability() >= cc_mi100));
  }

  if (is_metric_type_bitwise() && is_compute_method_gpu() &&
      TC::INT8 == tc_try) {
    // ISSUE: may need to adjust HIP compute capability here.
    // NOTE: pre-Turing can support INT8 but it is not necessarily fastest.
    result = result &&((BuildHas::CUDA && System::compute_capability() >= 750)
                  //|| (BuildHas::HIP && System::compute_capability() >= 1000));
                  //|| (BuildHas::HIP && System::compute_capability() >= 908));
                    || (BuildHas::HIP && System::compute_capability() >= cc_mi100));
  }

  if (is_metric_type_bitwise() && is_compute_method_gpu() && TC::B1 == tc_try) {
    // ISSUE: may need to adjust HIP compute capability here.
    // FIX: Temporary code below for testing xor mockup code on summit.
//  result = result && ((BuildHas::CUDA && System::compute_capability() >= 750)
    result = result && ((BuildHas::CUDA && System::compute_capability() >= 700)
                     || (BuildHas::HIP && System::compute_capability() >= cc_minone))
                    ; // && can_use_xor_(tc_try);
  }

  return result;
}

//-----------------------------------------------------------------------------
// \brief Select best tc value if AUTO specified.

int CEnv::tc_eff_compute_() const {

  if (TC::AUTO != tc_)
    return tc_;

  // NOTE: order is important here: fastest first.
  // TODO: move B1 to most favored status.
  //for (auto tc_try : {TC::B1, TC::INT8, TC::FP16, TC::FP32}) {
  for (auto tc_try : {TC::INT8, TC::FP16, TC::FP32, TC::B1}) {
    if (can_run(tc_try))
      return tc_try;
  }

  return TC::NO;
}

//-----------------------------------------------------------------------------
/// \brief MPI type to be used for metrics.

MPI_Datatype CEnv::metrics_mpi_type() const {

  if (metric_type() == MetricType::CZEK)
    return COMET_MPI_FLOAT;
  else if (metric_type() == MetricType::CCC)
    return MPI_DOUBLE_COMPLEX;
  else if (metric_type() == MetricType::DUO)
    return MPI_DOUBLE_COMPLEX;

  COMET_INSIST(false && "Invalid metric_type.");
  return MPI_DOUBLE_COMPLEX; // Should never get here.
}

//=============================================================================
// Counters

//-----------------------------------------------------------------------------
/// \brief Synchronize CPU to all GPU activity.

void CEnv::accel_sync_() const {

  if (! is_proc_active())
    return;

  if (!is_compute_method_gpu())
    return;

# if defined COMET_USE_CUDA
    cudaDeviceSynchronize();
# elif defined COMET_USE_HIP
    hipDeviceSynchronize();
# endif
  COMET_INSIST(System::accel_last_call_succeeded() &&
           "Failure in call to device synchronize.");
}

//-----------------------------------------------------------------------------
/// \brief Get system time, synced across all active processes and CPUs/GPUs.

double CEnv::synced_time() {

  if (! is_proc_active())
    return 0;

  accel_sync_();

  COMET_MPI_SAFE_CALL(MPI_Barrier(comm_));

  return System::time();
}

//-----------------------------------------------------------------------------
/// \brief Increment byte count of per-rank CPU memory used.

void CEnv::cpu_mem_local_inc(size_t n) {
 cpu_mem_local_ += n;
//printf("CPU %zu %zu\n", n, cpu_mem_local_); //FIX
 cpu_mem_max_local_ = utils::max(cpu_mem_max_local_, cpu_mem_local_);
}

//-----------------------------------------------------------------------------
/// \brief Increment byte count of per-rank GPU memory used.

void CEnv::gpu_mem_local_inc(size_t n) {
 gpu_mem_local_ += n;
//printf("                             GPU %zu %zu\n", n, gpu_mem_local_); //FIX
 gpu_mem_max_local_ = utils::max(gpu_mem_max_local_, gpu_mem_local_);
}

//-----------------------------------------------------------------------------
/// \brief Decrement byte count of per-rank CPU memory used.

void CEnv::cpu_mem_local_dec(size_t n) {
  COMET_INSIST(n <= cpu_mem_local_);
  cpu_mem_local_ -= n;
}

//-----------------------------------------------------------------------------
/// \brief Decrement byte count of per-rank GPU memory used.

void CEnv::gpu_mem_local_dec(size_t n) {
  COMET_INSIST(n <= gpu_mem_local_);
  gpu_mem_local_ -= n;
}

//-----------------------------------------------------------------------------
/// \brief Compute and return per-rank CPU memory (global) high water mark.

size_t CEnv::cpu_mem_max() const {
  size_t result = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&cpu_mem_max_local_, &result, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm()));
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Compute and return per-rank GPU memory (global) high water mark.

size_t CEnv::gpu_mem_max() const {
  size_t result = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&gpu_mem_max_local_, &result, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm()));
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Compute and return (global) number of operations performed.

double CEnv::ops() const {
  double result = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&ops_local_, &result, 1, MPI_DOUBLE,
    MPI_SUM, comm()));
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Compute and return (global) number of GEMM operations performed.

double CEnv::ops_gemm() const {
  double result = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&ops_gemm_local_, &result, 1, MPI_DOUBLE,
    MPI_SUM, comm()));
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Compute and return (global) sum GEMM timings.

double CEnv::gemmtime_sum() const {
  double result = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&gemmtime_, &result, 1, MPI_DOUBLE,
    MPI_SUM, comm()));
  return result;
}

//-----------------------------------------------------------------------------
/// \brief GEMM timer start.

void CEnv::gemmtime_start() {
# if defined COMET_USE_CUDA
    cudaEventRecord(start_event(), stream_compute());
# elif defined COMET_USE_HIP
    hipEventRecord(start_event(), stream_compute());
# endif
}

//-----------------------------------------------------------------------------
/// \brief GEMM timer end.

void CEnv::gemmtime_end() {
# if defined COMET_USE_CUDA
    cudaEventRecord(end_event(), stream_compute());
# elif defined COMET_USE_HIP
    hipEventRecord(end_event(), stream_compute());
# endif
  is_event_active(true);
}

//-----------------------------------------------------------------------------
/// \brief GEMM timer record.

void CEnv::gemmtime_record() {
  if (is_event_active()) {
#   if defined COMET_USE_CUDA
      cudaEventSynchronize(end_event());
      float time = 0;
      cudaEventElapsedTime(&time, start_event(), end_event());
      gemmtime_inc(time / 1000.);
#   elif defined COMET_USE_HIP
      hipEventSynchronize(end_event());
      float time = 0;
      hipEventElapsedTime(&time, start_event(), end_event());
      gemmtime_inc(time / 1000.);
#   endif
    is_event_active(false);
  }
}

//=============================================================================
// MPI comms

//-----------------------------------------------------------------------------
/// \brief Allocate MPI communicators needed for computation.

void CEnv::comms_initialize_() {

  if (are_comms_initialized_)
    return;

  if (! make_comms_)
    return;

  // Communicator for active procs to use for calculation.

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_, is_proc_active_, proc_num_,
    &comm_));

  // Communicator along repl / vector axis.

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_,
      is_proc_active_ ? proc_num_field() : num_proc_,
      is_proc_active_ ? proc_num_repl_vector() : proc_num_,
      &comm_repl_vector_));

  // Communicator along field axis.

  COMET_MPI_SAFE_CALL(MPI_Comm_split(comm_base_,
      is_proc_active_ ? proc_num_repl_vector() : num_proc_,
      is_proc_active_ ? proc_num_field() : proc_num_,
      &comm_field_));

  are_comms_initialized_ = true;
}

//-----------------------------------------------------------------------------
/// \brief Dellocate previously allocated MPI communicators.

void CEnv::comms_terminate_() {

  if (! are_comms_initialized_)
    return;

  // Destroy communicators.

  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_));
  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_repl_vector_));
  COMET_MPI_SAFE_CALL(MPI_Comm_free(&comm_field_));

  are_comms_initialized_ = false;
}

//=============================================================================
// Accelerator streams.

//-----------------------------------------------------------------------------
/// \brief Allocated accelerator streams needed for computation.

void CEnv::streams_initialize_() {

  if (are_streams_initialized_)
    return;

  if (!is_compute_method_gpu())
    return;

  for (AccelStream_t* const stream : {&stream_compute_, &stream_togpu_,
                                      &stream_fromgpu_}) {
#   if defined COMET_USE_CUDA
      cudaStreamCreate(stream);
#   elif defined COMET_USE_HIP
      hipStreamCreate(stream);
#   else
      if (stream) {}
#   endif
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to stream create.");
  }

# if defined COMET_USE_CUDA
    cudaEventCreate(&start_event_);
    cudaEventCreate(&end_event_);
# elif defined COMET_USE_HIP
    hipEventCreate(&start_event_);
    hipEventCreate(&end_event_);
# endif
  is_event_active_ = false;
  COMET_INSIST(System::accel_last_call_succeeded() &&
           "Failure in call to event create.");

  are_streams_initialized_ = true;
}

//-----------------------------------------------------------------------------
/// \brief Dellocate previously allocated accelerator streams.

void CEnv::streams_terminate_() {

  if (! are_streams_initialized_)
    return;

  for (const AccelStream_t stream : {stream_compute_, stream_togpu_,
                                     stream_fromgpu_}) {
#   if defined COMET_USE_CUDA
      cudaStreamDestroy(stream);
#   elif defined COMET_USE_HIP
      hipStreamDestroy(stream);
#   else
      if (stream) {}
#   endif
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to stream destroy.");
  }

# if defined COMET_USE_CUDA
    cudaEventDestroy(start_event_);
    cudaEventDestroy(end_event_);
# elif defined COMET_USE_HIP
    hipEventDestroy(start_event_);
    hipEventDestroy(end_event_);
# endif
  COMET_INSIST(System::accel_last_call_succeeded() &&
           "Failure in call to event destroy.");
  COMET_INSIST(!is_event_active_);

  are_streams_initialized_ = false;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator stream for kernel launches on accelerator.

AccelStream_t CEnv::stream_compute() {
  streams_initialize_(); // Lazy initialization.
  return stream_compute_;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator stream for transfers from CPU to GPU.

AccelStream_t CEnv::stream_togpu() {
  streams_initialize_(); // Lazy initialization.
  return stream_togpu_;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator stream for transfers to CPU from GPU.

AccelStream_t CEnv::stream_fromgpu() {
  streams_initialize_(); // Lazy initialization.
  return stream_fromgpu_;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator event for timing start.

AccelEvent_t CEnv::start_event() {
  streams_initialize_(); // Lazy initialization.
  return start_event_;
}

//-----------------------------------------------------------------------------
/// \brief Accelerator event for timing end.

AccelEvent_t CEnv::end_event() {
  streams_initialize_(); // Lazy initialization.
  return end_event_;
}

//-----------------------------------------------------------------------------
/// \brief CPU wait for accelerator stream to complete queued work.

void CEnv::stream_synchronize(AccelStream_t stream) const {

  if (!is_compute_method_gpu())
    return;

  COMET_INSIST(are_streams_initialized_);

# if defined COMET_USE_CUDA
    cudaStreamSynchronize(stream);
# elif defined COMET_USE_HIP
    hipStreamSynchronize(stream);
# endif
  COMET_INSIST(System::accel_last_call_succeeded() &&
           "Failure in call to stream synchronize.");
}

//=============================================================================
// MPI proc counts.

//-----------------------------------------------------------------------------
/// \brief Set up proc counts and proc numbers.

void CEnv::set_num_proc_(int num_proc_vector,
                      int num_proc_repl, int num_proc_field) {
  COMET_INSIST(num_proc_vector > 0);
  COMET_INSIST(num_proc_repl > 0);
  COMET_INSIST(num_proc_field > 0);

  if (make_comms_ && !BuildHas::MPI) {
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
  if (make_comms_) {
    COMET_INSIST(num_proc_ <= num_proc_base_ &&
                 "Number of procs requested exceeds number available.");
  }

  // Set proc nums

  proc_num_ = proc_num_base_;

  is_proc_active_ = proc_num_ < num_proc_;

#if 0
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
#endif
//    proc_num_field_ = proc_num_ % num_proc_field_;
//    proc_num_vector_ = (proc_num_ / num_proc_field_) % num_proc_vector_;
//    proc_num_repl_ = (proc_num_ / num_proc_field_) / num_proc_vector_;
#if 0
  }
#endif

//  proc_num_repl_vector_ = proc_num_repl_ + num_proc_repl_ * proc_num_vector_;

  // Destroy old communicators if necessary

  comms_terminate_();

  // Make new communicators

  comms_initialize_();
}

//=============================================================================
// Memory, arrays and floating point

void* gm_malloc(size_t n, CEnv* env) {
  COMET_INSIST(env);
  COMET_INSIST(n+1 >= 1);

  void* p = malloc(n);
  COMET_INSIST(p &&
           "Invalid pointer from malloc, possibly due to insufficient memory.");
  env->cpu_mem_local_inc(n);
  return p;
}

//-----------------------------------------------------------------------------

void gm_free(void* p, size_t n, CEnv* env) {
  COMET_INSIST(p && env);
  COMET_INSIST(n+1 >= 1);

  free(p);
  env->cpu_mem_local_dec(n);
}

//-----------------------------------------------------------------------------

GMFloat* GMFloat_malloc(size_t n, CEnv* env) {
  COMET_INSIST(env);
  COMET_INSIST(n+1 >= 1);

  GMFloat* p = (GMFloat*)gm_malloc(n * sizeof(GMFloat), env);
  GMFloat_fill_nan(p, n);
  return p;
}

//-----------------------------------------------------------------------------

void GMFloat_free(GMFloat* p, size_t n, CEnv* env) {
  COMET_INSIST(p && env);

  gm_free(p, n * sizeof(GMFloat), env);
}

//-----------------------------------------------------------------------------

void GMFloat_fill_nan(GMFloat* const a, size_t n) {
  COMET_INSIST(a);
  COMET_INSIST(n+1 >= 1);

# ifdef COMET_ASSERTIONS_ON
    GMFloat value = sqrt(-1);
    for (size_t i=0; i<n; ++i) {
      a[i] = value;
    }
# endif
}

//-----------------------------------------------------------------------------

void GMFloat_check(GMFloat* const a, size_t n) {
  COMET_INSIST(a);
  COMET_INSIST(n+1 >= 1);

# ifdef COMET_ASSERTIONS_ON
    bool no_nans_found = true;
    for (size_t i=0; i<n; ++i) {
      if (a[i] != a[i]) {
        no_nans_found = false;
      }
    }
    COMET_INSIST(no_nans_found);
# endif
}

//-----------------------------------------------------------------------------

size_t gm_array_cksum(unsigned char* a, size_t n) {
  COMET_INSIST(a);
  COMET_INSIST(n+1 >= 1);

  size_t result = 0;
  const size_t mask = (((size_t)1) << 32) - 1;

# pragma omp parallel for schedule(dynamic,1000) reduction(+:result)
  for (size_t i=0; i<n; ++i) {
    result += (a[i] * i) & mask;
  }

  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
