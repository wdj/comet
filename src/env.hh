//-----------------------------------------------------------------------------
/*!
 * \file   env.hh
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

#ifndef _COMET_ENV_HH_
#define _COMET_ENV_HH_

#include "cstdint"
#include "cstddef"
#include "cstring"
#include "assert.h"
#include "float.h"
#include "algorithm"
#include "vector"
#include "limits"
//#include "cstdio"  // for printf debugging

#include "mpi.h"

#if defined COMET_USE_CUDA
#  include "cuda.h"
#  include "cuda_runtime.h"
#elif defined COMET_USE_HIP
//#  include "hip/hip_runtime_api.h"
#  include "hip/hip_runtime.h"
//void hipdummy() {
//  void* dummy = (void*)s;
//  dummy++;
//}
#else
#  define __host__
#  define __device__
#  define __global__
#  define __forceinline__
static void dim3(size_t dim0, size_t dim1, size_t dim2) {}
#endif

#include "assertions.hh"
#include "types.hh"
#include "utils.hh"

//=============================================================================

#define COMET_MPI_SAFE_CALL(s) {int error_code = (s); \
                                COMET_INSIST(MPI_SUCCESS == error_code);}

//-----------------------------------------------------------------------------

#if defined(COMET_USE_CUDA)
# define COMET_LAUNCH_KERNEL(name, \
    numthreadblocks, threadblocksize, sharedmem, stream, ...) \
    name <<< numthreadblocks, threadblocksize, sharedmem, stream >>> \
      (__VA_ARGS__)
#elif defined(COMET_USE_HIP)
# define COMET_LAUNCH_KERNEL(name, \
    numthreadblocks, threadblocksize, sharedmem, stream, ...) \
    hipLaunchKernelGGL(name, \
      numthreadblocks, threadblocksize, sharedmem, stream, __VA_ARGS__)
#else
# define COMET_LAUNCH_KERNEL(name, \
    numthreadblocks, threadblocksize, sharedmem, stream, ...) \
    {numthreadblocks;} \
    {threadblocksize;} \
    {COMET_INSIST(false && "Attempt to launch kernel for non-accelerator build.");};
#endif

//-----------------------------------------------------------------------------

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Abstracted accelerator thread indexing/dimensions functions.

#if defined COMET_USE_CUDA && defined __CUDA_ARCH__
  __device__ static int threadIdx_x_() { return threadIdx.x; }

  __device__ static int blockIdx_x_() { return blockIdx.x; }
  __device__ static int blockIdx_y_() { return blockIdx.y; }
  __device__ static int blockIdx_z_() { return blockIdx.z; }

  __device__ static int blockDim_x_() { return blockDim.x; }

  __device__ static int gridDim_y_() { return gridDim.y; }
#elif defined COMET_USE_HIP && defined __HIPCC__
  __device__ static int threadIdx_x_() { return hipThreadIdx_x; }

  __device__ static int blockIdx_x_() { return hipBlockIdx_x; }
  __device__ static int blockIdx_y_() { return hipBlockIdx_y; }
  __device__ static int blockIdx_z_() { return hipBlockIdx_z; }

  __device__ static int blockDim_x_() { return hipBlockDim_x; }

  __device__ static int gridDim_y_() { return hipGridDim_y; }
#else
  __device__ static int threadIdx_x_() { return 0; }

  __device__ static int blockIdx_x_() { return 0; }
  __device__ static int blockIdx_y_() { return 0; }
  __device__ static int blockIdx_z_() { return 0; }

  __device__ static int blockDim_x_() { return 0; }

  __device__ static int gridDim_y_() { return 0; }
#endif

//-----------------------------------------------------------------------------
/// \brief Build options enums

struct BuildHas {
# ifdef COMET_USE_MPI
    enum {MPI = true};
# else
    enum {MPI = false};
# endif

# ifdef COMET_USE_CUDA
    enum {CUDA = true};
# else
    enum {CUDA = false};
# endif

# ifdef COMET_USE_HIP
    enum {HIP = true};
# else
    enum {HIP = false};
# endif

# ifdef COMET_USE_ACCEL
    enum {ACCEL = true};
# else
    enum {ACCEL = false};
# endif

# ifdef COMET_USE_MAGMA
    enum {MAGMA = true};
# else
    enum {MAGMA = false};
# endif

# ifdef COMET_USE_CPUBLAS
    enum {CPUBLAS = true};
# else
    enum {CPUBLAS = false};
# endif

# ifdef COMET_FP_PRECISION_DOUBLE
    enum {DOUBLE_PREC = true};
# else
    enum {DOUBLE_PREC = false};
# endif

# ifdef COMET_BUILD_TYPE_DEBUG
    enum {DEBUG = true};
# else
    enum {DEBUG = false};
# endif
};

//-----------------------------------------------------------------------------
/// \brief Helper class for system functions.

struct System {
  static int num_proc();
  static int proc_num();
  static bool is_proc_num_0() {return !proc_num();}
  static int compute_capability();
  static double time();
  static bool accel_last_call_succeeded();
};

//-----------------------------------------------------------------------------
/// \brief Helper class for metric type values.

struct MetricType {
  enum {INVALID = -1,
        //SORENSON = 0, // Not implemented
        CZEK = 1,
        CCC = 2,
        DUO = 3,
  };

  static const char* str(int metric_type) {
    return metric_type == CZEK ? "czekanowski" :
           metric_type == CCC  ? "ccc" :
           metric_type == DUO  ? "duo" :
                                 "(invalid)";
  }

  static int value(char* str) {
    return strcmp(str, "czekanowski") == 0 ? (int)CZEK :
           strcmp(str, "ccc") == 0         ? (int)CCC :
           strcmp(str, "duo") == 0         ? (int)DUO :
                                             (int)INVALID;
  }

#if 0
  static size_t metric_size(int metric_type) {
    COMET_INSIST(CZEK == metric_type || CCC == metric_type ||
                 DUO == metric_type);
    return   CZEK == metric_type    ? sizeof(GMFloat) :
             CCC  == metric_type    ? sizeof(GMTally2x2)
           /* DUO == metric_type */ : sizeof(GMTally4x2);
  }

  static size_t metric_item_coords_size() {return sizeof(MetricItemCoords_t);}
#endif
};

//-----------------------------------------------------------------------------
/// \brief Helper class to denote counted bits per element for metric type.

struct CountedBitsPerElement {
  enum {NONE = -1,
        DUO = 1,
        CCC = 2};
};

typedef CountedBitsPerElement CBPE;

//-----------------------------------------------------------------------------
/// \brief Helper class for compute method values.

struct ComputeMethod {
  enum {INVALID = -1,
        CPU = 0,
        GPU = 1,
        REF = 2,
  };

  static const char* str(int compute_method) {
    return compute_method == CPU ? "CPU" :
           compute_method == GPU ? "GPU" :
           compute_method == REF ? "REF" : "(invalid)";
  }

  static int value(char* str) {
    return strcmp(str, "CPU") == 0 ? (int)CPU :
           strcmp(str, "GPU") == 0 ? (int)GPU :
           strcmp(str, "REF") == 0 ? (int)REF : (int)INVALID;
  }

  static bool is_valid(int compute_method) {
    return compute_method >= 0 && compute_method < 3;
  }
};

//-----------------------------------------------------------------------------
/// \brief Helper class for num way values.

struct NumWay {
  enum {_2 = 2,
        _3 = 3};
  enum {MAX = 3};

  static bool is_valid(int num_way) {
    return num_way == _2 || num_way == _3;
  }
};

//-----------------------------------------------------------------------------
/// \brief Helper class for tc method values.

struct TC {
  enum {
    INVALID = -1,
    NO = 0,
    FP16 = 1,
    INT8 = 2,
    FP32 = 3,
    AUTO = 4,
    B1 = 5,
    NUM = 6
    //INT4 = 6,
  };

  static bool is_valid(int tc) {
    return tc >= 0 && tc <= NUM;
  }
};

//-----------------------------------------------------------------------------
/// \brief Helper class to manage metric item coordinates types.

struct CoordsType {
  enum {BY_METRIC = 1,
        BY_ENTRY = 2,
  };
};

//=============================================================================

class CEnv {
public:

  // NOTE: calling this class "Env" seems to cause strange runtime errors
  // on Lyra.

  //----------------------------------------
  // Constructor/destructor

  CEnv(MPI_Comm base_comm,
      int argc,
      char** argv,
      const char* const description = NULL);

  CEnv(MPI_Comm base_comm,
      const char* const options,
      const char* const description = NULL);

  CEnv(const char* const options,
      int num_proc = System::num_proc(),
      int proc_num = System::proc_num());

  ~CEnv();

  static void create_args(char* argstring, int* argc, char** argv);

  //----------------------------------------
  // CoMet Settings

  int metric_type() const {return metric_type_;}
  int all2all() const {return all2all_;}
  int compute_method() const {return compute_method_;}
  int num_way() const {return num_way_;}
  int pow2_num_way() const {return 1 << num_way_;}
  bool sparse() const {return sparse_;}
  int tc() const {return tc_;};
  int tc_eff() const {return tc_eff_;}
  int tc_eff_compute_() const;
  int num_tc_steps() const {return num_tc_steps_;};
  double threshold() const {return threshold_;}
  double metrics_shrink() const {return metrics_shrink_;}
  double ccc_multiplier() const {return ccc_multiplier_;}
  double duo_multiplier() const {return duo_multiplier_;}
  double ccc_param() const {return ccc_param_;}
  int num_stage() const {return num_stage_;}
  void num_stage(int value) {num_stage_ = value;}
  int stage_num() const {return stage_num_;}
  void stage_num(int value) {stage_num_ = value;}
  int num_phase() const {return num_phase_;}
  void num_phase(int value) {num_phase_ = value;}
  int phase_num() const {return phase_num_;}
  void phase_num(int value) {phase_num_ = value;}

  bool is_double_prec() const {return BuildHas::DOUBLE_PREC;}

  // CoMet Settings: threshold.

  static bool is_threshold(double t) {return t >= 0;}
  bool is_threshold() const {return CEnv::is_threshold(threshold_);}

  static double threshold_eff(double t) {
    return is_threshold(t) ? t : std::numeric_limits<double>::lowest();}
  double threshold_eff() const {return CEnv::threshold_eff(threshold_);}

  template<typename T>
  static __host__ __device__
  bool pass_threshold(T value, double threshold_eff) {
   return value > threshold_eff;
 }

  template<typename T>
  bool pass_threshold(T value) {
   //return CEnv::pass_threshold(value, threshold_eff());
   return CEnv::pass_threshold(value, threshold_eff_cache_);
 }

  // CoMet Settings: multiplier/param.

  static double ccc_multiplier_default() {return 9e0 / 2e0;}
  static double duo_multiplier_default() {return 4e0; }
  static double ccc_param_default() {return 2e0 / 3e0;}
  bool are_ccc_params_default() const {return are_ccc_params_default_;}

  // CoMet Settings: derived settings.

  bool is_compute_method_gpu() const {
    return ComputeMethod::GPU == compute_method_;
  }
  bool is_metric_type_bitwise() const {
    return MetricType::CCC == metric_type_ || MetricType::DUO == metric_type_;
  }
  int ijkE_max() const {return is_metric_type_bitwise() ? 2 : 1;}
  // Do we use TC package.
  static bool is_try_tc_(int tc_try) {return TC::NO != tc_try;}
  bool is_using_tc() const {return is_try_tc_(tc_eff());}
//  bool is_using_tc() const {
//    return is_try_tc_(tc_eff()) && is_metric_type_bitwise();
//  }
  // Do we use either MAGMA or TC.
  bool can_use_linalg_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    return ComputeMethod::GPU == compute_method_ ||
      (ComputeMethod::CPU == compute_method_ && is_metric_type_bitwise() &&
       is_try_tc_(tc_try));
  }
  bool is_using_linalg() const {return can_use_linalg_(tc_eff());}
//return ComputeMethod::GPU == compute_method_ ||
//    (ComputeMethod::CPU == compute_method_ && is_using_tc());
//  }
  // Do we form the X matrix in the TC package.
  bool form_matX_tc() const {return is_using_tc();}
  bool is_bitwise_3way_2step() const {return is_using_tc();}
  int num_step_2way_for_3way() const {
    return !(is_metric_type_bitwise() && is_using_linalg()) ? 1 :
           is_bitwise_3way_2step() ? 2 : 3;
  }
  bool does_3way_need_2way() const {
    return metric_type_ == MetricType::CZEK && is_using_linalg();
  }
  // CCC vs. DUO.
  int counted_bits_per_elt() const {
    COMET_INSIST(is_metric_type_bitwise());
    return MetricType::CCC == metric_type_ ? 2 : 1;
  }
  double ccc_duo_multiplier() {
    return counted_bits_per_elt() == CBPE::DUO ?
      duo_multiplier() : ccc_multiplier();
  }
  // Do we do final metrics calc and thresholding in TC package.
  bool can_threshold_tc_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    return is_try_tc_(tc_try) && sparse() && num_proc_field() == 1 && is_threshold()
      && !is_double_prec();
  }
  bool is_threshold_tc() const {return can_threshold_tc_(tc_eff());}
//    return is_using_tc() && sparse() && num_proc_field() == 1 && is_threshold()
//      && !is_double_prec();
//  }
  // Are 3-way metrics computed half block-plane at a time.
  bool is_vectors_halved() const {
    return NumWay::_3 == num_way() && is_threshold_tc() &&
          is_bitwise_3way_2step();
  }

  bool can_compress_enable_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    const bool try_compress = true;
    return
      try_compress &&
      is_try_tc_(tc_try) &&
      can_threshold_tc_(tc_try) &&
      num_way() == NumWay::_3 && // TODO: implement for 2-way
      BuildHas::ACCEL &&
      is_compute_method_gpu() &&
      !do_reduce();
  }

  bool is_compress_enabled() const {return can_compress_enable_(tc_eff());}
//    const bool try_compress = true;
//    return
//      try_compress &&
//      is_threshold_tc() &&
//      num_way() == NumWay::_3 && // TODO: implement for 2-way
//      BuildHas::ACCEL &&
//      is_compute_method_gpu() &&
//      !do_reduce();
//  }

  int metric_format() const {return is_threshold_tc() ?
    MetricFormat::SINGLE : MetricFormat::PACKED_DOUBLE;}

  // Accessors pertaining to metric sizes.

  int num_entries_per_metric() const {
    return MetricType::CZEK == metric_type_ ? 1 : pow2_num_way();
  }

  size_t metric_size() const {
    return MetricType::CZEK == metric_type_ ? sizeof(GMFloat) :
           NumWay::_2 == num_way_ ? sizeof(GMTally2x2) : sizeof(GMTally4x2);
  }

  // Shrink-related.

  bool can_shrink_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    const size_t storage_per_metric = metric_size() +
      sizeof(MetricItemCoords_t);
    const size_t storage_per_metric_shrink = metric_size() +
      sizeof(MetricItemCoords_t) * num_entries_per_metric();
    return can_threshold_tc_(tc_try) && can_compress_enable_(tc_try) &&
      is_try_tc_(tc_try) &&
      NumWay::_3 == num_way() && // TODO: implement 2-way
      storage_per_metric_shrink < metrics_shrink_ * storage_per_metric;
  }

  bool is_shrink() const {return can_shrink_(tc_eff());}
//    const size_t storage_per_metric = metric_size() +
//      sizeof(MetricItemCoords_t);
//    const size_t storage_per_metric_shrink = metric_size() +
//      sizeof(MetricItemCoords_t) * num_entries_per_metric();
//    return is_threshold_tc() && is_compress_enabled() &&
//      NumWay::_3 == num_way() && // TODO: implement 2-way
//      storage_per_metric_shrink < metrics_shrink_ * storage_per_metric;
//  }

  //size_t metric_entry_size() const {
  //  COMET_INSIST(metric_size() % num_entries_per_metric() == 0);
  //  return metric_size() / num_entries_per_metric();
  //}

  int num_entries_per_metric_item() const {
    return is_shrink() ? 1 : num_entries_per_metric();
  }

  int num_metric_items_per_metric() const {
    return is_shrink() ? num_entries_per_metric() : 1;
  }

  size_t metric_item_size() const {
    COMET_INSIST(metric_size() % num_metric_items_per_metric() == 0);
    return metric_size() / num_metric_items_per_metric();
  }

  size_t apply_shrink(size_t v) const {
    const double fuzz = 1.e-10;
    return is_shrink() ? (size_t)((v / metrics_shrink_) * (1.+fuzz)) : v;
  }

  int coords_type_compute_() const {
   return is_shrink() ? CoordsType::BY_ENTRY : CoordsType::BY_METRIC;
  }
  int coords_type() const {
   return coords_type_cache_;
   //return is_shrink() ? CoordsType::BY_ENTRY : CoordsType::BY_METRIC;
  }

  bool coords_type_by_metric() const {
   return CoordsType::BY_METRIC == coords_type_cache_;
  }

  // XOR GEMM-related.

  bool can_use_xor_(int tc_try) const {
    //return false;
    const bool try_use_xor_nonlinalg = true;
    const bool can_use_xor_nonlinalg =
      try_use_xor_nonlinalg &&
      NumWay::_2 == num_way_ &&
      ComputeMethod::CPU == compute_method_ &&
      !can_use_linalg_(tc_try);
    const bool is_not_using_threshold_detect =
      !is_threshold() || can_shrink_(tc_try) || can_threshold_tc_(tc_try);
    return
      MetricType::DUO == metric_type_ && // THIS LINE CURRENTLY REQUIRED
      sparse() && // THIS LINE CURRENTLY REQUIRED
      is_not_using_threshold_detect && // THIS LINE CURRENTLY REQUIRED
      //(can_use_xor_nonlinalg || is_try_tc_(tc_try));
      num_way() == NumWay::_2 && // TODO: implement for 3-way
      (can_use_xor_nonlinalg || TC::B1 == tc_try);
  }

  int is_using_xor() const {return can_use_xor_(tc_eff());}

  // Misc.

  int data_type_vectors() const;
  int data_type_metrics() const;
  int matrix_buf_elt_size() const {return MetricType::CZEK == metric_type_ ?
    sizeof(GMFloat) : sizeof(GMTally2x2);
  }

  MPI_Datatype metrics_mpi_type() const;

  // Is it possible to do run with this tc option on this platform / build.
  bool can_run(int tc) const;
  bool can_run() const {return can_run(tc_eff());};

  //----------------------------------------
  // Counters

  double ctime() const {return ctime_;}
  void ctime_inc(double t) {ctime_ += t;}
  double synced_time();
  size_t cpu_mem_local() const {return cpu_mem_local_;}
  size_t gpu_mem_local() const {return gpu_mem_local_;}
  void cpu_mem_local_inc(size_t n);
  void gpu_mem_local_inc(size_t n);
  void cpu_mem_local_dec(size_t n);
  void gpu_mem_local_dec(size_t n);
  size_t cpu_mem_max() const;
  size_t gpu_mem_max() const;
  void ops_local_inc(double n) {ops_local_ += n;}
  double ops() const;
  double entry_compares() const {return entry_compares_;}
  double metric_compares() const {return metric_compares_;}
  double vec_compares() const {return vec_compares_;}
  void entry_compares_inc(double n) {entry_compares_ += n;}
  void metric_compares_inc(double n) {metric_compares_ += n;}
  void vec_compares_inc(double n) {vec_compares_ += n;}
  size_t metric_entries() const {return metric_entries_;}
  size_t metric_entries_computed() const {return metric_entries_computed_;}
  void metric_entries_inc(size_t n) {metric_entries_ += n;}
  void metric_entries_computed_inc(size_t n) {metric_entries_computed_ += n;}
  double shrink_achieved() const {return shrink_achieved_;}
  void shrink_achieved_set(double v) {shrink_achieved_ = v;}

  //----------------------------------------
  // MPI comms

  MPI_Comm comm() const {return comm_;}
  MPI_Comm comm_repl_vector() const {return comm_repl_vector_;}
  MPI_Comm comm_field() const {return comm_field_;}

  //----------------------------------------
  // MPI proc counts

  int num_block_vector() const {return num_proc_vector_;}
  int num_proc_vector() const {return num_proc_vector_;}
  int num_proc_field() const {return num_proc_field_;}
  int num_proc_repl() const {return num_proc_repl_;}
  int num_proc_repl_vector() const {return num_proc_repl_*num_proc_vector_;}
  int num_proc() const {return num_proc_repl_vector()*num_proc_field_;}
  bool do_reduce() const {return num_proc_field_ > 1;}

  //----------------------------------------
  // MPI proc numbers

  int proc_num_vector() const {return proc_num_vector_;}
  int proc_num_field() const {return proc_num_field_;}
  int proc_num_repl() const {return proc_num_repl_;}
  int proc_num() const {return proc_num_;}
  bool is_proc_active() const {return is_proc_active_;}

  //----------------------------------------
  // Accelerator streams

# if defined COMET_USE_CUDA
    typedef cudaStream_t Stream_t;
# elif defined COMET_USE_HIP
    typedef hipStream_t Stream_t;
# else
    typedef int Stream_t;
# endif
  Stream_t stream_compute();
  Stream_t stream_togpu();
  Stream_t stream_fromgpu();
  void stream_synchronize(Stream_t stream) const;

//----------------------------------------
private:

  // Constructor/destructor
  void create_impl_(MPI_Comm base_comm, int argc,
    char** argv, const char* const description,
    bool make_comms, int num_proc, int proc_num);
  void set_defaults_();
  void parse_args_(int argc, char** argv);

  // CoMet Settings
  int metric_type_;
  bool all2all_;
  int compute_method_;
  int num_way_;
  bool sparse_;
  int tc_;
  int tc_eff_;
  int num_tc_steps_;
  double threshold_;
  double threshold_eff_cache_;
  double metrics_shrink_;
  double ccc_param_;
  double ccc_multiplier_;
  double duo_multiplier_;
  bool are_ccc_params_default_;
  void ccc_param_set_(double value) {
    ccc_param_ = value;
    are_ccc_params_default_ = ccc_multiplier_default() == ccc_multiplier_ &&
                              ccc_param_default() == ccc_param_;
  }
  void ccc_multiplier_set_(double value) {
    COMET_INSIST(value >= 0);
    ccc_multiplier_ = value;
    are_ccc_params_default_ = ccc_multiplier_default() == ccc_multiplier_ &&
                              ccc_param_default() == ccc_param_;
  }
  void duo_multiplier_set_(double value) {
    COMET_INSIST(value >= 0);
    duo_multiplier_ = value;
  }
  int num_stage_;
  int stage_num_;
  int num_phase_;
  int phase_num_;

  int coords_type_cache_;

  // Counters
  void accel_sync_() const;
  double ctime_;
  double ops_local_;
  size_t cpu_mem_local_;
  size_t gpu_mem_local_;
  size_t cpu_mem_max_local_;
  size_t gpu_mem_max_local_;
  double metric_compares_;
  double entry_compares_;
  double vec_compares_;
  size_t metric_entries_;
  size_t metric_entries_computed_;
  double shrink_achieved_;

  // MPI comms
  bool make_comms_;
  MPI_Comm comm_base_;
  MPI_Comm comm_;
  MPI_Comm comm_repl_vector_;
  MPI_Comm comm_field_;
  void comms_initialize_();
  void comms_terminate_();
  bool are_comms_initialized_;

  // MPI proc counts
  int num_proc_base_;
  int num_proc_;
  int num_proc_field_;
  int num_proc_repl_;
  int num_proc_vector_;
  int num_proc_repl_vector_;
  void set_num_proc_(int num_proc_vector, int num_proc_repl,
    int num_proc_field);

  // MPI proc numbers
  int proc_num_base_;
  int proc_num_;
  int proc_num_field_;
  int proc_num_repl_;
  int proc_num_vector_;
  int proc_num_repl_vector_;
  bool is_proc_active_;

  // Accelerator streams
  Stream_t stream_compute_;
  Stream_t stream_togpu_;
  Stream_t stream_fromgpu_;
  void streams_initialize_();
  void streams_terminate_();
  bool are_streams_initialized_;

  // Other
  const char* description_;

  // Disallowed methods.
  CEnv(const CEnv&);
  void operator=(const CEnv&);
};

//-----------------------------------------------------------------------------
/// \brief Null operation to avoid unused variable warning.

template<typename T>
static void no_unused_variable_warning(T& v) {}

//-----------------------------------------------------------------------------
/// \brief Templatized access to CCC or DUO front multiplier.

template<int COUNTED_BITS_PER_ELT = CBPE::NONE>
static double env_ccc_duo_multiplier(const CEnv& env);

template<>
double env_ccc_duo_multiplier<CBPE::CCC>(const CEnv& env) {
  return env.ccc_multiplier();
}

template<>
double env_ccc_duo_multiplier<CBPE::DUO>(const CEnv& env) {
  return env.duo_multiplier();
}

template<>
double env_ccc_duo_multiplier<CBPE::NONE>(const CEnv& env) {
  return env.counted_bits_per_elt() == CBPE::DUO ?
    env.duo_multiplier() : env.ccc_multiplier();
}

//=============================================================================
// Arrays and floating point

void* gm_malloc(size_t n, CEnv* env);
void gm_free(void* p, size_t n, CEnv* env);

GMFloat* GMFloat_malloc(size_t n, CEnv* env);
void GMFloat_free(GMFloat* p, size_t n, CEnv* env);

void GMFloat_fill_nan(GMFloat* const a, size_t n);
void GMFloat_check(GMFloat* const a, size_t n);

size_t gm_array_cksum(unsigned char* a, size_t n);

//=============================================================================


} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_ENV_HH_

//-----------------------------------------------------------------------------
