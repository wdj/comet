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
#include "errno.h"
//#include "cstdio"  // for printf debugging

#include "mpi.h"
# if COMET_USE_OPENMP
#   include "omp.h"
# endif

#if defined COMET_USE_CUDA
#  include "cuda.h"
#  include "cuda_runtime.h"
#  include "cublas_v2.h"
#  include "cusparse.h"
#elif defined COMET_USE_HIP
//#  include "hip/hip_runtime_api.h"
#  include "hip/hip_runtime.h"
#  include "rocblas.h"
#  include "rocsparse.h"
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
    {numthreadblocks; \
     threadblocksize; \
     COMET_INSIST(false && "Attempt to launch kernel for non-accelerator build.");};
#endif

//-----------------------------------------------------------------------------

namespace comet {

#if 0
//-----------------------------------------------------------------------------
/// \brief Helper class for system functions.

struct System {
  static int num_proc();
  static int proc_num();
  static bool is_proc_num_0() {return !proc_num();}
  static int compute_capability();
  static double time();
  static bool accel_last_call_succeeded();

  static bool is_in_parallel_region() {
#   if COMET_USE_OPENMP
      return omp_in_parallel();
#   else
      return false;
#   endif
  };
};
#endif

//-----------------------------------------------------------------------------
/// \brief Accelerator-related types.

# if defined COMET_USE_CUDA
    typedef cudaStream_t AccelStream_t;
    typedef cudaEvent_t AccelEvent_t;
//    typedef cublasHandle_t AccelBlasHandle_t;
//    typedef cusparseHandle_t AccelSparseHandle_t;
# elif defined COMET_USE_HIP
    typedef hipStream_t AccelStream_t;
    typedef hipEvent_t AccelEvent_t;
//    typedef rocblas_handle AccelBlasHandle_t;
//    typedef rocsparse_handle AccelSparseHandle_t;
# else
    typedef int AccelStream_t;
    typedef int AccelEvent_t;
//    typedef int AccelBlasHandle_t;
//    typedef int AccelSparseHandle_t;
# endif

//-----------------------------------------------------------------------------
/// \brief Struct to help manage (MAGMA) accelerator queues.

class QueueContainer {
public:
  bool is_initialized;
  void* magma_queue;
};

//-----------------------------------------------------------------------------
/// \brief Abstracted accelerator thread indexing/dimensions functions.

#if defined COMET_USE_CUDA && defined __CUDA_ARCH__
  __device__ static int threadIdx_x_() { return threadIdx.x; }

  __device__ static int blockIdx_x_() { return blockIdx.x; }
  __device__ static int blockIdx_y_() { return blockIdx.y; }
  __device__ static int blockIdx_z_() { return blockIdx.z; }

  __device__ static int blockDim_x_() { return blockDim.x; }

  __device__ static int gridDim_x_() { return gridDim.x; }
  __device__ static int gridDim_y_() { return gridDim.y; }
#elif defined COMET_USE_HIP && defined __HIPCC__
  __device__ static int threadIdx_x_() { return hipThreadIdx_x; }

  __device__ static int blockIdx_x_() { return hipBlockIdx_x; }
  __device__ static int blockIdx_y_() { return hipBlockIdx_y; }
  __device__ static int blockIdx_z_() { return hipBlockIdx_z; }

  __device__ static int blockDim_x_() { return hipBlockDim_x; }

  __device__ static int gridDim_x_() { return hipGridDim_x; }
  __device__ static int gridDim_y_() { return hipGridDim_y; }
#else
  __device__ static int threadIdx_x_() { return 0; }

  __device__ static int blockIdx_x_() { return 0; }
  __device__ static int blockIdx_y_() { return 0; }
  __device__ static int blockIdx_z_() { return 0; }

  __device__ static int blockDim_x_() { return 0; }

  __device__ static int gridDim_x_() { return 0; }
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

# ifdef COMET_USE_CUTLASS
    enum {CUTLASS = true};
# else
    enum {CUTLASS = false};
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

  enum {INT128 = sizeof(BasicTypes::BigUInt) == 128/8};

}; // BuildHas

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
    INT4 = 6,
    NUM = 7
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

//-----------------------------------------------------------------------------
// \brief Class to manage thresholds.

class Thresholds {

  enum {THRESHOLDS_COUNT_MAX = NumWay::MAX + 2}; // = 5

  typedef MetricFormatTraits<MetricFormat::SINGLE>::TypeIn MFTypeIn;

public:

  typedef double Elt_t;

  // Lightweight struct to pass thresholds to kernels.
  template<int NT>
  struct Data {
    Elt_t data[NT];
  };

  // Copy threshold values into lightweight struct.
  template<int NT> void copy(Data<NT>& data) {
    COMET_INSIST(NT == count_);
    for (int i=0; i<NT; ++i) {
      data.data[i] = thresholds_[i];
    }
  }

  // Set threshold values from string.
  void set(const char* input) {
    // TODO: implement for multi

    // Initialize (though technically not needed).
    for (int i=0; i<THRESHOLDS_COUNT_MAX; ++i) {
      thresholds_[i] = 0;
    }

    // parse string, store
    errno = 0;
    thresholds_[0] = (Elt_t)strtod(input, NULL);
    COMET_INSIST(0 == errno && "Invalid setting for threshold.");

    // obtain count_
    count_ = 1;



    // Check inputs.
    COMET_INSIST(1 == count_ || 4 == count_ || 5 == count_);
    if (count_ >= 4) {
      COMET_INSIST(is_active_(thresholds_[0]) == is_active_(thresholds_[1]) &&
                   is_active_(thresholds_[0]) == is_active_(thresholds_[2]) &&
                   is_active_(thresholds_[0]) == is_active_(thresholds_[3]));
    }
    if (count_ >= 5) {
      COMET_INSIST(is_active_(thresholds_[0]) == is_active_(thresholds_[4]));
    }
  }

  // TODO: cache this
  // Is this multi-threshlld case.
  bool is_multi() const {return 1 != count_;}

  // TODO: cache this
  // If multi threshold, is 2-way or 3-way implied.
  int num_way_multi() const {
    COMET_INSIST(is_multi());
    return 4 == count_ ? NumWay::_2 : NumWay::_3;
  }

  // TODO: cache this
  // Is thresholding active.
  // NOTE: all indiv. thresholds must be active for thresholding to be active. - CHECK
  bool is_active() const {
    if (!is_multi()) {
      return is_active_(thresholds_[0]);
    } else if (num_way_multi() == NumWay::_2) {
      return is_active_(thresholds_[0]) && is_active_(thresholds_[1]) &&
             is_active_(thresholds_[2]) && is_active_(thresholds_[3]);
    } else { // if (num_way_multi() == NumWay::_3)
      return is_active_(thresholds_[0]) && is_active_(thresholds_[1]) &&
             is_active_(thresholds_[2]) && is_active_(thresholds_[3]) &&
             is_active_(thresholds_[4]);
    }
  }

  // Is a value zero - for example from being thresholded out.
  template<typename T>
  static __host__ __device__
  bool is_zero(T value) {
   return value == (T)0;
  }

  // Does a single value pass a single threshold, active threshlds case.
  template<typename T>
  static __host__ __device__
  bool is_pass_active(T value, Elt_t threshold) {
   return value > threshold;
  }

  // Does a value pass the threshold, non-is_multi case.
  template<typename T>
  bool is_pass(T v) const {
    COMET_ASSERT(!is_multi());
    if (!is_active())
      return true;
    //return v > thresholds_[0]; // && v > 0;
    return is_pass_active(v, thresholds_[0]); // && v > 0;
  }

  // Does a value pass the threshold, for is_multi 2-way case.
  __host__ __device__
  static bool is_pass(Elt_t t0, Elt_t t1, Elt_t t2, Elt_t t3,
    // NOTE: t0=tLL, t1=tLH, t2=tHH, t3=TLLHH.
    int iE, int jE,
    MFTypeIn v00, MFTypeIn v01, MFTypeIn v10, MFTypeIn v11) {
    bool result = false;
    if (0 == iE) {
      if (0 == jE) {
        const bool is_pass_LL = v00 > t0; // && v00 > 0
        const bool is_pass_LLHH = v00 + v11 > t3;
        const bool is_pos_LL = v00 > 0;
        const bool is_pos_HH = v11 > 0;
        result = is_pass_LL || (is_pass_LLHH && is_pos_LL && is_pos_HH);
      } else {
        const bool is_pass_LH = v01 > t1; // && v01 > 0
        result = is_pass_LH;
      }
    } else {
      if (0 == jE) {
        const bool is_pass_HL = v10 > t1; // && v10 > 0
        result = is_pass_HL;
      } else {
        const bool is_pass_HH = v11 > t2; // && v11 > 0
        const bool is_pass_LLHH = v00 + v11 > t3;
        const bool is_pos_LL = v00 > 0;
        const bool is_pos_HH = v11 > 0;
        result = is_pass_HH || (is_pass_LLHH && is_pos_LL && is_pos_HH);
      }
    }
    return result;
  }

  // Does a value pass the threshold, for is_multi 2-way case.
  bool is_pass(Tally2x2<MetricFormat::SINGLE> ttable, int iE, int jE) const {
    COMET_ASSERT(0 == iE || 1 == iE);
    COMET_ASSERT(0 == jE || 1 == jE);
    if (!is_active())
      return true;
    typedef Tally2x2<MetricFormat::SINGLE> T;
    const auto v = T::get(ttable, iE, jE);
    bool result = false;
    if (!is_multi()) {
      result = is_pass_active(v, thresholds_[0]); // && v > 0;
    } else {
      COMET_ASSERT(num_way_multi() == NumWay::_2);
      result = is_pass(
        thresholds_[0], thresholds_[1], thresholds_[2], thresholds_[3], iE, jE,
        T::get(ttable, 0, 0), T::get(ttable, 0, 1),
        T::get(ttable, 1, 0), T::get(ttable, 1, 1));
    }
    COMET_ASSERT(!(is_active() && result && v <= 0));
    return result;
  }

  // Does a value pass the threshold, for is_multi 3-way case.
  __host__ __device__
  static bool is_pass(Elt_t t0, Elt_t t1, Elt_t t2, Elt_t t3, Elt_t t4,
    // NOTE: t0=tLLL, t1=tLLH, t2=tLHH, t3=tHHH t4=TLLLHHH.
    int iE, int jE, int kE,
    MFTypeIn v000, MFTypeIn v001, MFTypeIn v010, MFTypeIn v011,
    MFTypeIn v100, MFTypeIn v101, MFTypeIn v110, MFTypeIn v111) {
    bool result = false;
    if (0 == iE) {
      if (0 == jE) {
        if (0 == kE) {
          const bool is_pass_LLL = v000 > t0; // && v000 > 0
          const bool is_pass_LLLHHH = v000 + v111 > t4;
          const bool is_pos_LLL = v000 > 0;
          const bool is_pos_HHH = v111 > 0;
          result = is_pass_LLL || (is_pass_LLLHHH && is_pos_LLL && is_pos_HHH);
        } else {
          const bool is_pass_LLH = v001 > t1; // && v001 > 0
          result = is_pass_LLH;
        }
      } else {
        if (0 == kE) {
          const bool is_pass_LHL = v010 > t1; // && v010 > 0
          result = is_pass_LHL;
        } else {
          const bool is_pass_LHH = v011 > t2; // && v011 > 0
          result = is_pass_LHH;
        }
      }
    } else {
      if (0 == jE) {
        if (0 == kE) {
          const bool is_pass_HLL = v100 > t1; // && v100 > 0
          result = is_pass_HLL;
        } else {
          const bool is_pass_HLH = v101 > t2; // && v101 > 0
          result = is_pass_HLH;
        }
      } else {
        if (0 == kE) {
          const bool is_pass_HHL = v110 > t2; // && v110 > 0
          result = is_pass_HHL;
        } else {
          const bool is_pass_HHH = v111 > t3; // && v111 > 0
          const bool is_pass_LLLHHH = v000 + v111 > t4;
          const bool is_pos_LLL = v000 > 0;
          const bool is_pos_HHH = v111 > 0;
          result = is_pass_HHH || (is_pass_LLLHHH && is_pos_LLL && is_pos_HHH);
        }
      }
    }
    return result;
  }

  // Does a value pass the threshold, for is_multi 3-way case.
  template<typename TTable_t>
  bool is_pass(TTable_t ttable, int iE, int jE, int kE) const {
    COMET_ASSERT(0 == iE || 1 == iE);
    COMET_ASSERT(0 == jE || 1 == jE);
    COMET_ASSERT(0 == kE || 1 == kE);
    if (!is_active())
      return true;
    typedef TTable_t T;
    const auto v = T::get(ttable, iE, jE);
    bool result = false;
    if (!is_multi()) {
      result = is_pass_active(v, thresholds_[0]); // && v > 0;
    } else {
      COMET_ASSERT(num_way_multi() == NumWay::_3);
      result = is_pass(
        thresholds_[0], thresholds_[1], thresholds_[2], thresholds_[3],
        thresholds_[4], iE, jE, kE,
        T::get(ttable, 0, 0, 0), T::get(ttable, 0, 0, 1),
        T::get(ttable, 0, 1, 0), T::get(ttable, 0, 1, 1),
        T::get(ttable, 1, 0, 0), T::get(ttable, 1, 0, 1),
        T::get(ttable, 1, 1, 0), T::get(ttable, 1, 1, 1));
    }
    COMET_ASSERT(!(is_active() && result && v <= 0));
    return result;
  }

  //---------------------------------------------------------------------------

private:

  // Is a particular threshold value active.
  static bool is_active_(Elt_t t) {return t >= (Elt_t)0;}

  int count_;
  Elt_t thresholds_[THRESHOLDS_COUNT_MAX];
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
  //double threshold() const {return threshold_;}
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

  Thresholds& thresholds() {return thresholds_;}

  // Does a value pass the threhold.
  //template<typename T>
  //static __host__ __device__
  //bool pass_threshold_active(T value, double threshold) {
  //  return value > threshold;
  //}

  // Does a value pass the threhold.
  //template<typename T>
  //bool pass_threshold(T value) {
  //  return thresholds_.is_active() ? CEnv::pass_threshold_active(value, threshold_) : true;
  //}

  // CoMet Settings: multiplier/param.

  // Enable certain optimizations if using default settings.
  static double ccc_multiplier_default() {return 9e0 / 2e0;}
  static double duo_multiplier_default() {return 4e0; }
  static double ccc_param_default() {return 2e0 / 3e0;}
  bool are_ccc_params_default() const {return are_ccc_params_default_;}

  // CoMet Settings: derived settings.

  // Convenience function.
  bool is_compute_method_gpu() const {
    return ComputeMethod::GPU == compute_method_;
  }

  // Convenience function.
  bool is_metric_type_bitwise() const {
    return MetricType::CCC == metric_type_ || MetricType::DUO == metric_type_;
  }

  // Metric table dimension along each axis (=1 if scalar).
  int ijkE_max() const {return is_metric_type_bitwise() ? 2 : 1;}

  // Does tc value indicate a request to use the tc package.
  static bool is_try_tc_(int tc_try) {return TC::NO != tc_try;}

  // Is it requested AND possible to use the tc package.
  bool is_using_tc() const {return is_try_tc_(tc_eff());}

  // Can we use the linalg package for requested tc value
  // (enables MAGMA, CUTLASS, TC).  (= all GPU methods, or CPU via BLAS).
  bool can_use_linalg_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    return ComputeMethod::GPU == compute_method_ ||
      (ComputeMethod::CPU == compute_method_ && is_metric_type_bitwise() &&
       is_try_tc_(tc_try));
  }

  // Do we use the linalg packag.
  bool is_using_linalg() const {return can_use_linalg_(tc_eff());}

  // Do we form the X matrix (for 3-way) in the TC package.
  // (otherwise we form it on the CPU).
  bool form_matX_tc() const {return is_using_tc();}

  // Do 3-way bitwise methods use new faster method to
  // compute in 2 (sub)steps instead of 3.
  bool is_bitwise_3way_2step() const {return is_using_tc();}

  // How many (sub)steps actually needed for 3-way case.
  int num_step_2way_for_3way() const {
    return !(is_metric_type_bitwise() && is_using_linalg()) ? 1 :
           is_bitwise_3way_2step() ? 2 : 3;
  }

  // Does 3-way require prior computation of several 2-way results.
  bool does_3way_need_2way() const {
    return metric_type_ == MetricType::CZEK && is_using_linalg();
  }

  // Does each 2-bit value represent 2 (CCC) or 1 (DUO) bits counted.
  int counted_bits_per_elt() const {
    COMET_INSIST(is_metric_type_bitwise());
    return MetricType::CCC == metric_type_ ? 2 : 1;
  }

  // Convenience function.
  double ccc_duo_multiplier() {
    return counted_bits_per_elt() == CBPE::DUO ?
      duo_multiplier() : ccc_multiplier();
  }

  // Can we do metrics calc and thresholding in TC package for requested tc.
  bool can_threshold_tc_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    return is_try_tc_(tc_try) && sparse() && num_proc_field() == 1 &&
           thresholds_.is_active()
      && !is_double_prec();
  }

  // Do we do metrics calc and thresholding in TC package.
  bool is_threshold_tc() const {return can_threshold_tc_(tc_eff());}

  // Are 3-way metrics computed half block-plane at a time.
  // (each substep computes all 8 table entries for half of results).
  bool is_vectors_halved() const {
    return NumWay::_3 == num_way() && is_threshold_tc() &&
          is_bitwise_3way_2step();
  }

  // Can we enable metrics compression for requested tc value.
  // (enable compression means it will compress the data if compressible).
  // This is a lossless compression done right after thresholding.
  bool can_compress_enable_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    // Change this line to short circuit what follows, if desired.
    const bool try_compress = true;
    return
      try_compress &&
      is_try_tc_(tc_try) &&
      can_threshold_tc_(tc_try) &&
      //num_way() == NumWay::_3 && // TODO: implement for 2-way
      BuildHas::ACCEL &&
      is_compute_method_gpu() &&
      !do_reduce();
  }

  // Do we enable metrics compression.
  bool is_compress_enabled() const {return can_compress_enable_(tc_eff());}

  // Do we store 2X2 metrics tables a 4 ints or 2 packed doubles.
  int metric_format() const {return is_threshold_tc() ?
    MetricFormat::SINGLE : MetricFormat::PACKED_DOUBLE;}

  // Accessors pertaining to metric sizes.

  // Nomenclature:
  // "metric": a table of 4 or 8 entries (or trivially 1, if non-bitwise).
  // "entry": one of the entries of this table.
  // "metric item": an entry if using the "shrink" method, else a metric table.

  int num_entries_per_metric() const {
    return MetricType::CZEK == metric_type_ ? 1 : pow2_num_way();
  }

  // sizeof for a single table.
  size_t metric_size() const {
    return MetricType::CZEK == metric_type_ ? sizeof(GMFloat) :
           NumWay::_2 == num_way_ ? sizeof(GMTally2x2) : sizeof(GMTally4x2);
  }

  // Shrink-related.

  // Can we shrink metrics storage allocated on the CPU for the requested tc.
  // This only makes sense if the tc package compresses the metrics.

  bool can_shrink_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    const size_t storage_per_metric = metric_size() +
      sizeof(MetricItemCoords_t);
    const size_t storage_per_metric_shrink = metric_size() +
      sizeof(MetricItemCoords_t) * num_entries_per_metric();
    const bool is_shrinking_helpful =
      storage_per_metric_shrink < metrics_shrink_ * storage_per_metric;
    return can_threshold_tc_(tc_try) && can_compress_enable_(tc_try) &&
      is_try_tc_(tc_try) &&
      //NumWay::_3 == num_way() && // TODO: implement 2-way
      is_shrinking_helpful;
  }

  // Do we shrink metrics storage allocated on the CPU.
  bool is_shrink() const {return can_shrink_(tc_eff());}

  int num_entries_per_metric_item() const {
    return is_shrink() ? 1 : num_entries_per_metric();
  }

  int num_metric_items_per_metric() const {
    return is_shrink() ? num_entries_per_metric() : 1;
  }

  // sizeof for a metric item.
  size_t metric_item_size() const {
    COMET_INSIST(metric_size() % num_metric_items_per_metric() == 0);
    return metric_size() / num_metric_items_per_metric();
  }

  // Compute shrunk value based on shrinking factor.
  size_t apply_shrink(size_t v) const {
    const double fuzz = 1.e-10;
    return is_shrink() ? (size_t)((v / metrics_shrink_) * (1.+fuzz)) : v;
  }

  // Are metrics stored on CPU by metric or by metric table entry.
  int coords_type_compute_() const {
   return is_shrink() ? CoordsType::BY_ENTRY : CoordsType::BY_METRIC;
  }

  int coords_type() const {
   return coords_type_cache_;
   //return is_shrink() ? CoordsType::BY_ENTRY : CoordsType::BY_METRIC;
  }

  // Cache precomputed coords_type value as performance optimization.
  bool coords_type_by_metric() const {
   return CoordsType::BY_METRIC == coords_type_cache_;
  }

  // Can we use fast pre-check of all 4 (or 8) metric values for requested tc.
  bool can_use_threshold_detector(int tc_try) const {
      // Don't do if not doing thresholding at all.
      // Won't make sense if shrink because not stored as table.
      // Not done if thresholding is done early, in the tc package.
      return !(!thresholds_.is_active() ||
               can_shrink_(tc_try) ||
               can_threshold_tc_(tc_try));
  }

  // Do we use fast pre-check of all 4 (or 8) metric values.
  bool is_using_threshold_detector() const {
    return can_use_threshold_detector(tc_eff());
  }

  // XOR GEMM-related.

  // Can we use 1-bit xor gemm for requested tc value.
  bool can_use_xor_(int tc_try) const {
    COMET_INSIST(TC::AUTO != tc_try);
    // Change this line to short circuit what follows, if desired.
    const bool try_use_xor_nonlinalg = true;
    // Is xor method available outside of linalg package (i.e., on CPU).
    const bool can_use_xor_nonlinalg =
      try_use_xor_nonlinalg &&
      ComputeMethod::CPU == compute_method_ &&
      !can_use_linalg_(tc_try);
    // Change this line to short circuit what follows, if desired.
    const bool try_use_xor = true;
    const bool result =
      try_use_xor &&
      // 1-bit xor gemm currently only implemented for duo.
      MetricType::DUO == metric_type_ && 
      // 1-bit xor gemm currently only implemented for sparse.
      sparse() &&
      // Threshold detector currently doesn't support xor gemm.
      !can_use_threshold_detector(tc_try) &&
      // xor 3-way requires is_vectors_halved, thus also can_threshold_tc.
      (num_way() == NumWay::_2 ||
        // Implementation limits.
        // TODO: (possibly) implement more cases for 3-way
        (num_way() == NumWay::_3 && (can_threshold_tc_(tc_try) ||
                                     ComputeMethod::CPU == compute_method_))) &&
      // Have Ampere use standard GEMM.
      ((BuildHas::CUDA && compute_capability_cache_ <= 750) ||
       //((BuildHas::CUDA && compute_capability_cache_ <= 800) ||
       !(ComputeMethod::GPU == compute_method_ && TC::B1 == tc_try)) &&
      // Can only do if using 1-bit TC (check HW elsewhere) or if nonlinalg.
      (can_use_xor_nonlinalg || TC::B1 == tc_try);
    return result;
  }

  // Do we use 1-bit xor gemm.
  bool is_using_xor() const {return can_use_xor_(tc_eff());}

  // Do we use GPUDirect or similar.
  bool is_comm_gpu() const {
    const bool try_comm_gpu = false; // true; // change this to en/disable.
    return try_comm_gpu && ComputeMethod::GPU == compute_method_ &&
           num_way() == NumWay::_2 && !do_reduce();
  }

  // Do we use CUTLASS.

  bool can_use_cutlass(int tc_try) const {
    const bool use_b1 = tc_try == TC::B1 &&
      System::compute_capability() > 700 && MetricType::DUO == metric_type_;
    const bool use_int4 = tc_try == TC::INT4 &&
      System::compute_capability() > 700;
    //const bool use_int8 = false;
    const bool use_int8 = tc_try == TC::INT8 &&
      System::compute_capability() > 700;
    return BuildHas::CUDA && BuildHas::CUTLASS &&
      ComputeMethod::GPU == compute_method_ &&
      is_metric_type_bitwise() &&
      (use_b1 || use_int4 || use_int8);
  }

  bool is_using_cutlass() const {return can_use_cutlass(tc_eff());}

  bool can_use_cutlass_mockup(int tc_try) const {
    //return true;
    return (tc_try == TC::B1 || tc_try == TC::INT4) &&
      BuildHas::CUDA && BuildHas::CUTLASS &&
      System::compute_capability() <= 700;
  }

  bool is_using_cutlass_mockup() const {
    return can_use_cutlass_mockup(tc_eff());
  }

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
  double gemmtime_sum() const;
  void gemmtime_start();
  void gemmtime_end();
  void gemmtime_record();
  void ctime_inc(double t) {ctime_ += t;}
  void gemmtime_inc(double t) {gemmtime_ += t;}
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
  void ops_gemm_local_inc(double n) {ops_gemm_local_ += n;}
  double ops() const;
  double ops_gemm() const;

  double vec_compares() const {return vec_compares_;}
  void vec_compares_inc(double n) {vec_compares_ += n;}
  double vec_active_compares() const {return vec_active_compares_;}
  void vec_active_compares_inc(double n) {vec_active_compares_ += n;}
  double metric_compares() const {return metric_compares_;}
  void metric_compares_inc(double n) {metric_compares_ += n;}
  double metric_active_compares() const {return metric_active_compares_;}
  void metric_active_compares_inc(double n) {metric_active_compares_ += n;}
  double entry_compares() const {return entry_compares_;}
  void entry_compares_inc(double n) {entry_compares_ += n;}
  double entry_active_compares() const {return entry_active_compares_;}
  void entry_active_compares_inc(double n) {entry_active_compares_ += n;}

  size_t metric_entries() const {return metric_entries_;}
  void metric_entries_inc(size_t n) {metric_entries_ += n;}
  size_t metric_entries_computed() const {return metric_entries_computed_;}
  void metric_entries_computed_inc(size_t n) {metric_entries_computed_ += n;}
  double shrink_achieved() const {return shrink_achieved_;}
  void shrink_achieved_set(double v) {shrink_achieved_ = v;}

  //----------------------------------------
  // MPI comms

  MPI_Comm comm() const {return comm_;}
  MPI_Comm comm_repl_vector() const {return comm_repl_vector_;}
  MPI_Comm comm_field() const {return comm_field_;}

  bool is_comm_ring() const {return NumWay::_2 == num_way();}
  //bool is_comm_ring() const {return false;}

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

  int proc_num() const {return proc_num_;}
  bool is_proc_active() const {return is_proc_active_;}

  // This is FVR-order mapping of proc axes to MPI_COMM_WORLD ranks.

  int proc_num_field() const {return proc_num_ % num_proc_field_;}

  int proc_num_vector() const {
    return (proc_num_ / num_proc_field_) % num_proc_vector_;}

  int proc_num_repl() const {
    return (proc_num_ / num_proc_field_) / num_proc_vector_;}

  // This is a virtualized mapping, ordering is indep of map to HW (above).

  int proc_num_repl_vector(int proc_num_repl, int proc_num_vector) const {
    return proc_num_repl + num_proc_repl_ * proc_num_vector;}

  int proc_num_repl_vector() const {
    return proc_num_repl_vector(proc_num_repl(), proc_num_vector());}

  //----------------------------------------
  // Accelerator streams

  AccelStream_t stream_compute();
  AccelStream_t stream_togpu();
  AccelStream_t stream_fromgpu();
  AccelEvent_t start_event();
  AccelEvent_t end_event();
  bool is_event_active() const {return is_event_active_;}
  void is_event_active(bool value) {is_event_active_ = value;}
  void stream_synchronize(AccelStream_t stream) const;

  //----------------------------------------
  // Accelerator queues

  template<class MagmaQueue_t>
  void queues_initialize() {
    queue_compute_.magma_queue = (void*)new MagmaQueue_t(stream_compute(), *this);
    queue_togpu_.magma_queue = (void*)new MagmaQueue_t(stream_togpu(), *this);
    queue_fromgpu_.magma_queue = (void*)new MagmaQueue_t(stream_fromgpu(), *this);
    queue_compute_.is_initialized = true;
    queue_togpu_.is_initialized = true;
    queue_fromgpu_.is_initialized = true;
  }

  template<class MagmaQueue_t>
  void queues_terminate() {
    delete (MagmaQueue_t*)queue_compute_.magma_queue;
    delete (MagmaQueue_t*)queue_togpu_.magma_queue;
    delete (MagmaQueue_t*)queue_fromgpu_.magma_queue;
    queue_compute_.magma_queue = NULL;
    queue_togpu_.magma_queue = NULL;
    queue_fromgpu_.magma_queue = NULL;
    queue_compute_.is_initialized = false;
    queue_togpu_.is_initialized = false;
    queue_fromgpu_.is_initialized = false;
  }

  template<class MagmaQueue_t>
  typename MagmaQueue_t::queue_t queue_compute() {
    return ((MagmaQueue_t*)queue_compute_.magma_queue)->queue();
  }

  template<class MagmaQueue_t>
  typename MagmaQueue_t::queue_t queue_togpu() {
    return ((MagmaQueue_t*)queue_togpu_.magma_queue)->queue();
  }

  template<class MagmaQueue_t>
  typename MagmaQueue_t::queue_t queue_fromgpu() {
    return ((MagmaQueue_t*)queue_fromgpu_.magma_queue)->queue();
  }

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
  //double threshold_;
  Thresholds thresholds_;
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
  double gemmtime_;
  double ops_local_;
  double ops_gemm_local_;
  size_t cpu_mem_local_;
  size_t gpu_mem_local_;
  size_t cpu_mem_max_local_;
  size_t gpu_mem_max_local_;
  double vec_compares_;
  double metric_compares_;
  double entry_compares_;
  double vec_active_compares_;
  double metric_active_compares_;
  double entry_active_compares_;
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
  //int proc_num_field_;
  //int proc_num_repl_;
  //int proc_num_vector_;
  //int proc_num_repl_vector_;
  bool is_proc_active_;

  // Accelerator streams
  AccelStream_t stream_compute_;
  AccelStream_t stream_togpu_;
  AccelStream_t stream_fromgpu_;
  AccelEvent_t start_event_;
  AccelEvent_t end_event_;
  bool is_event_active_;
  void streams_initialize_();
  void streams_terminate_();
  bool are_streams_initialized_;

//FIX
public:
  // Accelerator queues
  QueueContainer queue_compute_;
  QueueContainer queue_togpu_;
  QueueContainer queue_fromgpu_;
  void queues_initialize_();
  void queues_terminate_();
  bool are_queues_initialized_;
private:

  // Other
  const char* description_;
  int compute_capability_cache_;

  // Disallowed methods.
  CEnv(const CEnv&);
  void operator=(const CEnv&);
};

//-----------------------------------------------------------------------------
/// \brief Null operation to avoid unused variable warning.

template<typename T>
static void no_unused_variable_warning(T& v) {}

template<typename T>
static void no_unused_type_warning() {}

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
