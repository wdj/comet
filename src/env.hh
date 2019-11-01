//-----------------------------------------------------------------------------
/*!
 * \file   env.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Basic environment - settings, MPI communicators, etc.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_env_hh_
#define _gm_env_hh_

#include "cstdint"
#include "cstddef"
#include "cstring"
#include "assert.h"
#include "float.h"
#include "algorithm"
#include "cstdio"  //FIX

#include "mpi.h"

#if defined USE_CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#elif defined USE_HIP
#include "hip/hip_runtime.h"
#endif

#include "assertions.hh"
#include "types.hh"

//=============================================================================

#define COMET_MPI_SAFE_CALL(s) {int error_code = (s); \
                                GMInsist(MPI_SUCCESS == error_code);}

//-----------------------------------------------------------------------------

namespace comet {

//-----------------------------------------------------------------------------
// Build options enums

struct BuildHas {
# ifdef USE_MPI
  enum {MPI = true};
# else
  enum {MPI = false};
# endif

# ifdef USE_CUDA
  enum {CUDA = true};
# else
  enum {CUDA = false};
# endif

# ifdef USE_HIP
  enum {HIP = true};
# else
  enum {HIP = false};
# endif

# ifdef USE_ACCEL
  enum {ACCEL = true};
# else
  enum {ACCEL = false};
# endif

# ifdef USE_MAGMA
  enum {MAGMA = true};
# else
  enum {MAGMA = false};
# endif

# ifdef USE_CPUBLAS
  enum {CPUBLAS = true};
# else
  enum {CPUBLAS = false};
# endif
};

//-----------------------------------------------------------------------------

struct System {
  static int num_proc();
  static int proc_num();
  static bool is_proc_num_0() {return !proc_num();}
  static int compute_capability();
  static double time();
};

//-----------------------------------------------------------------------------

struct MetricType {
  enum {INVALID = -1,
        //SORENSON = 0, // Not implemented
        CZEK = 1,
        CCC = 2,
        DUO = 3,
  };

  static const char* str(int metric_type) {
    return metric_type == CZEK ? "czekanowski" :
           metric_type == CCC ? "ccc" :
           metric_type == DUO ? "duo" : "(invalid)";
  }

  static int value(char* str) {
    return strcmp(str, "czekanowski") == 0 ? (int)CZEK :
           strcmp(str, "ccc") == 0 ? (int)CCC :
           strcmp(str, "duo") == 0 ? (int)DUO : (int)INVALID;
  }
};

//-----------------------------------------------------------------------------

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

struct NUM_WAY {
  enum {_2 = 2,
        _3 = 3};
  enum {MAX = 3};

  static bool is_valid(int num_way) {
    return num_way == _2 || num_way == _3;
  }
};

//-----------------------------------------------------------------------------

struct TC {
  enum {
    INVALID = -1,
    NONE = 0,
    //GM_TC_METHOD_BEST = 1,
    FP16 = 1,
    INT8 = 2,
    FP32 = 3,
    //GM_TC_METHOD_INT4 = 3,
    //GM_TC_METHOD_INT1 = 4,
  };

  static bool is_valid(int tc) {
    return tc >= 0 && tc < 4;
  }
};

//-----------------------------------------------------------------------------

class Env {
public:

  //----------------------------------------
  // CoMet Settings
  int num_stage;
  int stage_num;
  int num_phase;
  int phase_num;
  // Counters
  double compares;
  double eltcompares;
  double veccompares;
  double ops_local;
  double ops;
  size_t cpu_mem_local;
  size_t cpu_mem_max_local;
  size_t cpu_mem_max;
  size_t gpu_mem_local;
  size_t gpu_mem_max_local;
  size_t gpu_mem_max;

  //----------------------------------------
  // Constructor/destructor

  Env(MPI_Comm base_comm,
      int argc,
      char** argv,
      const char* const description = NULL);

  Env(MPI_Comm base_comm,
      const char* const options,
      const char* const description = NULL);

  Env(const char* const options,
      int num_proc = System::num_proc(),
      int proc_num = System::proc_num());

  ~Env();

  void set_defaults();
  void parse_args(int argc, char** argv);

  bool can_run() const;

  //----------------------------------------
  // CoMet Settings

  int metric_type() const {return metric_type_;}
  int all2all() const {return all2all_;}
  int compute_method() const {return compute_method_;}
  int num_way() const {return num_way_;}
  bool does_3way_need_2way() const {return metric_type_ == MetricType::CZEK;}
  bool is_metric_type_bitwise() const {
    return metric_type_==MetricType::CCC ||
           metric_type_==MetricType::DUO;
  }
  bool sparse() const {return sparse_;}
  int tc_eff() const;
  int num_tc_steps() const {return num_tc_steps_;};
  static double ccc_multiplier_default() {return ((double) 9) / ((double) 2);}
  static double duo_multiplier_default() {return (double) 4; }
  static double ccc_param_default() {return ((double) 2) / ((double) 3);}
  double ccc_param() const {return ccc_param_;}
  double ccc_multiplier() const {return ccc_multiplier_;}
  double duo_multiplier() const {return duo_multiplier_;}
  bool are_ccc_params_default() const {return are_ccc_params_default_;}

  // Counters
  double ctime() const {return ctime_;}
  void ctime_inc(double t) {ctime_ += t;}
  double synced_time();
  void cpu_mem_inc(size_t n) {
   cpu_mem_local += n;
   cpu_mem_max_local = std::max(cpu_mem_max_local, cpu_mem_local);
  }
  void gpu_mem_inc(size_t n) {
   gpu_mem_local += n;
   gpu_mem_max_local = std::max(gpu_mem_max_local, gpu_mem_local);
  }
  void cpu_mem_dec(size_t n) {cpu_mem_local -= n;}
  void gpu_mem_dec(size_t n) {gpu_mem_local -= n;}

  // MPI comms
  MPI_Comm comm() const {return comm_;}
  MPI_Comm comm_repl_vector() const {return comm_repl_vector_;}
  MPI_Comm comm_field() const {return comm_field_;}

  // MPI proc counts
  int num_block_vector() const {return num_proc_vector_;}
  int num_proc_vector() const {return num_proc_vector_;}
  int num_proc_field() const {return num_proc_field_;}
  int num_proc_repl() const {return num_proc_repl_;}
  int num_proc_repl_vector() const {return num_proc_repl_*num_proc_vector_;}
  int num_proc() const {return num_proc_repl_vector()*num_proc_field_;}
  bool do_reduce() const {return num_proc_field_ > 1;}

  // MPI proc numbers
  int proc_num_vector() const {return proc_num_vector_;}
  int proc_num_field() const {return proc_num_field_;}
  int proc_num_repl() const {return proc_num_repl_;}
  int proc_num() const {return proc_num_;}
  bool is_proc_active() const {return is_proc_active_;}

  // accelerator streams
# if defined USE_CUDA
  typedef cudaStream_t Stream_t;
# elif defined USE_HIP
  typedef hipStream_t Stream_t;
# else
  typedef int Stream_t;
# endif
  void stream_synchronize(Stream_t stream) const;
  Stream_t stream_compute();
  Stream_t stream_togpu();
  Stream_t stream_fromgpu();
  void streams_terminate_();

  bool are_streams_initialized_;

private:

  void create_impl_(MPI_Comm base_comm, int argc,
    char** argv, const char* const description,
    bool make_comms, int num_proc, int proc_num);

  // CoMet Settings
  int metric_type_;
  bool all2all_;
  int compute_method_;
  int num_way_;
  bool sparse_;
  int tc_;
  int num_tc_steps_;
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
    GMInsist(value >= 0);
    ccc_multiplier_ = value;
    are_ccc_params_default_ = ccc_multiplier_default() == ccc_multiplier_ &&
                              ccc_param_default() == ccc_param_;
  }
  void duo_multiplier_set_(double value) {
    GMInsist(value >= 0);
    duo_multiplier_ = value;
  }

  // Counters
  double ctime_;

  // MPI comms
  bool make_comms_;
  MPI_Comm comm_base_;
  MPI_Comm comm_;
  MPI_Comm comm_repl_vector_;
  MPI_Comm comm_field_;
  bool are_comms_initialized_;
  void comms_initialize_();
  void comms_terminate_();

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

  // accelerator streams
  void streams_initialize_();
  Stream_t stream_compute_;
  Stream_t stream_togpu_;
  Stream_t stream_fromgpu_;

  // OTHER
  const char* description_;

  // Disallowed methods.
  Env(  const Env&);
  void operator=(const Env&);
};

//----------

typedef Env GMEnv;

//=============================================================================
// Initialize environment

void gm_create_args(char* argstring, int* argc, char** argv);

//=============================================================================
// Accessors: general

//-----------------------------------------------------------------------------

int GMEnv_data_type_vectors(const GMEnv* const env);
int GMEnv_data_type_metrics(const GMEnv* const env);

//=============================================================================
// Math utility functions

// TODO: maybe put in separate file.

static int gm_min_i(const int i, const int j) {
  return i < j ? i : j;
}

//-----------------------------------------------------------------------------

static size_t gm_min_i8(const size_t i, const size_t j) {
  return i < j ? i : j;
}

//-----------------------------------------------------------------------------

static int gm_max_i(const int i, const int j) {
  return i > j ? i : j;
}

//-----------------------------------------------------------------------------

static size_t gm_max_i8(const size_t i, const size_t j) {
  return i > j ? i : j;
}

//-----------------------------------------------------------------------------

static int gm_floor_i(const int i, const int j) {
  GMAssert(j > 0);

  return i >= 0 ? i / j : (i - j + 1) / j;
}

//-----------------------------------------------------------------------------

static size_t gm_floor_i8(const size_t i, const size_t j) {
  GMAssert(i + 1 >= 1);
  GMAssert(j + 1 > 1);

  return i / j;
}

//-----------------------------------------------------------------------------

static int gm_ceil_i(const int i, const int j) {
  GMAssert(j > 0);

  return -gm_floor_i(-i, j);
}

//-----------------------------------------------------------------------------

static size_t gm_ceil_i8(const size_t i, const size_t j) {
  GMAssert(i + 1 >= 1);
  GMAssert(j + 1 > 1);

  return (i + j - 1) / j;
}

//-----------------------------------------------------------------------------

static int gm_mod_i(const int i, const int j) {
  GMAssert(j > 0);

  return i - j * gm_floor_i(i, j);
}

//-----------------------------------------------------------------------------

static size_t gm_randomize_max() {
  const size_t im = 714025;

  return im;
}

//-----------------------------------------------------------------------------

static size_t gm_randomize(size_t i) {
  const size_t im = 714025;
  const size_t ia = 4096;
  const size_t ic = 150889;

  return (i * ia + ic) % im;
}

//-----------------------------------------------------------------------------

static int gm_log2(size_t n) {
  if (n == 0) {
    return 0;
  }
 
  int result = 0; 
  for (result = 0, n--; result <= 8 * (int)sizeof(size_t); ++result) {
    if (n == 0) {
      break;
    }
    n >>= 1;
  }

  return result;
}

//-----------------------------------------------------------------------------

static size_t gm_nchoosek(int n, int k) {
  GMAssert(n >= 0);
  GMAssert(k >= 0 && k <= n);
  size_t numer = 1;
  size_t denom = 1;
  for (int i = 0; i < k; ++i) {
    numer *= (n - i);
    denom *= (i + 1);
  }
  return numer / denom;
}

//-----------------------------------------------------------------------------

static int gm_popcount64(uint64_t x) {
  // Adapted from https://en.wikipedia.org/wiki/Hamming_weight
  const uint64_t m1 = 0x5555555555555555;
  const uint64_t m2 = 0x3333333333333333;
  const uint64_t m4 = 0x0f0f0f0f0f0f0f0f;
  const uint64_t h01 = 0x0101010101010101;
  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;
  return (x * h01) >> 56;
}

//-----------------------------------------------------------------------------

static void GMFloat_sort_3(GMFloat* const __restrict__ min,
                           GMFloat* const __restrict__ mid,
                           GMFloat* const __restrict__ max,
                           const GMFloat* const __restrict__ a,
                           const GMFloat* const __restrict__ b,
                           const GMFloat* const __restrict__ c) {
  if (*a > *b) {
    if (*a > *c) {
      *max = *a;
      if (*b > *c) {
        *mid = *b;
        *min = *c;
      } else {
        *mid = *c;
        *min = *b;
      }
    } else {
      *mid = *a;
      *max = *c;
      *min = *b;
    }
  } else {
    if (*b > *c) {
      *max = *b;
      if (*a > *c) {
        *mid = *a;
        *min = *c;
      } else {
        *mid = *c;
        *min = *a;
      }
    } else {
      *mid = *b;
      *max = *c;
      *min = *a;
    }
  }
  GMAssert(*min <= *mid);
  GMAssert(*mid <= *max);
}

//=============================================================================
// Arrays and floating point

void* gm_malloc(size_t n, GMEnv* env);
void gm_free(void* p, size_t n, GMEnv* env);

bool GMEnv_is_ppc64();

GMFloat* GMFloat_malloc(size_t n, GMEnv* env);
void GMFloat_free(GMFloat* p, size_t n, GMEnv* env);

void GMFloat_fill_nan(GMFloat* const a, size_t n);
void GMFloat_check(GMFloat* const a, size_t n);

template<typename T> int gm_mant_dig();

MPI_Datatype gm_mpi_type(const GMEnv* const env);

size_t gm_array_cksum(unsigned char* a, size_t n);

//=============================================================================
// Misc.

bool GMEnv_accel_last_call_succeeded(const GMEnv* const env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _gm_env_hh_

//-----------------------------------------------------------------------------
