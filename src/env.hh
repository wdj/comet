//-----------------------------------------------------------------------------
/*!
 * \file   env.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Basic environment - settings, MPI communicators, etc.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_env_hh_
#define _comet_env_hh_

#include "cstdint"
#include "cstddef"
#include "cstring"
#include "assert.h"
#include "float.h"
#include "algorithm"
#include "vector"
#include "cstdio"  //FIX

#include "mpi.h"

#if defined USE_CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#elif defined USE_HIP
#include "hip/hip_runtime.h"
#else
#define __host__
#define __device__
#endif

#include "assertions.hh"
#include "types.hh"
#include "utils.hh"

//=============================================================================

#define COMET_MPI_SAFE_CALL(s) {int error_code = (s); \
                                COMET_INSIST(MPI_SUCCESS == error_code);}

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

struct NUM_WAY {
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
    //GM_TC_METHOD_INT4 = 5,
    //GM_TC_METHOD_INT1 = 6,
    NUM = 5
  };

  static bool is_valid(int tc) {
    return tc >= 0 && tc <= NUM;
  }
};

//=============================================================================

class Env {
public:

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

  static void create_args(char* argstring, int* argc, char** argv);

  //----------------------------------------
  // CoMet Settings

  int metric_type() const {return metric_type_;}
  int all2all() const {return all2all_;}
  int compute_method() const {return compute_method_;}
  int num_way() const {return num_way_;}
  bool sparse() const {return sparse_;}
  int tc() const {return tc_;};
  int tc_eff() const;
  int num_tc_steps() const {return num_tc_steps_;};
  static double ccc_multiplier_default() {return ((double) 9) / ((double) 2);}
  static double duo_multiplier_default() {return (double) 4; }
  static double ccc_param_default() {return ((double) 2) / ((double) 3);}
  double ccc_param() const {return ccc_param_;}
  double ccc_multiplier() const {return ccc_multiplier_;}
  double duo_multiplier() const {return duo_multiplier_;}
  bool are_ccc_params_default() const {return are_ccc_params_default_;}
  int num_stage() const {return num_stage_;}
  void num_stage(int value) {num_stage_ = value;}
  int stage_num() const {return stage_num_;}
  void stage_num(int value) {stage_num_ = value;}
  int num_phase() const {return num_phase_;}
  void num_phase(int value) {num_phase_ = value;}
  int phase_num() const {return phase_num_;}
  void phase_num(int value) {phase_num_ = value;}

  bool is_metric_type_bitwise() const {
    return metric_type_==MetricType::CCC || metric_type_==MetricType::DUO;
  }
  bool is_using_tc() const {return tc_eff() != TC::NO &&
    is_metric_type_bitwise();
  }
  bool is_using_linalg() const {return ComputeMethod::GPU == compute_method_ ||
    (ComputeMethod::CPU == compute_method_ && is_using_tc());
  }
  int num_step_2way_for_3way() const {
    return is_metric_type_bitwise() && is_using_linalg() ? 3 : 1;
  }
  bool does_3way_need_2way() const {
    return metric_type_ == MetricType::CZEK && is_using_linalg();
  }

  int data_type_vectors() const;
  int data_type_metrics() const;
  MPI_Datatype metrics_mpi_type() const;
  bool form_matX_on_accel() const {return is_using_tc();}
  //bool form_matX_on_accel() const {return false;}

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
  double compares() const {return compares_;}
  double eltcompares() const {return eltcompares_;}
  double veccompares() const {return veccompares_;}
  void compares_inc(double n) {compares_ += n;}
  void eltcompares_inc(double n) {eltcompares_ += n;}
  void veccompares_inc(double n) {veccompares_ += n;}

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

# if defined USE_CUDA
    typedef cudaStream_t Stream_t;
# elif defined USE_HIP
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

  // Counters
  void accel_sync_() const;
  double ctime_;
  double ops_local_;
  size_t cpu_mem_local_;
  size_t gpu_mem_local_;
  size_t cpu_mem_max_local_;
  size_t gpu_mem_max_local_;
  double compares_;
  double eltcompares_;
  double veccompares_;

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
  Env(const Env&);
  void operator=(const Env&);
};

//----------

typedef Env GMEnv;

//=============================================================================
// Arrays and floating point

void* gm_malloc(size_t n, GMEnv* env);
void gm_free(void* p, size_t n, GMEnv* env);

GMFloat* GMFloat_malloc(size_t n, GMEnv* env);
void GMFloat_free(GMFloat* p, size_t n, GMEnv* env);

void GMFloat_fill_nan(GMFloat* const a, size_t n);
void GMFloat_check(GMFloat* const a, size_t n);

size_t gm_array_cksum(unsigned char* a, size_t n);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_env_hh_

//-----------------------------------------------------------------------------
