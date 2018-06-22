//-----------------------------------------------------------------------------
/*!
 * \file   env.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Environment settings and general utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_env_hh_
#define _gm_env_hh_

#include "stddef.h"
#include "assert.h"
#include "float.h"
#include "stdio.h"  //FIX

#include "mpi.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"

#include "assertions.hh"
#include "types.hh"

//=============================================================================
/*---Environment struct declarations---*/

typedef struct {
  /*---Settings---*/
  int metric_type_;
  int num_way_;
  bool all2all_;
  int compute_method_;
  int num_stage;
  int stage_num;
  int num_phase;
  int phase_num;
  double ccc_param_;
  double ccc_multiplier_;
  bool are_ccc_params_default;
  bool sparse;
  /*---Counters---*/
  double time;
  double compares;
  double eltcompares;
  double ops_local;
  double ops;
  size_t cpu_mem;
  size_t cpu_mem_max;
  size_t gpu_mem;
  size_t gpu_mem_max;
  /*---MPI---*/
  bool make_comms_;
  MPI_Comm mpi_comm_base_;
  MPI_Comm mpi_comm_;
  MPI_Comm mpi_comm_repl_vector_;
  MPI_Comm mpi_comm_field_;
  /*---*/
  int num_proc_base_;
  int num_proc_;
  int num_proc_field_;
  int num_proc_repl_;
  int num_proc_vector_i_;
  int num_proc_repl_vector_;
  bool are_mpi_comms_initialized_;
  /*---*/
  int proc_num_base_;
  int proc_num_;
  int proc_num_field_;
  int proc_num_repl_;
  int proc_num_vector_i_;
  int proc_num_repl_vector_;
  /*---*/
  bool is_proc_active_;
  /*---CUDA---*/
  cudaStream_t stream_compute_;
  cudaStream_t stream_togpu_;
  cudaStream_t stream_fromgpu_;
  bool are_cuda_streams_initialized_;
  /*---HELPERS---*/
  bool do_reduce;
  bool need_2way; // does the 3-way calc require 2-way metrics
  /*---OTHER---*/
  const char* description;
  int tc;
  void* tc_buf_left;
  void* tc_buf_right;
  size_t tc_buf_size;
  int npvfl_step_max;
  int num_tc_steps;
  cublasHandle_t cublas_handle;
} GMEnv;

enum {
  GM_METRIC_TYPE_SORENSON = 0, // Not implemented
  GM_METRIC_TYPE_CZEK = 1,
  GM_METRIC_TYPE_CCC = 2,
  GM_NUM_METRIC_TYPE = 3
};

enum {
  GM_COMPUTE_METHOD_CPU = 0,
  GM_COMPUTE_METHOD_GPU = 1,
  GM_COMPUTE_METHOD_REF = 2,
  GM_NUM_COMPUTE_METHOD = 3
};

enum { GM_NUM_WAY_2 = 2, GM_NUM_WAY_3 = 3, GM_NUM_NUM_WAY = 2 };

//=============================================================================
/*---Null object---*/

GMEnv GMEnv_null(void);

//=============================================================================
/*---Initialize environment---*/

void gm_create_args(char* argstring, int* argc, char** argv);

void GMEnv_create(GMEnv* const env, MPI_Comm base_comm, int argc, char** argv,
                  const char* const description);

void GMEnv_create(GMEnv* const env, MPI_Comm base_comm,
                  const char* const options,
                  const char* const description);

void GMEnv_create_no_comms(GMEnv* const env, int argc, char** argv,
                           const char* const description,
                           int num_proc, int proc_num);

void GMEnv_create_no_comms(GMEnv* const env, const char* const options,
                           const char* const description,
                           int num_proc, int proc_num);

//=============================================================================
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* const env);

//=============================================================================
/*---Manage cuda streams---*/

void GMEnv_initialize_streams(GMEnv* const env);
void GMEnv_terminate_streams(GMEnv* const env);

//=============================================================================
/*---Accessors: general---*/

static int GMEnv_metric_type(const GMEnv* const env) {
  GMAssert(env);
  return env->metric_type_;
}

//-----------------------------------------------------------------------------

static double GMEnv_ccc_multiplier_default() {
  return ((double) 9) / ((double) 2);
}

//-----------------------------------------------------------------------------

static double GMEnv_ccc_param_default() {
  return ((double) 2) / ((double) 3);
}

//-----------------------------------------------------------------------------

static double GMEnv_ccc_param(const GMEnv* const env) {
  GMAssert(env);
  return env->ccc_param_;
}

//-----------------------------------------------------------------------------

static double GMEnv_ccc_multiplier(const GMEnv* const env) {
  GMAssert(env);
  return env->ccc_multiplier_;
}

//-----------------------------------------------------------------------------

static void GMEnv_ccc_param_set(double value, GMEnv* const env) {
  GMInsist(env);
  env->ccc_param_ = value;
  env->are_ccc_params_default =
     GMEnv_ccc_multiplier_default() == env->ccc_multiplier_ &&
     GMEnv_ccc_param_default() == env->ccc_param_;
}

//-----------------------------------------------------------------------------

static void GMEnv_ccc_multiplier_set(double value, GMEnv* const env) {
  GMInsist(env);
  GMInsist(value >= 0);
  env->ccc_multiplier_ = value;
  env->are_ccc_params_default =
     GMEnv_ccc_multiplier_default() == env->ccc_multiplier_ &&
     GMEnv_ccc_param_default() == env->ccc_param_;
}

//-----------------------------------------------------------------------------

static int GMEnv_num_way(const GMEnv* const env) {
  GMAssert(env);
  return env->num_way_;
}

//-----------------------------------------------------------------------------

static bool GMEnv_all2all(const GMEnv* const env) {
  GMAssert(env);
  return env->all2all_;
}

//-----------------------------------------------------------------------------

static int GMEnv_compute_method(const GMEnv* const env) {
  GMAssert(env);
  return env->compute_method_;
}

//-----------------------------------------------------------------------------

static MPI_Comm GMEnv_mpi_comm(const GMEnv* const env) {
  GMAssert(env);
  return env->mpi_comm_;
}

//-----------------------------------------------------------------------------

static MPI_Comm GMEnv_mpi_comm_repl_vector(const GMEnv* const env) {
  GMAssert(env);
  return env->mpi_comm_repl_vector_;
}

//-----------------------------------------------------------------------------

static MPI_Comm GMEnv_mpi_comm_field(const GMEnv* const env) {
  GMAssert(env);
  return env->mpi_comm_field_;
}

//-----------------------------------------------------------------------------

static int GMEnv_is_proc_active(const GMEnv* const env) {
  GMAssert(env);
  return env->is_proc_active_;
}

//-----------------------------------------------------------------------------

void GMEnv_set_compute_method(GMEnv* const env, int compute_method);
int GMEnv_data_type_vectors(const GMEnv* const env);
int GMEnv_data_type_metrics(const GMEnv* const env);

void GMEnv_set_num_proc(GMEnv* const env, int num_proc_vector_i,
                      int num_proc_repl, int num_proc_field);

cudaStream_t GMEnv_stream_compute(GMEnv* const env);
cudaStream_t GMEnv_stream_togpu(GMEnv* const env);
cudaStream_t GMEnv_stream_fromgpu(GMEnv* const env);

//=============================================================================
/*---Accessors: num proc---*/

static int GMEnv_num_block_vector(const GMEnv* const env) {
  GMAssert(env);
  return env->num_proc_vector_i_;
}

//-----------------------------------------------------------------------------

static int GMEnv_num_proc_vector_i(const GMEnv* const env) {
  GMAssert(env);
  return env->num_proc_vector_i_;
}

//-----------------------------------------------------------------------------

static int GMEnv_num_proc_repl(const GMEnv* const env) {
  GMAssert(env);
  return env->num_proc_repl_;
}

//-----------------------------------------------------------------------------

static int GMEnv_num_proc_repl_vector(const GMEnv* const env) {
  GMAssert(env);
  return env->num_proc_repl_vector_;
}

//-----------------------------------------------------------------------------

static int GMEnv_num_proc_field(const GMEnv* const env) {
  GMAssert(env);
  return env->num_proc_field_;
}

//-----------------------------------------------------------------------------

static int GMEnv_num_proc(const GMEnv* const env) {
  GMAssert(env);
  return env->num_proc_;
}

//=============================================================================
/*---Accessors: proc_num---*/

static int GMEnv_proc_num(const GMEnv* const env) {
  GMAssert(env);
  return env->proc_num_;
}

//-----------------------------------------------------------------------------

static int GMEnv_proc_num_vector_i(const GMEnv* const env) {
  GMAssert(env);
  return env->proc_num_vector_i_;
}

//-----------------------------------------------------------------------------

static int GMEnv_proc_num_repl(const GMEnv* const env) {
  GMAssert(env);
  return env->proc_num_repl_;
}

//-----------------------------------------------------------------------------

static int GMEnv_proc_num_field(const GMEnv* const env) {
  GMAssert(env);
  return env->proc_num_field_;
}

//=============================================================================
/*---Timer functions---*/

double GMEnv_get_time(const GMEnv* const env);
double GMEnv_get_synced_time(const GMEnv* const env);

//=============================================================================
/*---Math utility functions---*/

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

static int gm_popcount64(GMUInt64 x) {
  /*---Adapted from https://en.wikipedia.org/wiki/Hamming_weight---*/
  const GMUInt64 m1 = 0x5555555555555555;
  const GMUInt64 m2 = 0x3333333333333333;
  const GMUInt64 m4 = 0x0f0f0f0f0f0f0f0f;
  const GMUInt64 h01 = 0x0101010101010101;
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
/*---Misc.---*/

bool GMEnv_cuda_last_call_succeeded(const GMEnv* const env);

int gm_gpu_compute_capability();

#if 0
static bool gm_assertions_on() {
#ifdef GM_ASSERTIONS_ON
  return true;
#else
  return false;
#endif
}
#endif

//-----------------------------------------------------------------------------

void* gm_malloc(size_t n, GMEnv* env);
void gm_free(void* p, size_t n, GMEnv* env);

bool GMEnv_is_ppc64();

GMFloat* GMFloat_malloc(size_t n, GMEnv* env);
void GMFloat_free(GMFloat* p, size_t n, GMEnv* env);

void GMFloat_fill_nan(GMFloat* const a, size_t n);
void GMFloat_check(GMFloat* const a, size_t n);
//int GMFloat_mant_dig();

template<typename T> int gm_mant_dig();

//-----------------------------------------------------------------------------

void gm_tc_bufs_malloc(GMEnv* const env, int num_vector_local,
                       int num_field_local);
void gm_tc_bufs_free(GMEnv* const env);

size_t gm_num_vector_local_required(size_t num_vector_active,
                                    GMEnv* const env);

size_t gm_gemm_size_required(size_t size_requested, GMEnv* const env);

//-----------------------------------------------------------------------------

MPI_Datatype gm_mpi_type(const GMEnv* const env);

//-----------------------------------------------------------------------------

size_t gm_array_cksum(unsigned char* a, size_t n);

//=============================================================================

#endif /*---_gm_env_hh_---*/

//-----------------------------------------------------------------------------
