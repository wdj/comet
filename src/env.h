/*---------------------------------------------------------------------------*/
/*!
 * \file   env.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Environment settings and general utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================


=============================================================================*/

#ifndef _env_h_
#define _env_h_

#include <stddef.h>
#include <assert.h>

#include "mpi.h"
#include "cuda.h"
#include "cuda_runtime.h"

/*===========================================================================*/
/*---Assertions---*/

#define GMAssertAlways(condition) \
  (void)((condition) || (gm_assert(#condition, __FILE__, __LINE__), 0))

#define GMInsist(env, condition) \
  (void)((condition) || (gm_insist(env, #condition, __FILE__, __LINE__), 0))

#ifndef NDEBUG
#define GM_ASSERTIONS_ON
#define GMAssert(condition) GMAssertAlways(condition)
#else
#define GMAssert(condition)
#endif

#ifndef NDEBUG
/*---Fail compilation and (hopefully) give a filename/line number---*/
#define GMStaticAssert(condition) \
  {                               \
    int a[(condition) ? 1 : -1];  \
    (void) a;                     \
  }
#else
#define GM_StaticAssert(condition)
#endif

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

void gm_assert(const char* condition_string,
               const char* file,
               int line);

void gm_insist(const void* env,
               const char* condition_string,
               const char* file,
               int line);

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/
/*---Types---*/

#ifdef __cplusplus
extern "C" {
#endif

/*---Boolean type---*/

#ifdef __cplusplus
typedef bool _Bool;
#endif

enum { GM_BOOL_TRUE = (1 == 1), GM_BOOL_FALSE = (1 == 0) };

/*---Types for packed bits objects---*/

typedef double GMBits;

typedef unsigned long long int GMULInt;

/*---Default floating point type---*/

#ifdef FP_PRECISION_SINGLE
typedef float GMFloat;
enum { GM_MPI_FLOAT = MPI_FLOAT };
#endif

#ifdef FP_PRECISION_DOUBLE
typedef double GMFloat;
enum { GM_MPI_FLOAT = MPI_DOUBLE };
#endif

/*---Table for CCC---*/

typedef unsigned short int GMUSInt;

typedef struct { GMUSInt data[4][4]; } GMTable4X4;

/*---Type ids---*/

enum {
  GM_DATA_TYPE_FLOAT = 1,
  GM_DATA_TYPE_BIT = 2,
  GM_DATA_TYPE_TABLE4X4 = 3
};

/*---Dual CPU/GPU pointer---*/

typedef struct {
  void* __restrict__ h;
  void* __restrict__ d;
} GMMirroredPointer;

GMMirroredPointer GMMirroredPointer_null(void);

/*---Struct with checksum info---*/

enum { GM_CHECKSUM_SIZE = 3 };

typedef struct {
  size_t data[GM_CHECKSUM_SIZE];
} GMChecksum;

static _Bool gm_are_checksums_equal(GMChecksum c1, GMChecksum c2) {
  _Bool result = GM_BOOL_TRUE;
  int i = 0;
  for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    result = result && c1.data[i] == c2.data[i];
  }
  return result;
}

/*===========================================================================*/
/*---Environment struct declarations---*/

typedef struct {
  /*---Settings---*/
  int metric_type_;
  int num_way_;
  _Bool all2all_;
  int compute_method_;
  /*---MPI---*/
  int mpi_comm_;
  int num_proc_world_;
  int num_proc_;
  int proc_num_;
  _Bool is_proc_active_;
  /*---CUDA---*/
  cudaStream_t stream_compute_;
  cudaStream_t stream_togpu_;
  cudaStream_t stream_fromgpu_;
  _Bool are_cuda_streams_initialized_;
} GMEnv;

enum {
  GM_METRIC_TYPE_SORENSON = 0,
  GM_METRIC_TYPE_CZEKANOWSKI = 1,
  GM_METRIC_TYPE_CCC = 2,
  GM_NUM_METRIC_TYPE = 3
};

enum {
  GM_COMPUTE_METHOD_CPU = 0,
  GM_COMPUTE_METHOD_GPU = 1,
  GM_COMPUTE_METHOD_REF = 2,
  GM_NUM_COMPUTE_METHOD = 3
};

/*===========================================================================*/
/*---Null object---*/

GMEnv GMEnv_null(void);

/*===========================================================================*/
/*---Initialize environment---*/

void GMEnv_create(GMEnv* env);
void GMEnv_create_from_args(GMEnv* env, int argc, const char** argv);

/*===========================================================================*/
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* env);

/*===========================================================================*/
/*---Manage cuda streams---*/

void GMEnv_initialize_streams(GMEnv* env);
void GMEnv_terminate_streams(GMEnv* env);

/*===========================================================================*/
/*---Accessors---*/

static int Env_metric_type(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->metric_type_;
}

/*---------------------------------------------------------------------------*/

static int Env_num_way(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->num_way_;
}

/*---------------------------------------------------------------------------*/

static _Bool Env_all2all(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->all2all_;
}

/*---------------------------------------------------------------------------*/

static int Env_compute_method(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->compute_method_;
}

/*---------------------------------------------------------------------------*/

static int Env_mpi_comm(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->mpi_comm_;
}

/*---------------------------------------------------------------------------*/

static int Env_num_proc(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->num_proc_;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_num(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->proc_num_;
}

/*---------------------------------------------------------------------------*/

static int Env_is_proc_active(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->is_proc_active_;
}

/*---------------------------------------------------------------------------*/

void Env_set_compute_method(GMEnv* env, int compute_method);
int Env_data_type(const GMEnv* env);

void Env_set_num_proc(GMEnv* env, int num_proc);

cudaStream_t Env_stream_compute(const GMEnv* env);
cudaStream_t Env_stream_togpu(const GMEnv* env);
cudaStream_t Env_stream_fromgpu(const GMEnv* env);

/*===========================================================================*/
/*---Timer functions---*/

double GMEnv_get_time(GMEnv* env);
double GMEnv_get_synced_time(GMEnv* env);

/*===========================================================================*/
/*---Misc utility functions---*/

static int gm_min_i(const int i, const int j) {
  return i < j ? i : j;
}

/*---------------------------------------------------------------------------*/

static int gm_max_i(const int i, const int j) {
  return i > j ? i : j;
}

/*---------------------------------------------------------------------------*/

static int gm_floor_i(const int i, const int j) {
  GMAssert(j > 0);

  return i >= 0 ? i / j : (i - j + 1) / j;
}

/*---------------------------------------------------------------------------*/

static int gm_ceil_i(const int i, const int j) {
  GMAssert(j > 0);

  return -gm_floor_i(-i, j);
}

/*---------------------------------------------------------------------------*/

static size_t gm_ceil_i8(const size_t i, const size_t j) {
  GMAssert(i + 1 > 1);
  GMAssert(j + 1 > 1);

  return (i + j - 1) / j;
}

/*---------------------------------------------------------------------------*/

static size_t gm_randomize_max() {
  const size_t im = 714025;

  return im;
}

/*---------------------------------------------------------------------------*/

static size_t gm_randomize(size_t i) {
  const size_t im = 714025;
  const size_t ia = 4096;
  const size_t ic = 150889;

  return (i * ia + ic) % im;
}

/*---------------------------------------------------------------------------*/

static int gm_log2(size_t n) {
  if (n == 0) {
    return 0;
  }

  int result = 0;

  for (result = 0, n--; result <= 8*(int)sizeof(size_t); ++result) {
    if (n == 0) {
      break;
    }
    n >>= 1;
  }

  return result;
}

/*---------------------------------------------------------------------------*/

static size_t gm_nchoosek(int n, int k) {
  GMAssert(n >= 1);
  GMAssert(k >= 0 && k <= n);
  int i;
  size_t num = 1;
  size_t denom = 1;
  for (i = 0; i < k; ++i) {
    num *= (n - i);
    denom *= (i + 1);
  }
  return num / denom;
}

/*===========================================================================*/
/*---Misc.---*/

_Bool GMEnv_cuda_last_call_succeeded(GMEnv* env);

GMFloat* GMFloat_malloc(size_t n);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_env_h_---*/

/*---------------------------------------------------------------------------*/
