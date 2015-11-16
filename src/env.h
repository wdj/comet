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

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Types---*/

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

/*===========================================================================*/
/*---Environment struct declarations---*/

typedef struct {
  int metric_type;
  int num_way;
  _Bool all2all;
  int compute_method;
  int mpi_comm;
  int num_proc;
  int proc_num;
  cudaStream_t stream_compute;
  cudaStream_t stream_togpu;
  cudaStream_t stream_fromgpu;
  _Bool are_cuda_streams_initialized;
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

void GMEnv_create_from_args(GMEnv* env, int argc, char** argv);

void GMEnv_initialize_streams(GMEnv* env);

void GMEnv_finish_initializations(GMEnv* env);

/*===========================================================================*/
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* env);

/*===========================================================================*/
/*---Timer functions---*/

double GMEnv_get_time(GMEnv* env);
double GMEnv_get_synced_time(GMEnv* env);

/*===========================================================================*/
/*---Assertions---*/

#define GMAssert(v) assert(v)

#ifndef GMInsist
#define GMInsist(env, condition) \
  (void)((condition) || (gm_insist(env, #condition, __FILE__, __LINE__), 0))
#endif

void gm_insist(GMEnv* env,
               const char* condition_string,
               const char* file,
               int line);

#ifndef NDEBUG
#define GM_ASSERTIONS_ON
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

static size_t gm_randomize(size_t i) {
  const size_t im = 714025;
  const size_t ia = 4096;
  const size_t ic = 150889;

  return (i * ia + ic) % im;
}

/*---------------------------------------------------------------------------*/

static size_t gm_randomize_max() {
  const size_t im = 714025;

  return im;
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

int gm_data_type_from_metric_type(int metric_type, GMEnv* env);

_Bool GMEnv_cuda_last_call_succeeded(GMEnv* env);

GMFloat* GMFloat_malloc(size_t n);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_env_h_---*/

/*---------------------------------------------------------------------------*/
