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
#include "stdio.h"  //FIX

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
#define GMStaticAssert(condition)
#endif

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

void gm_assert(const char* condition_string, const char* file, int line);

void gm_insist(const void* env,
               const char* condition_string,
               const char* file,
               int line);

/*===========================================================================*/
/*---Types: general---*/

/*---Boolean type---*/

#ifdef __cplusplus
typedef bool _Bool;
#endif

enum { GM_BOOL_TRUE = (1 == 1), GM_BOOL_FALSE = (1 == 0) };

/*---Floating point of explicit precision---*/

typedef double GMFp64;

static void gm_check_fp64() {
  GMStaticAssert(sizeof(GMFp64) == 64/8);
}

/*===========================================================================*/
/*---Types for Sorenson metric---*/

//---TODO: complete the design for orenson metric.

/*---For Vectors: single 1-bit value: use unsigned int---*/
typedef unsigned int GMBits1;

typedef unsigned long long int GMBits1x64;

typedef unsigned long long int GMUInt64;

/*---Return null value; also use static asserts to check sizes---*/
static GMBits1x64 GMBits1x64_null() {
  GMStaticAssert(sizeof(GMBits1) * 8 >= 1);
  GMStaticAssert(sizeof(GMBits1x64) == 8);
  GMStaticAssert(sizeof(GMUInt64) == 8);
  GMStaticAssert(sizeof(GMUInt64) == sizeof(GMBits1x64)); /*---for Magma---*/

  GMBits1x64 value;
  value = 0;
  return value;
}

// TODO: remove the following.

typedef double GMBits;

typedef unsigned long long int GMULInt;

/*===========================================================================*/
/*---Types for Czekanowski metric---*/

//---TODO: revise nomenclature to be different from "GMFloat2" ...

#ifdef FP_PRECISION_SINGLE
  typedef float GMFloat;
  enum { GM_MPI_FLOAT = MPI_FLOAT };
#ifdef FP_PRECISION_DOUBLE
#error Cannot set both FP_PRECISION_SINGLE and FP_PRECISION_DOUBLE.
#endif
#else
#ifdef FP_PRECISION_DOUBLE
  typedef double GMFloat;
  enum { GM_MPI_FLOAT = MPI_DOUBLE };
#else
#error Must set FP_PRECISION_SINGLE or FP_PRECISION_DOUBLE.
#endif
#endif

/*===========================================================================*/
/*---Types for CCC metric---*/

/*---For Vectors: single 2-bit value (semi-nibble):
     use unsigned int as a container for a single item---*/
typedef unsigned int GMBits2;

/*---For Vectors: packed: 2 long integers, used to store 64 semi-nibbles---*/
typedef unsigned long long int GMBits1_2x64;
typedef struct { GMBits1_2x64 data[2]; } GMBits2x64;

/*---For Vectors: largest allowed size of a data value---*/
enum { GM_BITS2_MAX_VALUE_BITS = 2 };

/*---For metrics: single integer to store a tally result---*/
typedef unsigned int GMTally1;

/*---For Metrics: 2 (4) doubles to store 4 (8) packed tally results:
     use 25 bits of each 52-bit mantissa to store a result---*/
typedef struct { GMFp64 data[2]; } GMTally2x2;
typedef struct { GMFp64 data[4]; } GMTally4x2;

/*---For Metrics: largest allowed size of a data value---*/
enum { GM_TALLY1_MAX_VALUE_BITS = 25 };

/*---For Metrics: for packing of multipliers---*/
typedef GMFp64 GMFloat2;
typedef struct { GMFp64 data[2]; } GMFloat3;

/*===========================================================================*/
/*---Types for CCC metric: functions---*/

/*---Return null value; also use static asserts to check sizes---*/
static GMBits2x64 GMBits2x64_null() {
  GMStaticAssert(sizeof(GMBits2) * 8 >= GM_BITS2_MAX_VALUE_BITS);
  GMStaticAssert(sizeof(GMBits1_2x64) == 8);
  GMStaticAssert(sizeof(GMBits2x64) == 2 * sizeof(GMBits1_2x64));
  GMStaticAssert(sizeof(GMBits2x64) == 16);
  GMStaticAssert(sizeof(GMTally2x2) == sizeof(GMBits2x64)); /*---for Magma---*/

  GMBits2x64 value;
  value.data[0] = 0;
  value.data[1] = 0;
  return value;
}

/*---Return null value; also use static asserts to check sizes---*/
static GMTally2x2 GMTally2x2_null() {
  GMStaticAssert(sizeof(GMTally1) * 8 >= GM_TALLY1_MAX_VALUE_BITS);
  GMStaticAssert(sizeof(GMTally2x2) == 16);
  GMStaticAssert(sizeof(GMTally2x2) == sizeof(GMBits2x64)); /*---for Magma---*/

  GMTally2x2 value;
  value.data[0] = 0;
  value.data[1] = 0;
  return value;
}

/*---Return null value; also use static asserts to check sizes---*/
static GMTally4x2 GMTally4x2_null() {
  GMStaticAssert(sizeof(GMTally4x2) == 32);

  GMTally4x2 value;
  value.data[0] = 0;
  value.data[1] = 0;
  value.data[2] = 0;
  value.data[3] = 0;
  return value;
}

/*---Encode/decode between float and pair of tally values---*/

static void GMTally1_decode(GMTally1* __restrict__ val0,
                            GMTally1* __restrict__ val1,
                            GMFloat v) {
  GMAssert(val0 != NULL);
  GMAssert(val1 != NULL);
  const GMUInt64 tally2 = (GMUInt64)v;
  GMAssert(v == (GMFloat)tally2);
  const GMTally1 v0 =
      tally2 & ((((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS) - 1);
  const GMTally1 v1 = tally2 >> GM_TALLY1_MAX_VALUE_BITS;
  *val0 = v0;
  *val1 = v1;
  GMAssert(v ==
           (GMFloat)(v0 + (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS) * v1));
  GMAssert(v0 >= 0);
  GMAssert(v1 >= 0);
  GMAssert(v0 < (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS));
  GMAssert(v1 < (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS));
}

/*-----*/

static GMFloat GMTally1_encode(GMTally1 val0, GMTally1 val1) {
  const GMUInt64 tally2 =
      val0 + (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS) * val1;
  const GMFloat result = (GMFloat)tally2;
  GMAssert(val0 == (((GMUInt64)result) &
                    ((((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS) - 1)));
  GMAssert(val1 == ((GMUInt64)result) >> GM_TALLY1_MAX_VALUE_BITS);
  return result;
}

/*---Encode/decode for multipliers/sums---*/

static void GMFloat2_decode(GMTally1* __restrict__ val0,
                            GMTally1* __restrict__ val1,
                            GMFloat2 v) {
  GMTally1_decode(val0, val1, v);
}

/*-----*/

static void GMFloat3_decode(GMTally1* __restrict__ val0,
                            GMTally1* __restrict__ val1,
                            GMTally1* __restrict__ val2,
                            GMFloat3 v) {
  GMTally1_decode(val0, val1, v.data[0]);
  GMTally1 dummy;
  GMTally1_decode(val2, &dummy, v.data[1]);
}

/*-----*/

static GMFloat2 GMFloat2_encode(GMTally1 val0, GMTally1 val1) {
  return GMTally1_encode(val0, val1);
}

/*-----*/

static GMFloat3 GMFloat3_encode(GMTally1 val0, GMTally1 val1, GMTally1 val2) {
  GMFloat3 result; /*---here we should set = null to be super cautious---*/
  const GMTally1 dummy = 0;
  result.data[0] = GMTally1_encode(val0, val1);
  result.data[1] = GMTally1_encode(val2, dummy);
  return result;
}

/*---Get an entry---*/

static GMTally1 GMTally2x2_get(GMTally2x2 tally2x2, int i0, int i1) {
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);

  const GMUInt64 tally2 = tally2x2.data[i0];

  const GMTally1 result =
      i1 == 0 ? tally2 % (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS)
              : tally2 / (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS);
  GMAssert(result >= 0);
  GMAssert(result < (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS));
  return result;
}

/*---Get an entry---*/

static GMTally1 GMTally4x2_get(GMTally4x2 tally4x2, int i0, int i1, int i2) {
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);
  GMAssert(i2 >= 0 && i2 < 2);

  const GMUInt64 tally2 = tally4x2.data[i1 + 2 * i0];

  const GMTally1 result =
      i2 == 0 ? tally2 % (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS)
              : tally2 / (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS);
  GMAssert(result >= 0);
  GMAssert(result < (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS));
  return result;
}

/*===========================================================================*/
/*---Type ids---*/

enum {
  GM_DATA_TYPE_FLOAT = 1,
  GM_DATA_TYPE_BITS1 = 2,
  GM_DATA_TYPE_UINT64 = 3,  //---(design of this entry is not complete)
  GM_DATA_TYPE_BITS2 = 4,
  GM_DATA_TYPE_TALLY2X2 = 5,
  GM_DATA_TYPE_TALLY4X2 = 6
};

/*===========================================================================*/
/*---Dual CPU/GPU pointer---*/

typedef struct {
  void* __restrict__ h;
  void* __restrict__ d;
} GMMirroredPointer;

GMMirroredPointer GMMirroredPointer_null(void);

/*===========================================================================*/
/*---Checksums---*/

/*---Struct with checksum info---*/

enum { GM_CHECKSUM_SIZE = 3 };

typedef struct { size_t data[GM_CHECKSUM_SIZE]; } GMChecksum;

/*---------------------------------------------------------------------------*/

static GMChecksum GMChecksum_null() {
  GMChecksum result;
  int i = 0;
  for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    result.data[i] = 0;
  }
  return result;
}

/*---------------------------------------------------------------------------*/

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
  /*---Counters---*/
  double time;
  double ops;
  /*---MPI---*/
  int mpi_comm_;
  int mpi_comm_vector_;
  int mpi_comm_field_;
  /*---*/
  int num_proc_world_;
  int num_proc_;
  int num_proc_field_;
  int num_proc_repl_;
  int num_proc_vector_i_;
  int num_proc_vector_j_;
  int num_proc_vector_k_;
  int num_proc_vector_;
  /*---*/
  int proc_num_;
  int proc_num_field_;
  int proc_num_vector_i_;
  int proc_num_vector_j_;
  int proc_num_vector_k_;
  int proc_num_vector_;
  /*---*/
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

enum { GM_NUM_WAY_2 = 2, GM_NUM_WAY_3 = 3, GM_NUM_NUM_WAY = 2 };

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
/*---Accessors: general---*/

//---TODO: fix Env_* functions to follow namespace guidelines.

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

static int Env_mpi_comm_vector(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->mpi_comm_vector_;
}

/*---------------------------------------------------------------------------*/

static int Env_mpi_comm_field(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->mpi_comm_field_;
}

/*---------------------------------------------------------------------------*/

static int Env_is_proc_active(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->is_proc_active_;
}

/*---------------------------------------------------------------------------*/

void Env_set_compute_method(GMEnv* env, int compute_method);
int Env_data_type_vectors(const GMEnv* env);
int Env_data_type_metrics(const GMEnv* env);

void Env_set_num_proc(GMEnv* env, int num_proc_vector, int num_proc_repl,
                      int num_proc_field);

cudaStream_t Env_stream_compute(GMEnv* env);
cudaStream_t Env_stream_togpu(GMEnv* env);
cudaStream_t Env_stream_fromgpu(GMEnv* env);

/*===========================================================================*/
/*---Accessors: num proc---*/

static int Env_num_block_vector(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->num_proc_vector_i_;
}

/*---------------------------------------------------------------------------*/

static int Env_num_proc_vector_i(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->num_proc_vector_i_;
}

/*---------------------------------------------------------------------------*/

static int Env_num_proc_vector_j(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->num_proc_vector_j_;
}

/*---------------------------------------------------------------------------*/

static int Env_num_proc_vector_k(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->num_proc_vector_k_;
}

/*---------------------------------------------------------------------------*/

static int Env_num_proc_vector(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->num_proc_vector_;
}

/*---------------------------------------------------------------------------*/

static int Env_num_proc_field(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->num_proc_field_;
}

/*===========================================================================*/
/*---Accessors: proc_num---*/

static int Env_proc_num(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->proc_num_;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_num_vector_i(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->proc_num_vector_i_;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_num_vector_j(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->proc_num_vector_j_;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_num_vector_k(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->proc_num_vector_k_;
}

/*---------------------------------------------------------------------------*/

static int Env_proc_num_field(const GMEnv* env) {
  GMAssert(env != NULL);
  return env->proc_num_field_;
}

/*===========================================================================*/
/*---Timer functions---*/

double GMEnv_get_time(GMEnv* env);
double GMEnv_get_synced_time(GMEnv* env);

/*===========================================================================*/
/*---Math utility functions---*/

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

  for (result = 0, n--; result <= 8 * (int)sizeof(size_t); ++result) {
    if (n == 0) {
      break;
    }
    n >>= 1;
  }

  return result;
}

/*---------------------------------------------------------------------------*/

static size_t gm_nchoosek(int n, int k) {
  GMAssert(n >= 0);
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

/*---------------------------------------------------------------------------*/

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
