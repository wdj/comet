/*---------------------------------------------------------------------------*/
/*!
 * \file   env.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Environment settings and general utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================
=============================================================================*/

#ifndef _gm_env_hh_
#define _gm_env_hh_

#include <stddef.h>
#include <assert.h>
#include <float.h>
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

void gm_assert(const char* condition_string, const char* file, int line);

void gm_insist(void const * const env,
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


typedef unsigned int GMUInt32;

typedef signed long long int GMInt64;
typedef unsigned long long int GMUInt64;

/*---Floating point of explicit (double) precision---*/

typedef float GMFp32;
typedef double GMFp64;

static void gm_check_type_sizes() {
  GMStaticAssert(sizeof(GMFp32) == 32/8);
  GMStaticAssert(sizeof(GMFp64) == 64/8);
  GMStaticAssert(sizeof(int) == 4);
  GMStaticAssert(sizeof(size_t) == 8);
  GMStaticAssert(sizeof(GMUInt32) == 4);
  GMStaticAssert(sizeof(GMInt64) == 8);
  GMStaticAssert(sizeof(GMUInt64) == 8);
}

/*===========================================================================*/
/*---Types for Sorenson metric---*/

//---TODO: complete the design for Sorenson metric.

/*---For Vectors: single 1-bit value: use unsigned int---*/
typedef unsigned int GMBits1;

typedef unsigned long long int GMBits1x64;

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

/*===========================================================================*/
/*---Types for Czekanowski metric---*/

//---TODO: revise nomenclature to be different from "GMFloat2" ...

#ifdef FP_PRECISION_SINGLE
  typedef float GMFloat;
  enum { GM_MPI_FLOAT = MPI_FLOAT };
  enum { GM_FP_PRECISION_DOUBLE = GM_BOOL_FALSE };
#ifdef FP_PRECISION_DOUBLE
#error Cannot set both FP_PRECISION_SINGLE and FP_PRECISION_DOUBLE.
#endif
#else
#ifdef FP_PRECISION_DOUBLE
  typedef double GMFloat;
  enum { GM_MPI_FLOAT = MPI_DOUBLE };
  enum { GM_FP_PRECISION_DOUBLE = GM_BOOL_TRUE };
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
                            GMFp64 v) {
  GMAssert(val0 != NULL);
  GMAssert(val1 != NULL);
  const GMUInt64 tally2 = (GMUInt64)v;
  GMAssert(v == (GMFp64)tally2);
  const GMTally1 v0 =
      tally2 & ((((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS) - 1);
  const GMTally1 v1 = tally2 >> GM_TALLY1_MAX_VALUE_BITS;
  *val0 = v0;
  *val1 = v1;
  GMAssert(v ==
           (GMFp64)(v0 + (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS) * v1));
  GMAssert(v0 >= 0);
  GMAssert(v1 >= 0);
  GMAssert(v0 < (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS));
  GMAssert(v1 < (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS));
}

/*-----*/

static GMFp64 GMTally1_encode(GMTally1 val0, GMTally1 val1) {
  const GMUInt64 tally2 =
      val0 + (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS) * val1;
  const GMFp64 result = (GMFp64)tally2;
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
  size_t size;
} GMMirroredPointer;

GMMirroredPointer GMMirroredPointer_null(void);

/*===========================================================================*/
/*---Checksums---*/

/*---Multiprecision integers---*/

enum { GM_MULTIPREC_INT_SIZE = 16 };

typedef struct {
  size_t data[GM_MULTIPREC_INT_SIZE];
} GMMultiprecInt;

/*---Struct with checksum info---*/

enum { GM_CHECKSUM_SIZE = 3 };

typedef struct {
  size_t data[GM_CHECKSUM_SIZE];
  _Bool is_overflowed;
  double value_max;
  GMMultiprecInt sum;
  double sum_d;
  _Bool is_started;
} GMChecksum;

/*---------------------------------------------------------------------------*/

static GMChecksum GMChecksum_null() {
  GMChecksum result;
  int i = 0;
  for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    result.data[i] = 0;
  }
  result.is_overflowed = GM_BOOL_FALSE;
  result.value_max = -DBL_MAX;
  GMMultiprecInt sum = {0};
  result.sum = sum;
  result.sum_d = 0;
  return result;
}

/*---------------------------------------------------------------------------*/

static _Bool gm_are_checksums_equal(GMChecksum c1, GMChecksum c2) {
  _Bool result = GM_BOOL_TRUE;
  int i = 0;
  for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    result = result && c1.data[i] == c2.data[i];
  }
  result = result && c1.is_overflowed == c2.is_overflowed;
  result = result && c1.value_max == c2.value_max;
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
  int num_stage;
  int stage_num;
  int num_phase;
  int phase_num;
  double ccc_param_;
  /*---Counters---*/
  double time;
  double compares;
  double ops_local;
  double ops;
  size_t cpu_mem;
  size_t cpu_mem_max;
  size_t gpu_mem;
  size_t gpu_mem_max;
  /*---MPI---*/
  _Bool make_comms_;
  int mpi_comm_base_;
  int mpi_comm_;
  int mpi_comm_vector_;
  int mpi_comm_field_;
  /*---*/
  int num_proc_base_;
  int num_proc_;
  int num_proc_field_;
  int num_proc_repl_;
  int num_proc_vector_i_;
  int num_proc_vector_total_;
  _Bool are_mpi_comms_initialized_;
  /*---*/
  int proc_num_base_;
  int proc_num_;
  int proc_num_field_;
  int proc_num_repl_;
  int proc_num_vector_i_;
  int proc_num_vector_;
  /*---*/
  _Bool is_proc_active_;
  /*---CUDA---*/
  cudaStream_t stream_compute_;
  cudaStream_t stream_togpu_;
  cudaStream_t stream_fromgpu_;
  _Bool are_cuda_streams_initialized_;
  /*---OTHER---*/
  const char* description;
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

void gm_create_args(char* argstring, int* argc, char** argv);

void GMEnv_create(GMEnv* const env, MPI_Comm comm, int argc, char** argv,
                  const char* const description);

void GMEnv_create(GMEnv* const env, MPI_Comm comm, const char* const options,
                  const char* const description);

void GMEnv_create_no_comms(GMEnv* const env, int argc, char** argv,
                           const char* const description,
                           int num_proc, int proc_num);

void GMEnv_create_no_comms(GMEnv* const env, const char* const options,
                           const char* const description,
                           int num_proc, int proc_num);

/*===========================================================================*/
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* const env);

/*===========================================================================*/
/*---Manage cuda streams---*/

void GMEnv_initialize_streams(GMEnv* const env);
void GMEnv_terminate_streams(GMEnv* const env);

/*===========================================================================*/
/*---Accessors: general---*/

static int GMEnv_metric_type(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->metric_type_;
}

/*---------------------------------------------------------------------------*/

static double GMEnv_ccc_param(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->ccc_param_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_num_way(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->num_way_;
}

/*---------------------------------------------------------------------------*/

static _Bool GMEnv_all2all(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->all2all_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_compute_method(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->compute_method_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_mpi_comm(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->mpi_comm_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_mpi_comm_vector(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->mpi_comm_vector_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_mpi_comm_field(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->mpi_comm_field_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_is_proc_active(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->is_proc_active_;
}

/*---------------------------------------------------------------------------*/

void GMEnv_set_compute_method(GMEnv* const env, int compute_method);
int GMEnv_data_type_vectors(const GMEnv* const env);
int GMEnv_data_type_metrics(const GMEnv* const env);

void GMEnv_set_num_proc(GMEnv* const env, int num_proc_vector_i,
                      int num_proc_repl, int num_proc_field);

cudaStream_t GMEnv_stream_compute(GMEnv* const env);
cudaStream_t GMEnv_stream_togpu(GMEnv* const env);
cudaStream_t GMEnv_stream_fromgpu(GMEnv* const env);

/*===========================================================================*/
/*---Accessors: num proc---*/

static int GMEnv_num_block_vector(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->num_proc_vector_i_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_num_proc_vector_i(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->num_proc_vector_i_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_num_proc_repl(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->num_proc_repl_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_num_proc_vector_total(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->num_proc_vector_total_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_num_proc_field(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->num_proc_field_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_num_proc(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->num_proc_;
}

/*===========================================================================*/
/*---Accessors: proc_num---*/

static int GMEnv_proc_num(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->proc_num_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_proc_num_vector_i(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->proc_num_vector_i_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_proc_num_repl(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->proc_num_repl_;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_proc_num_field(const GMEnv* const env) {
  GMAssert(env != NULL);
  return env->proc_num_field_;
}

/*===========================================================================*/
/*---Timer functions---*/

double GMEnv_get_time(const GMEnv* const env);
double GMEnv_get_synced_time(const GMEnv* const env);

/*===========================================================================*/
/*---Math utility functions---*/

static int gm_min_i(const int i, const int j) {
  return i < j ? i : j;
}

/*---------------------------------------------------------------------------*/

static size_t gm_min_i8(const size_t i, const size_t j) {
  return i < j ? i : j;
}

/*---------------------------------------------------------------------------*/

static int gm_max_i(const int i, const int j) {
  return i > j ? i : j;
}

/*---------------------------------------------------------------------------*/

static size_t gm_max_i8(const size_t i, const size_t j) {
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

static int gm_mod_i(const int i, const int j) {
  GMAssert(j > 0);

  return i - j * gm_floor_i(i, j);
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

/*---------------------------------------------------------------------------*/

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

/*===========================================================================*/
/*---Misc.---*/

_Bool GMEnv_cuda_last_call_succeeded(const GMEnv* const env);

/*---------------------------------------------------------------------------*/

void* gm_malloc(size_t n, GMEnv* env);
void gm_free(void* p, size_t n, GMEnv* env);

GMFloat* GMFloat_malloc(size_t n, GMEnv* env);
void GMFloat_free(GMFloat* p, size_t n, GMEnv* env);

void GMFloat_fill_nan(GMFloat* const a, size_t n);
void GMFloat_check(GMFloat* const a, size_t n);

/*---------------------------------------------------------------------------*/
/*---NOTE: the following does not specialize based on part1/2/3---*/

static int GMEnv_num_section_steps(const GMEnv* const env, int part_num) {
  GMAssert(env != NULL);
  GMAssert(part_num >= 1 && part_num <= 3);
  return (!GMEnv_all2all(env)) || GMEnv_num_proc_repl(env) == 1 || part_num == 3 ? 1 : 6;
  //return part_num == 3 ? 1 : 1;
}

/*---------------------------------------------------------------------------*/

static int GMEnv_num_sections(const GMEnv* const env, int part_num) {
  GMAssert(env != NULL);
  GMAssert(part_num >= 1 && part_num <= 3);
  return part_num == 3 ? 6 : GMEnv_num_section_steps(env, part_num);
}

/*---------------------------------------------------------------------------*/

int gm_mpi_type(const GMEnv* const env);

/*===========================================================================*/

//#ifdef __cplusplus
//extern "C" {
//#endif
//
//#ifdef __cplusplus
//} /*---extern "C"---*/
//#endif

/*===========================================================================*/

#endif /*---_gm_env_hh_---*/

/*---------------------------------------------------------------------------*/
